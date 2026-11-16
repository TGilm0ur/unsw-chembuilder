"""
chem_engine.py - all the chemistry behind UNSW ChemBuilder.

Everything here is deliberately SIMPLE (teaching-level) and uses only RDKit + numpy:

  * geometry optimisation : MMFF94 molecular mechanics (bond-stretch/bend/torsion springs
                            + simple partial-charge electrostatics). We record every few
                            optimiser steps so the front end can animate the molecule
                            "relaxing" into place.
  * partial charges       : Gasteiger-Marsili charges (electronegativity equalisation).
  * dipole moment         : mu = sum(q_i * r_i)   (point charges on the optimised geometry)
  * IR spectrum           : each bond is a harmonic spring, nu = (1/2 pi c) sqrt(k/mu_reduced)
                            with typical force constants; intensity ~ bond polarity.
  * thermodynamics        : average bond enthalpies (Hess's Law) for dHf and dHc, and a
                            statistical-mechanics estimate of entropy from mass, shape,
                            symmetry and vibrations; dG = dH - T dS.

None of these are research-grade numbers. They are meant to show the right trends.
"""
from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem import rdForceFieldHelpers as FFH
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

RDLogger.DisableLog("rdApp.*")

ALLOWED = ("C", "H", "N", "O")
MAX_HEAVY = 4
VALENCE = {"H": 1, "C": 4, "N": 3, "O": 2}
DEBYE_PER_EANG = 4.80320471  # 1 e*Angstrom in debye

BOND_TYPES = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}


# --------------------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------------------
def validate(graph: dict) -> tuple[list[str], list[str]]:
    """Return (errors, warnings). Errors block the calculation, warnings do not."""
    atoms = graph.get("atoms", []) or []
    bonds = graph.get("bonds", []) or []
    errors: list[str] = []
    warnings: list[str] = []
    if not atoms:
        return ["Add some atoms to start building."], warnings

    ids = {a["id"]: a for a in atoms}
    for a in atoms:
        if a["el"] not in ALLOWED:
            errors.append(f"Element {a['el']} is not supported (use C, H, N, O).")
    heavy = [a for a in atoms if a["el"] != "H"]
    if len(heavy) > MAX_HEAVY:
        errors.append(f"Too many non-hydrogen atoms ({len(heavy)}). The maximum is {MAX_HEAVY}.")

    used = {i: 0 for i in ids}
    seen = set()
    for b in bonds:
        a, c, o = b["a"], b["b"], int(b["o"])
        if a not in ids or c not in ids or a == c:
            errors.append("A bond is connected to a missing atom.")
            continue
        key = (min(a, c), max(a, c))
        if key in seen:
            errors.append("Two bonds join the same pair of atoms.")
        seen.add(key)
        used[a] += o
        used[c] += o
    for i, u in used.items():
        el = ids[i]["el"]
        if el in VALENCE and u > VALENCE[el]:
            errors.append(f"A {el} atom has too many bonds ({u}; max {VALENCE[el]}).")
    if len(atoms) == 1 and atoms[0]["el"] == "H":
        errors.append("A single hydrogen atom is not a molecule - add a C, N or O.")
    if errors:
        return errors, warnings

    # ---- friendly warnings (calculation still runs) ----
    parent = {i: i for i in ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for b in bonds:
        parent[find(b["a"])] = find(b["b"])
    n_frag = len({find(i) for i in ids})
    if n_frag > 1:
        warnings.append(
            f"Your drawing has {n_frag} separate pieces. They are not bonded, so "
            "the results will treat them as one cluster - connect them for a single molecule."
        )
    for b in bonds:
        if ids[b["a"]]["el"] == "O" and ids[b["b"]]["el"] == "O" and b["o"] == 1:
            warnings.append("O-O single bonds (peroxides) are unusual and weak - results are rough.")
            break
    # ring size via simple DFS on heavy+H graph (rings of 3-4 are strained)
    adj: dict[int, list[int]] = {i: [] for i in ids}
    for b in bonds:
        adj[b["a"]].append(b["b"])
        adj[b["b"]].append(b["a"])
    if _has_small_ring(adj, ids):
        warnings.append("Your molecule has a 3- or 4-membered ring. These are strained, so expect distorted geometry.")
    return errors, warnings


def _has_small_ring(adj, ids, max_size=4) -> bool:
    for start in ids:
        stack = [(start, [start])]
        while stack:
            node, path = stack.pop()
            if len(path) > max_size:
                continue
            for nb in adj[node]:
                if nb == start and len(path) >= 3:
                    return True
                if nb not in path and len(path) < max_size:
                    stack.append((nb, path + [nb]))
    return False


def hill_formula(counts: dict) -> str:
    """C first then H then alphabetical (Hill); without carbon: N, H, O (gives NH3, H2O)."""
    order = ["C", "H", "N", "O"] if counts.get("C") else ["N", "H", "O"]
    return "".join(f"{e}{'' if counts[e] == 1 else counts[e]}" for e in order if counts.get(e))


def topology_key(graph: dict) -> str:
    """Hash of atoms + bonds only (not the 2D drawing positions)."""
    atoms = [(a["id"], a["el"]) for a in graph.get("atoms", [])]
    bonds = [(b["a"], b["b"], int(b["o"])) for b in graph.get("bonds", [])]
    return hashlib.md5(json.dumps([atoms, bonds]).encode()).hexdigest()[:12]


# --------------------------------------------------------------------------------------
# Building the RDKit molecule
# --------------------------------------------------------------------------------------
def build_mol(graph: dict):
    """Return (mol_with_Hs, order_map). Atom i of the graph keeps index i."""
    rw = Chem.RWMol()
    idx = {}
    for a in graph["atoms"]:
        idx[a["id"]] = rw.AddAtom(Chem.Atom(a["el"]))
    order_map = {}
    for b in graph["bonds"]:
        i, j, o = idx[b["a"]], idx[b["b"]], int(b["o"])
        rw.AddBond(i, j, BOND_TYPES[o])
        order_map[(min(i, j), max(i, j))] = o
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    molH = Chem.AddHs(mol)
    return molH, order_map


def _embed(molH, seed: int):
    ps = AllChem.ETKDGv3()
    ps.randomSeed = seed
    cid = AllChem.EmbedMolecule(molH, ps)
    if cid < 0:
        ps.useRandomCoords = True
        cid = AllChem.EmbedMolecule(molH, ps)
    if cid < 0:
        raise RuntimeError("Could not generate a 3D structure for this molecule.")
    return cid


# typical bond lengths (Angstrom) used to restrain the generic fallback force field
IDEAL_LENGTH = {
    (("C", "H"), 1): 1.09, (("H", "N"), 1): 1.01, (("H", "O"), 1): 0.96, (("H", "H"), 1): 0.74,
    (("C", "C"), 1): 1.54, (("C", "C"), 2): 1.34, (("C", "C"), 3): 1.20,
    (("C", "N"), 1): 1.47, (("C", "N"), 2): 1.28, (("C", "N"), 3): 1.16,
    (("C", "O"), 1): 1.43, (("C", "O"), 2): 1.21, (("C", "O"), 3): 1.13,
    (("N", "N"), 1): 1.45, (("N", "N"), 2): 1.25, (("N", "N"), 3): 1.10,
    (("N", "O"), 1): 1.40, (("N", "O"), 2): 1.20,
    (("O", "O"), 1): 1.48, (("O", "O"), 2): 1.21,
}


def _mmff_unreliable(molH, order_map) -> bool:
    """MMFF gives poor bond lengths for CO2-like cumulated systems and bonds between N/O atoms."""
    om = order_map or {}
    for a in molH.GetAtoms():
        i = a.GetIdx()
        n_double = sum(1 for nb in a.GetNeighbors() if om.get((min(i, nb.GetIdx()), max(i, nb.GetIdx())), 1) == 2)
        if n_double >= 2:
            return True
    for b in molH.GetBonds():
        e = {b.GetBeginAtom().GetSymbol(), b.GetEndAtom().GetSymbol()}
        if e <= {"N", "O"}:
            return True
    return False


def _forcefield(molH, cid=0, order_map=None):
    """MMFF94 when it knows every atom type; otherwise UFF restrained to typical bond lengths."""
    if FFH.MMFFHasAllMoleculeParams(molH) and not _mmff_unreliable(molH, order_map):
        props = FFH.MMFFGetMoleculeProperties(molH)
        return AllChem.MMFFGetMoleculeForceField(molH, props, confId=cid)
    ff = AllChem.UFFGetMoleculeForceField(molH, confId=cid)
    for b in molH.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        order = (order_map or {}).get((min(i, j), max(i, j)), 1)
        pair = tuple(sorted((molH.GetAtomWithIdx(i).GetSymbol(), molH.GetAtomWithIdx(j).GetSymbol())))
        L = IDEAL_LENGTH.get((pair, order))
        if pair == ("C", "O") and order == 2:
            for k in (i, j):
                at = molH.GetAtomWithIdx(k)
                if at.GetSymbol() == "C" and sum(1 for nb in at.GetNeighbors() if nb.GetSymbol() == "O"
                        and (order_map or {}).get((min(k, nb.GetIdx()), max(k, nb.GetIdx())), 1) == 2) >= 2:
                    L = 1.16
        if L:
            ff.UFFAddDistanceConstraint(i, j, False, L, L, 20000.0)
    return ff


def _kabsch_align(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Rotate/translate P onto Q (both Nx3). Returns aligned P."""
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return (R @ Pc.T).T + Q.mean(0)


def trajectory(graph: dict, max_frames: int = 40) -> dict:
    """
    Build the 3D structure and the optimisation movie.
    Returns {elements, bonds:[[i,j,order]], frames:[Nx3 lists], ...}; frame 0 is the rough
    start, the last frame is the optimised geometry.
    """
    molH, order_map = build_mol(graph)
    seed = int(topology_key(graph), 16) % 2**31
    n = molH.GetNumAtoms()
    cid = _embed(molH, seed)
    base = np.array(molH.GetConformer(cid).GetPositions())

    # Rough start: well-formed topology, but stretched/squashed and jittered
    rng = np.random.default_rng(seed)
    start = base + rng.normal(0.0, 0.32, base.shape)
    if n > 2:
        start = start * rng.uniform(0.9, 1.15, (1, 3))
    conf = molH.GetConformer(cid)
    for i in range(n):
        conf.SetAtomPosition(i, start[i].tolist())

    frames = [start.copy()]
    if n > 1:
        ff = _forcefield(molH, cid, order_map)
        ff.Initialize()
        total = 0
        while total < 1500:
            res = ff.Minimize(maxIts=3)
            total += 3
            frames.append(np.array(conf.GetPositions()))
            if res == 0:
                break
        ff.Minimize(maxIts=2000)  # polish
        frames.append(np.array(conf.GetPositions()))

    final = frames[-1]
    if len(frames) > max_frames:
        pick = np.linspace(0, len(frames) - 1, max_frames).round().astype(int)
        frames = [frames[k] for k in pick]
    aligned = []
    for f in frames:
        f = f if len(f) == 1 else _kabsch_align(f, final)
        aligned.append((f - final.mean(0)).round(4).tolist())

    bonds = []
    for b in molH.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bonds.append([i, j, order_map.get((min(i, j), max(i, j)), 1)])
    return {
        "elements": [a.GetSymbol() for a in molH.GetAtoms()],
        "bonds": bonds,
        "frames": aligned,
        "key": topology_key(graph),
    }


# --------------------------------------------------------------------------------------
# Partial charges, dipole
# --------------------------------------------------------------------------------------
def partial_charges(molH) -> list[float]:
    """Gasteiger-Marsili charges: electronegativity equalisation along the bonds (no quantum mechanics).
    MMFF charges are not used because MMFF sets every C-H bond to exactly zero."""
    AllChem.ComputeGasteigerCharges(molH)
    return [float(a.GetProp("_GasteigerCharge")) for a in molH.GetAtoms()]


# --------------------------------------------------------------------------------------
# IR: harmonic springs
# --------------------------------------------------------------------------------------
# typical stretching force constants, N/m  (key: sorted element pair, order)
FORCE_K = {
    (("C", "H"), 1): 490, (("H", "N"), 1): 650, (("H", "O"), 1): 780, (("H", "H"), 1): 575,
    (("C", "C"), 1): 450, (("C", "C"), 2): 960, (("C", "C"), 3): 1560,
    (("C", "N"), 1): 500, (("C", "N"), 2): 1000, (("C", "N"), 3): 1800,
    (("C", "O"), 1): 490, (("C", "O"), 2): 1200, (("C", "O"), 3): 1900,
    (("N", "N"), 1): 400, (("N", "N"), 2): 1300, (("N", "N"), 3): 2300,
    (("N", "O"), 1): 450, (("N", "O"), 2): 1400,
    (("O", "O"), 1): 350, (("O", "O"), 2): 1180,
}
AROMATIC_K = 800
MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999}
C_CM = 2.99792458e10  # cm/s
AMU = 1.66053907e-27


def _wavenumber(k: float, m1: float, m2: float) -> float:
    mu = (m1 * m2 / (m1 + m2)) * AMU
    return math.sqrt(k / mu) / (2 * math.pi * C_CM)


def _force_k(e1, e2, order, aromatic):
    if aromatic:
        return AROMATIC_K
    key = (tuple(sorted((e1, e2))), order)
    if key in FORCE_K:
        return FORCE_K[key]
    base = FORCE_K.get((tuple(sorted((e1, e2))), 1), 450)
    return base * (1 + 0.9 * (order - 1) * 1.0) * (1 if order == 1 else 1.0)


def stretch_modes(molH, bond_info, charges):
    """One entry per bond: wavenumber, intensity, atoms."""
    modes = []
    for b in bond_info:
        i, j, order, arom = b["a"], b["b"], b["order"], b["aromatic"]
        e1, e2 = b["els"]
        k = _force_k(e1, e2, order, arom)
        nu = _wavenumber(k, MASS[e1], MASS[e2])
        inten = 0.08 + abs(charges[i] - charges[j])
        modes.append({"i": i, "j": j, "nu": nu, "I": inten, "els": (e1, e2),
                      "order": order, "arom": arom})
    return modes


def ir_spectrum(molH, bond_info, charges) -> dict:
    # A homonuclear diatomic (H2, N2, O2) has no dipole at all, and stretching it does not create one,
    # so it cannot absorb infrared light.
    atoms_el = [a.GetSymbol() for a in molH.GetAtoms()]
    if len(atoms_el) == 2 and atoms_el[0] == atoms_el[1]:
        x = np.arange(4000, 399, -4.0)
        return {"peaks": [], "inactive": True, "x": x.astype(int).tolist(), "t": [100.0] * len(x)}
    modes = stretch_modes(molH, bond_info, charges)
    sym = {1: "-", 2: "=", 3: "#"}
    # group X-H bonds on the same heavy atom -> symmetric/asymmetric pair
    groups: dict[tuple, list] = {}
    peaks_raw = []
    for m in modes:
        e1, e2 = m["els"]
        if "H" in (e1, e2) and e1 != e2:
            heavy = m["i"] if e1 != "H" else m["j"]
            groups.setdefault((heavy, "XH"), []).append(m)
        else:
            peaks_raw.append(([m], m["nu"], m["I"], None))
    for (heavy, _), ms in groups.items():
        if len(ms) >= 2:
            nu = float(np.mean([m["nu"] for m in ms]))
            tot = sum(m["I"] for m in ms)
            peaks_raw.append((ms, nu - 45, tot * 0.45, "symmetric"))
            peaks_raw.append((ms, nu + 55, tot * 0.55, "asymmetric"))
        else:
            peaks_raw.append((ms, ms[0]["nu"], ms[0]["I"], None))

    def label(ms, paired, k):
        a, b = sorted(ms[0]["els"], key=lambda e: "CNOH".index(e))
        return f"{a}{'~' if ms[0]['arom'] else sym[ms[0]['order']]}{b}"

    merged: dict[tuple, dict] = {}
    for ms, nu, inten, paired in peaks_raw:
        lab = label(ms, paired, 0)
        key = (lab, round(nu / 4), paired)
        entry = merged.setdefault(key, {"wn": nu, "I": 0.0, "label": lab, "bonds": [], "groups": [], "mode": paired})
        entry["I"] += inten
        # one group = the X-H bonds on the same heavy atom (they vibrate together)
        entry["groups"].append([[m["i"], m["j"]] for m in ms])
        for m in ms:
            if [m["i"], m["j"]] not in entry["bonds"]:
                entry["bonds"].append([m["i"], m["j"]])
    peaks = sorted(merged.values(), key=lambda p: -p["wn"])
    for p in peaks:
        p["I"] = round(min(1.0, p["I"]), 3)
        p["wn"] = round(p["wn"])
        width = 70 if p["label"] == "O-H" else 35 if p["label"] == "N-H" else 14
        p["width"] = width

    x = np.arange(4000, 399, -4.0)
    A = np.zeros_like(x)
    for p in peaks:
        w = p["width"]
        A += p["I"] * (w * w) / ((x - p["wn"]) ** 2 + w * w)
    T = 100 * np.exp(-1.6 * A)
    return {"peaks": peaks, "inactive": False, "x": x.astype(int).tolist(), "t": T.round(2).tolist()}


# --------------------------------------------------------------------------------------
# Thermodynamics
# --------------------------------------------------------------------------------------
# average bond enthalpies, kJ/mol (typical textbook values)
BOND_ENERGY = {
    (("H", "H"), 1): 436,
    (("C", "H"), 1): 414, (("H", "N"), 1): 391, (("H", "O"), 1): 463,
    (("C", "C"), 1): 346, (("C", "C"), 2): 614, (("C", "C"), 3): 839,
    (("C", "N"), 1): 305, (("C", "N"), 2): 615, (("C", "N"), 3): 891,
    (("C", "O"), 1): 358, (("C", "O"), 2): 745, (("C", "O"), 3): 1072,
    (("N", "N"), 1): 163, (("N", "N"), 2): 418, (("N", "N"), 3): 945,
    (("N", "O"), 1): 201, (("N", "O"), 2): 607,
    (("O", "O"), 1): 146, (("O", "O"), 2): 498,
}
CO2_C_O = 799  # C=O in carbon dioxide
ATOMISATION = {"C": 717, "H": 218, "N": 473, "O": 249}  # kJ per mol of gaseous atoms
# standard molar entropies of the elements in their standard states, J/mol/K
S_C_GRAPHITE, S_H2, S_N2, S_O2 = 5.74, 130.68, 191.6, 205.2

KB = 1.380649e-23
H_PLANCK = 6.62607015e-34
R_GAS = 8.314462618


def bond_energy(e1, e2, order, co2_like=False):
    if order == 2 and co2_like and sorted((e1, e2)) == ["C", "O"]:
        return CO2_C_O, False
    key = (tuple(sorted((e1, e2))), order)
    if key in BOND_ENERGY:
        return BOND_ENERGY[key], False
    base = BOND_ENERGY.get((tuple(sorted((e1, e2))), 1), 300)
    return base * (1 + 0.8 * (order - 1)), True


def symmetry_number(elements, coords, tol=0.28) -> int:
    """Number of proper rotations mapping the molecule onto itself (rotational symmetry number)."""
    X = np.array(coords, float)
    n = len(X)
    if n == 1:
        return 1
    X = X - (X * np.array([MASS[e] for e in elements])[:, None]).sum(0) / sum(MASS[e] for e in elements)
    radii = np.linalg.norm(X, axis=1)
    # linear molecule?
    u, s, vt = np.linalg.svd(X - 0, full_matrices=False)
    if s[1] < 0.05:
        flipped = -X
        ok = all(
            any(elements[j] == elements[i] and np.linalg.norm(flipped[i] - X[j]) < tol for j in range(n))
            for i in range(n)
        )
        return 2 if ok else 1
    a0 = int(np.argmax(radii))
    # second reference atom: largest component perpendicular to a0
    e1 = X[a0] / np.linalg.norm(X[a0])
    perp = np.linalg.norm(X - np.outer(X @ e1, e1), axis=1)
    a1 = int(np.argmax(perp))

    def frame(p0, p1):
        f1 = p0 / np.linalg.norm(p0)
        f3 = np.cross(f1, p1)
        f3 = f3 / np.linalg.norm(f3)
        f2 = np.cross(f3, f1)
        return np.stack([f1, f2, f3], axis=1)

    Fu = frame(X[a0], X[a1])
    seen = set()
    for b0 in range(n):
        if elements[b0] != elements[a0] or abs(radii[b0] - radii[a0]) > tol:
            continue
        for b1 in range(n):
            if b1 == b0 or elements[b1] != elements[a1] or abs(radii[b1] - radii[a1]) > tol:
                continue
            if np.linalg.norm(np.cross(X[b0], X[b1])) < 1e-6:
                continue
            R = frame(X[b0], X[b1]) @ Fu.T
            Y = (R @ X.T).T
            ok = True
            for i in range(n):
                d = np.linalg.norm(X - Y[i], axis=1)
                d[[elements[j] != elements[i] for j in range(n)]] = 9
                if d.min() > tol:
                    ok = False
                    break
            if ok:
                seen.add(tuple(np.round(R, 1).flatten()))
    return max(1, len(seen))


def entropy_gas(elements, coords, stretch_nus, T=298.15, P=1e5):
    """Standard molar entropy (J/mol/K) of an ideal gas from mass, shape, symmetry, vibrations."""
    masses = np.array([MASS[e] for e in elements])
    M = masses.sum()
    n = len(elements)
    m = M * AMU
    s_trans = R_GAS * (math.log((2 * math.pi * m * KB * T / H_PLANCK**2) ** 1.5 * KB * T / P) + 2.5)
    s_rot = 0.0
    sigma = 1
    nvib = 0
    if n > 1:
        X = np.array(coords, float)
        com = (X * masses[:, None]).sum(0) / M
        Y = (X - com) * 1e-10
        I = np.zeros((3, 3))
        for mi, r in zip(masses * AMU, Y):
            I += mi * (np.dot(r, r) * np.eye(3) - np.outer(r, r))
        ev = np.sort(np.linalg.eigvalsh(I))
        sigma = symmetry_number(elements, coords)
        theta = [H_PLANCK**2 / (8 * math.pi**2 * x * KB) for x in ev if x > 1e-50]
        if ev[0] < 1e-3 * ev[-1] or len(theta) < 3:  # linear
            th = theta[-1] if theta else 1.0
            s_rot = R_GAS * (math.log(T / (sigma * th)) + 1)
            nvib = 3 * n - 5
        else:
            s_rot = R_GAS * (math.log(math.sqrt(math.pi) / sigma * T**1.5 / math.sqrt(theta[0] * theta[1] * theta[2])) + 1.5)
            nvib = 3 * n - 6
    nus = list(stretch_nus)[:nvib]
    nus += [800.0] * max(0, nvib - len(nus))  # bends / torsions: one soft generic frequency
    s_vib = 0.0
    for nu in nus:
        x = 1.438777 * nu / T
        s_vib += R_GAS * (x / math.expm1(x) - math.log(1 - math.exp(-x)))
    return s_trans + s_rot + s_vib, sigma, (s_trans, s_rot, s_vib)


def thermodynamics(molH, bond_info, stretch_nus, coords) -> dict:
    els = [a.GetSymbol() for a in molH.GetAtoms()]
    cnt = {e: els.count(e) for e in "CHNO"}
    x, y, z, w = cnt["C"], cnt["H"], cnt["N"], cnt["O"]

    # carbons with two C=O bonds are "CO2-like" -> use the stronger CO2 bond enthalpy
    co2_c = set()
    for atom in molH.GetAtoms():
        if atom.GetSymbol() == "C":
            n_co = sum(1 for b in bond_info if b["order"] == 2 and atom.GetIdx() in (b["a"], b["b"]) and set(b["els"]) == {"C", "O"})
            if n_co >= 2:
                co2_c.add(atom.GetIdx())
    total_bonds = 0.0
    estimated = False
    for b in bond_info:
        co2 = b["a"] in co2_c or b["b"] in co2_c
        e, est = bond_energy(b["els"][0], b["els"][1], b["order"], co2)
        b["energy"] = round(e)
        estimated |= est
        total_bonds += e
    atomisation = sum(ATOMISATION[e] * cnt[e] for e in "CHNO")
    hf = atomisation - total_bonds

    out = {
        "formula": hill_formula(cnt),
        "molar_mass": round(Descriptors.MolWt(molH), 2),
        "sum_bonds": round(total_bonds),
        "atomisation": round(atomisation),
        "hf": round(hf),
        "estimated_bonds": bool(estimated),
    }
    s_mol, sigma, parts = entropy_gas(els, coords, stretch_nus)
    out["entropy"] = round(s_mol, 1)
    out["sigma"] = sigma
    out["entropy_parts"] = [round(p, 1) for p in parts]

    # Formation reaction: elements in their standard states -> 1 mol of the molecule (gas)
    s_elements = x * S_C_GRAPHITE + (y / 2) * S_H2 + (z / 2) * S_N2 + (w / 2) * S_O2
    dsf = s_mol - s_elements

    def coeff(n2):  # n2 = number of half-units, e.g. 3 -> "1\u00bd"
        f = Fraction(n2, 2)
        if f == 0:
            return None
        if f.denominator == 1:
            return "" if f == 1 else str(f.numerator)
        whole, rem = divmod(f.numerator, f.denominator)
        return f"{whole if whole else ''}\u00bd"

    parts_eq = []
    for n2, formula in ((2 * x, "C(s)"), (y, "H2(g)"), (z, "N2(g)"), (w, "O2(g)")):
        c = coeff(n2)
        if c is not None:
            parts_eq.append(f"{c} {formula}".strip())
    out["formation"] = {
        "equation": f"{' + '.join(parts_eq)} \u2192 {out['formula']}(g)",
        "ds": round(dsf, 1),
        "delta_n_gas": float(1 - (y / 2 + z / 2 + w / 2)),
        "ds_sign": "increases" if dsf > 5 else "decreases" if dsf < -5 else "barely changes",
    }
    return out


# --------------------------------------------------------------------------------------
# Everything for the data viewer
# --------------------------------------------------------------------------------------
def compute_all(graph: dict, traj: dict) -> dict:
    molH, order_map = build_mol(graph)
    coords = np.array(traj["frames"][-1])
    n = molH.GetNumAtoms()
    cid = _embed(molH, 1)  # placeholder conformer, then overwrite with the optimised coords
    conf = molH.GetConformer(cid)
    for i in range(n):
        conf.SetAtomPosition(i, coords[i].tolist())

    q = partial_charges(molH)
    els = [a.GetSymbol() for a in molH.GetAtoms()]

    bond_info = []
    for b in molH.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        order = order_map.get((min(i, j), max(i, j)), 1)
        length = float(np.linalg.norm(coords[i] - coords[j]))
        bond_info.append({
            "a": i, "b": j, "order": order, "aromatic": b.GetIsAromatic(),
            "els": (els[i], els[j]), "length": round(length, 3),
            "dq": round(q[i] - q[j], 3),
        })

    # dipole: mu = sum q r
    mu_vec = (np.array(q)[:, None] * coords).sum(0) * DEBYE_PER_EANG
    mu = float(np.linalg.norm(mu_vec))
    masses = np.array([MASS[e] for e in els])
    com = (coords * masses[:, None]).sum(0) / masses.sum()
    # points to where negative charge is concentrated; draw from the centre of mass
    dipole = {"vec": mu_vec.round(3).tolist(), "mag": round(mu, 2), "origin": com.round(3).tolist()}

    ir = ir_spectrum(molH, bond_info, q)
    nus = [_wavenumber(_force_k(b["els"][0], b["els"][1], b["order"], b["aromatic"]), MASS[b["els"][0]], MASS[b["els"][1]]) for b in bond_info]
    thermo = thermodynamics(molH, bond_info, nus, coords)

    return {
        "formula": thermo["formula"],
        "elements": els,
        "coords": coords.round(4).tolist(),
        "charges": [round(c, 3) for c in q],
        "bonds": [{k: (list(v) if isinstance(v, tuple) else v) for k, v in b.items()} for b in bond_info],
        "dipole": dipole,
        "ir": ir,
        "thermo": thermo,
    }
