# UNSW ChemBuilder

A simple IQMol-inspired molecule builder and data viewer for Year 11 Chemistry.

```bash
pip install -r requirements.txt
streamlit run app.py
```

1. **Builder** – sketch a molecule (C, H, N, O; max 4 non-H atoms; hydrogens fill in automatically).
2. **Geometry Optimise** – watch the structure relax (molecular mechanics energy minimisation, visual only).
3. **Calculate** – opens the **Data Viewer**: bond lengths, partial charges, dipole moment, IR spectrum, thermodynamics.
4. **Connect to the Syllabus** – explainer panel for each output (edit `syllabus.py`).

## Files
| File | Purpose |
|---|---|
| `app.py` | Streamlit app, page flow, retro styling |
| `chem_engine.py` | All chemistry (RDKit + numpy) - swap in literature values here later |
| `syllabus.py` | Year 11 syllabus text for each output |
| `components/viewer/` | One component: 2D sketcher, 3D builder view and data viewer (3Dmol.js; needs internet for the CDN script) |
| `components/methods/` | "Approximations used" page (KaTeX equations; needs internet for the CDN script) |

Brand colour: `#FFD100` (search for it to change; it is in `app.py` and both component files).

## How the calculations work

No quantum chemistry is used. Everything is in `chem_engine.py` (RDKit + numpy), at 298.15 K.
Results are computed on the optimised geometry unless stated otherwise, and cached by molecule topology in `app.py`.

### Validation and structure
- Only C, H, N, O; at most 4 heavy atoms (`MAX_HEAVY`); `validate` returns errors (block the buttons) or warnings (note only).
- 3D structure: RDKit ETKDGv3 embedding.

### Geometry Optimise (force fields)
- Both force fields are **built into RDKit**; this code only chooses between them and calls the minimiser.
- Default: **MMFF94** (`MMFFGetMoleculeForceField`), used when RDKit has parameters for every atom.
- Fallback: **UFF** (`UFFGetMoleculeForceField`) with a distance restraint on every bond toward `IDEAL_LENGTH` (force constant 20000). Used when MMFF lacks parameters **or** `_mmff_unreliable` is true.
- `_mmff_unreliable` is a hand-written rule on the bond graph (it never inspects MMFF output): true if any atom has two or more double bonds (CO2-like), or any bond is between two N/O atoms. Chosen because MMFF gave CO2 a C=O of 1.405 A and N2 a bond of 1.46 A. It is crude: some caught molecules might be fine with MMFF, and some bad cases may slip through.
- `IDEAL_LENGTH` (typical bond lengths) was entered from memory; a special case sets 1.16 A for C=O on a CO2-like carbon.
- Animation: the start geometry is jittered (sigma 0.32 A) and stretched, then minimised 3 iterations at a time (up to 1500), saving each frame, then polished. At most 40 frames are kept, Kabsch-aligned to the final geometry.

### Bond lengths
Measured from the optimised geometry. For molecules that use the UFF fallback, lengths are largely set by the `IDEAL_LENGTH` restraints.

### Partial charges
**Gasteiger-Marsili** (`ComputeGasteigerCharges`, built into RDKit): electronegativity equalisation along bonds, iterated with damping (RDKit defaults). Electronegativity parameters depend on atom type/hybridisation. Uses the bond graph only, so charges do not change on optimisation. MMFF charges are not used because MMFF sets every C-H bond to exactly zero.

### Dipole moment
mu = sum(q_i r_i) x 4.8032 D per e.A, from the Gasteiger charges and optimised positions. The arrow is drawn from the centre of mass. mu points toward the positive end; the viewer flips it so the arrow points toward the negative end (chemistry/syllabus convention). Below 0.05 D it is shown as "no net dipole"; homonuclear diatomics (H2, N2, O2) have none. Individual bond dipoles are shown for bonds with a charge difference of at least 0.12 e. Point charges only (no lone pairs), so water and ammonia dipoles are underestimated.

### IR spectrum
Harmonic oscillator per bond: wavenumber = (1/2 pi c) sqrt(k / mu), mu = reduced mass (masses from the `MASS` table in the code, not RDKit).
- Force constants `FORCE_K` (N/m, by element pair and bond order; aromatic 800): **written from memory as typical textbook-range values, not from a cited source and not fitted to spectra.** Unlisted bond orders scale the single-bond value; fallback base 450.
- Intensity = 0.08 + |q1 - q2| (Gasteiger charge difference as a proxy for bond polarity). **A heuristic, not derived**: real intensities depend on the dipole derivative. The 0.08 floor keeps non-polar bonds visible.
- Two or more X-H bonds on one heavy atom become a symmetric/asymmetric pair: peaks at the mean minus 45 and plus 55 cm^-1, sharing 45% / 55% of the intensity (fixed, hand-picked numbers). The animation moves the bonds in phase / out of phase.
- Peaks within about 4 cm^-1 with the same label are merged. Line shape is Lorentzian (width 70 cm^-1 for O-H, 35 for N-H, 14 otherwise); transmittance = 100 exp(-1.6 A) from 4000 to 400 cm^-1 (appearance constants).
- Stretches only: no bends, twists or coupling. Known errors: CO2 asymmetric stretch ~1723 cm^-1 (real ~2349), ethanol C-C stretch too high.
- Homonuclear diatomics return no peaks and `inactive: True` ("IR inactive").

### Thermodynamics
- **Delta H_f** (Hess's Law) = sum of atomisation enthalpies - sum of bond enthalpies. Atomisation (kJ/mol): C 717, H 218, N 473, O 249. Bond enthalpies from `BOND_ENERGY` (e.g. C-H 414, O-H 463, C=O 745); a carbon with two C=O bonds uses 799 (CO2). Unlisted bonds: single-bond value scaled by +80% per extra order (300 if no single-bond entry) and flagged `estimated_bonds`. **All of these tables were entered from memory as typical textbook values, with no cited source; not yet cross-checked.**
- **S** = S_trans + S_rot + S_vib, ideal gas at 1 bar. S_trans is the translational part of the Sackur-Tetrode equation using the total molecular mass (the molecule's centre of mass as one particle); S_rot is a rigid rotor from the moments of inertia (linear or non-linear formula, divided by the rotational symmetry number); S_vib is a harmonic oscillator sum over 3N-5 (linear) or 3N-6 modes, using the IR stretch frequencies and **one generic 800 cm^-1 for all bends/torsions** (a chosen number). Missing electronic term: O2 (triplet) is about 9 J/mol/K too low (195.9 vs 205.2).
- **Delta S_f** = S(molecule) - standard entropies of the elements: C(graphite) 5.74, H2 130.68, N2 191.6, O2 205.2 J/mol/K (from memory).
- **Delta G_f** = Delta H_f - T Delta S_f, assuming Delta H and Delta S do not change with temperature; the viewer slider recomputes it in the browser. Combustion enthalpy is deliberately not shown.

### To do before relying on the numbers
Cross-check `FORCE_K`, `BOND_ENERGY`, `ATOMISATION`, element entropies and `IDEAL_LENGTH` against a cited source (e.g. NIST CCCBDB / CRC Handbook), and decide whether to fix the O2 spin entropy and the CO2 asymmetric stretch.

All numbers are teaching-level estimates and are labelled as such in the UI.
