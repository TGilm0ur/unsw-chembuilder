"""
UNSW ChemBuilder - a simple molecule builder + quantum-chemistry-style data viewer.

Run with:   streamlit run app.py
"""
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

import chem_engine as ce
from syllabus import SYLLABUS

HERE = Path(__file__).parent
UNSW_YELLOW = "#FFD100"  # change here (and in the two component CSS files) to retune the brand colour

st.set_page_config(page_title="UNSW ChemBuilder", page_icon="⚛", layout="wide", initial_sidebar_state="collapsed")

_methods = components.declare_component("unsw_methods", path=str(HERE / "components" / "methods"))
_viewer = components.declare_component("unsw_viewer", path=str(HERE / "components" / "viewer"))

CSS = f"""
<style>
  #MainMenu, footer, header[data-testid="stHeader"], [data-testid="stToolbar"], [data-testid="stDecoration"] {{ display:none !important; }}
  .stApp {{ background:#a3a199 repeating-linear-gradient(45deg, #a3a199 0 6px, #aaa8a0 6px 12px); }}
  .block-container {{ padding:14px 22px 10px 22px !important; max-width:1500px; }}
  html, body, .stApp, .stApp *:not(code):not(pre) {{ font-family: Tahoma, Verdana, "MS Sans Serif", sans-serif; }}
  .appbar {{ background:linear-gradient(90deg,{UNSW_YELLOW},#ffe566); border:2px solid; border-color:#fff #404040 #404040 #fff; box-shadow:1px 1px 0 #000;
            padding:6px 12px; display:flex; justify-content:space-between; align-items:baseline; margin-bottom:18px; color:#000; }}
  .appbar b {{ font-size:21.33px; letter-spacing:.5px; }}
  .appbar span {{ font-size:13.33px; }}
  .tbar {{ background:#d4d0c8; border:2px solid; border-color:#fff #404040 #404040 #fff; padding:6px 8px; margin-bottom:8px; }}
  div.stButton > button {{ background:#d4d0c8; color:#000; border:2px solid; border-color:#fff #404040 #404040 #fff; border-radius:0;
            font-family:inherit; font-weight:700; padding:6px 14px; font-size:calc(1rem + 1pt); box-shadow:none; }}
  div.stButton > button:hover {{ background:#e6e2d9; color:#000; border-color:#fff #404040 #404040 #fff; }}
  div.stButton > button:active {{ border-color:#404040 #fff #fff #404040; }}
  div.stButton > button[kind="primary"] {{ background:{UNSW_YELLOW}; color:#000; }}
  div.stButton > button[kind="primary"]:hover {{ background:#ffe14d; }}
  .note {{ background:#fff3b0; border:2px solid #a58b00; padding:5px 10px; margin:6px 0; font-size:14.33px; color:#000; }}
  .err  {{ background:#ffd9d9; border:2px solid #a00000; padding:5px 10px; margin:6px 0; font-size:14.33px; color:#000; }}
  iframe {{ border:none; }}
  .intro {{ background:#d4d0c8; border:2px solid; border-color:#fff #404040 #404040 #fff; box-shadow:2px 2px 0 #000; color:#000; max-width:1150px; margin:0 auto 10px auto; }}
  .intro .it {{ background:linear-gradient(90deg,{UNSW_YELLOW},#ffe566); border-bottom:2px solid #000; padding:4px 10px; font-weight:bold; font-size:16.33px; }}
  .intro .ib {{ padding:12px 20px 16px 20px; font-size:15.33px; line-height:1.55; }}
  .intro h3 {{ margin:14px 0 6px 0; font-size:18.33px; border-bottom:2px solid {UNSW_YELLOW}; padding-bottom:2px; }}
  .intro table {{ border-collapse:collapse; width:100%; margin:6px 0; font-size:14.33px; }}
  .intro th {{ background:#bdb9b0; border:1px solid #808080; padding:4px 8px; text-align:left; }}
  .intro td {{ border:1px solid #a8a49c; padding:4px 8px; vertical-align:top; background:#e6e2d9; }}
  .intro .disc {{ background:#fff3b0; border:2px solid #a58b00; padding:8px 12px; margin:8px 0; }}
  .intro ol li, .intro ul li {{ margin-bottom:4px; }}
  .intro .cols {{ display:grid; grid-template-columns:1fr 1fr; gap:26px; padding-top:4px; padding-bottom:12px; }}
  .intro .cols p {{ margin:4px 0; }}
  .intro .cols h3 {{ margin-top:6px; }}
  .intro .cols .disc {{ margin:4px 0; padding:6px 10px; }}
  .intro .cols ol {{ margin:4px 0; padding-left:22px; }}
  [data-testid="stVerticalBlock"] {{ gap:0.5rem; }}
</style>
"""

ss = st.session_state
ss.setdefault("page", "intro")
ss.setdefault("return_to", None)   # where the About pages should send the student back to
ss.setdefault("graph", {"atoms": [], "bonds": []})
ss.setdefault("opt", None)       # {"key", "play"} once the student has pressed Geometry Optimise
ss.setdefault("play", 0)
ss.setdefault("results", None)
ss.setdefault("auto_opt", False)


@st.cache_data(show_spinner=False, max_entries=64)
def _trajectory(graph_json: str):
    import json
    return ce.trajectory(json.loads(graph_json))


@st.cache_data(show_spinner=False, max_entries=64)
def _results(graph_json: str):
    import json
    g = json.loads(graph_json)
    return ce.compute_all(g, ce.trajectory(g))


def _topology_json(g: dict) -> str:
    import json
    # positions in the 2D sketch do not affect the chemistry, so leave them out of the cache key
    return json.dumps({"atoms": [{"id": a["id"], "el": a["el"]} for a in g["atoms"]],
                       "bonds": [{"a": b["a"], "b": b["b"], "o": b["o"]} for b in g["bonds"]]}, sort_keys=True)


def appbar(title: str, right: str):
    st.markdown(f'<div class="appbar"><b>⚛ UNSW ChemBuilder</b><span>{title} &nbsp;|&nbsp; {right}</span></div>', unsafe_allow_html=True)


# --------------------------------------------------------------------------------------


def _continue_button(col):
    """Primary button at the bottom of the About pages: returns to wherever the student came from."""
    back = ss.return_to if (ss.return_to == "builder" or (ss.return_to == "viewer" and ss.results)) else None
    label = {"viewer": "\u2190  Back to data viewer", "builder": "\u2190  Back to builder"}.get(back, "\u25b6  Start building")
    if col.button(label, type="primary", use_container_width=True):
        ss.page = back or "builder"
        ss.return_to = None
        st.rerun()


def intro_page():
    appbar("Welcome", "A simple test run of quantum chemistry workflow")
    st.markdown("""
<div class="intro">
<div class="it">About UNSW ChemBuilder</div>
<div class="ib cols">
<div>
<h3>What is this?</h3>
<p>Programs like <b>IQMol</b> are used by chemists to predict molecular properties with quantum chemistry.
<b>ChemBuilder is a simple, accessible test run of that workflow:</b></p>
<p style="text-align:center;font-weight:bold;margin:4px 0">Build &rarr; Optimise geometry &rarr; Explore the data</p>
<h3>Please note</h3>
<div class="disc"><b>No quantum chemistry is used here.</b> Results come from fast, classical approximations, so they show the right
trends but should be treated as estimates. For example, <b>Geometry Optimise</b> actually runs a molecular mechanics <i>energy minimisation</i>.
&lsquo;Geometry optimisation&rsquo; usually refers to a quantum calculation done on a supercomputer.</div>
</div>
<div>
<h3>How to use it</h3>
<ol>
<li><b>Build</b> a small molecule (C, H, N, O; up to 4 non-hydrogen atoms).</li>
<li>Press <b>Geometry Optimise</b> to watch it move into a lower-energy shape.</li>
<li>Press <b>Calculate</b> and explore the data, then use <b>Connect to the Syllabus</b> to link each output to Year 11 Chemistry.</li>
<li><b>Check the accuracy:</b> search for experimental or higher-accuracy data for your molecule (e.g. the
<a href="https://webbook.nist.gov/chemistry/" target="_blank">NIST Chemistry WebBook</a>) and compare. Which outputs are closest, and why?</li>
</ol>
</div>
</div></div>
""", unsafe_allow_html=True)
    _, b1, b2, _ = st.columns([1.2, 1.3, 1.3, 1.2])
    if b1.button("\u2139  Approximations used", use_container_width=True):
        ss.page = "methods"
        st.rerun()
    _continue_button(b2)


def methods_page():
    appbar("Approximations used", "How each value is estimated")
    _methods(height=1500, default=None, key="methods_doc")

    _, b1, b2, _ = st.columns([1.2, 1.3, 1.3, 1.2])
    if b1.button("\u2190  Back to about", use_container_width=True):
        ss.page = "intro"
        st.rerun()
    _continue_button(b2)


def builder_page():
    appbar("Molecule Builder", "Build &rarr; Geometry Optimise &rarr; Calculate")

    # The sketcher and the two buttons live inside the viewer component. Its latest value is already
    # in session_state when the script starts, so everything is computed in ONE run (no st.rerun,
    # no page flash).
    val = ss.get("builder")
    action = None
    if isinstance(val, dict) and "atoms" in val:
        ss.graph = {"atoms": val["atoms"], "bonds": val["bonds"]}
        if val.get("action") and val.get("seq") != ss.get("last_seq"):
            action = val["action"]
        ss.last_seq = val.get("seq")

    errors, warnings = ce.validate(ss.graph)
    key = ce.topology_key(ss.graph) if not errors else None
    if ss.opt and ss.opt["key"] != key:
        ss.opt = None

    if action == "about":
        ss.return_to = "builder"
        ss.page = "intro"
        st.rerun()
    if action == "optimise" and not errors:
        ss.play += 1
        ss.opt = {"key": key, "play": ss.play}
    if action == "calculate" and not errors:
        ss.auto_opt = ss.opt is None
        ss.results = _results(_topology_json(ss.graph))
        ss.page = "viewer"
        st.rerun()

    shown_errors = errors if ss.graph["atoms"] else []
    common = dict(initial=ss.graph, errors=shown_errors, warnings=warnings, disabled=bool(errors),
                  height=800, default=None, key="builder")
    if errors:
        _viewer(mode="builder", empty=True, **common)
    else:
        traj = _trajectory(_topology_json(ss.graph))
        optimised = ss.opt is not None
        status = ("Geometry optimised (molecular mechanics)." if optimised
                  else "Press 'Geometry Optimise' to minimise the energy of the structure.")
        _viewer(mode="builder", empty=False, elements=traj["elements"], bonds=traj["bonds"], frames=traj["frames"],
                mol_key=traj["key"], optimised=optimised, play_id=ss.opt["play"] if optimised else 0,
                status=status, **common)


def viewer_page():
    res = ss.results
    appbar("Data Viewer", res["formula"])
    c1, c2, c3, c4 = st.columns([1.3, 1.3, 1, 4])
    if c1.button("←  Back to Builder", use_container_width=True):
        ss.page = "builder"
        st.rerun()
    if c2.button("✖  New Molecule", use_container_width=True):
        ss.graph = {"atoms": [], "bonds": []}
        ss.opt = None
        ss.results = None
        ss.return_to = None
        ss.page = "builder"
        st.rerun()
    if c3.button("\u2139  About", use_container_width=True):
        ss.return_to = "viewer"
        ss.page = "intro"
        st.rerun()
    c4.empty()
    _, warnings = ce.validate(ss.graph)
    for w in warnings:
        st.markdown(f'<div class="note">⚠ {w}</div>', unsafe_allow_html=True)
    res = dict(res, auto_optimised=ss.auto_opt)
    _viewer(mode="data", data=res, syllabus=SYLLABUS, height=800, mol_key=ce.topology_key(ss.graph), default=None, key="dataview")


st.markdown(CSS, unsafe_allow_html=True)
if ss.page == "intro":
    intro_page()
elif ss.page == "methods":
    methods_page()
elif ss.page == "viewer" and ss.results:
    viewer_page()
else:
    ss.page = "builder"
    builder_page()
