# UNSW ChemBuilder

A simple WebMO-inspired molecule builder and data viewer for Year 11 Chemistry.

```bash
pip install -r requirements.txt
streamlit run app.py
```

1. **Builder** – sketch a molecule (C, H, N, O; max 4 non-H atoms; hydrogens fill in automatically).
2. **Geometry Optimise** – watch the structure relax (MMFF94 molecular mechanics, visual only).
3. **Calculate** – opens the **Data Viewer**: bond lengths, partial charges, dipole moment, IR spectrum, thermodynamics.
4. **Connect to the Syllabus** – explainer panel for each output (edit `syllabus.py`).

## Files
| File | Purpose |
|---|---|
| `app.py` | Streamlit app, page flow, retro styling |
| `chem_engine.py` | All chemistry (RDKit + numpy) - swap in literature values here later |
| `syllabus.py` | Year 11 syllabus text for each output |
| `components/sketcher/` | Custom 2D click-and-drag molecule sketcher |
| `components/viewer/` | 3Dmol.js viewer + data viewer window (needs internet for the 3Dmol CDN script) |

Brand colour: `#FFD100` (search for it to change; it is in `app.py` and both component files).

All numbers are teaching-level estimates (Gasteiger charges, spring-model IR, average bond
enthalpies, statistical-mechanics entropy) and are labelled as such in the UI.
