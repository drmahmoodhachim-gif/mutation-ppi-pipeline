# LOCKED VERSION — Mutation PPI Pipeline

## 🔒 PERFECT WORKING VERSION — NEVER EVER CHANGE

**Commit:** `aa8f251`  
**Tag:** `v1.0-perfect` (if created)  
**Date locked:** 2026-02-27

**This version is perfect. Do not modify app.py, config.py, predictors.py, or visualization.py.**

---

## Critical rules (DO NOT CHANGE)

1. **app.py**
   - Import: `from visualization import fetch_pdb, render_py3dmol_html`
   - Import: `get_ppi_ddg_predictions` from predictors (REQUIRED for Section 4)
   - **NEVER** use `create_mol_viewer` in app.py (causes NameError)
   - Welcome banner: **Hi Inna Aleksandrova** · مرحباً إينا ألكساندروفا · Привет, Инна Александрова (MUST remain after sub-header)
   - Welcome banner: "Hi Inna Aleksandrova" (Arabic, English, Russian) must remain
   - Section 6: PPI ΔΔG (wild vs mutant) for known interactors

2. **predictors.py**
   - AlphaMissense: parse `benign` / `pathogenic` / `ambiguous` (API format changed)
   - Import: NO `MMCSM_PPI_API` (not in config)
   - Must export: `get_ppi_ddg_predictions`

3. **config.py**
   - Must have: `PPI_PDB_COMPLEXES`, `VCL` in GENE_UNIPROT and PROTEIN_PDB

4. **visualization.py**
   - Must have: `render_py3dmol_html`, `_extract_chain` for chain A only
   - Must have: `fetch_pdb` with timeout=30

5. **requirements.txt**
   - NO `stmol` (causes Streamlit Cloud build issues)
   - Include: streamlit, py3Dmol, requests, pandas, biopython

## Included features

- Welcome: Hi Inna Aleksandrova | مرحباً إينا ألكساندروفا | Привет, Инна Александрова
- AlphaMissense pathogenicity (benign/pathogenic/ambiguous API)
- Structural impact heuristic
- Tissue-specific PPIs
- Side-by-side 3D: wild-type (blue) vs mutant (red), chain A only
- PPI ΔΔG: wild vs mutant for each known interactor
- VCL (Vinculin) in config
- mCSM-PPI2 link for manual predictions

## DO NOT CHANGE — This script works perfectly

- Example variant dropdown (SCN5A p.R526H, p.Ser1054Ala, VCL p.Gln836Arg)
- Conclusion dropdown (Decreased / Increased / No change)
- PPI ΔΔG table with Interacting protein, Role, Wild-type, Mutant ΔΔG, Conclusion, Pathway effect
- 3-letter mutation parsing (p.Ser1054Ala, p.Gln836Arg)

## Recovery

If broken: `git checkout aa8f251` or `git checkout v1.0-perfect` to restore this perfect version.
