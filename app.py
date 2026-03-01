"""
Mutation PPI Prediction Pipeline — Interactive Streamlit App
Input: gene, mutation, tissue → Variant QC, pathogenicity, structural impact, PPI (interface-aware)
"""

import streamlit as st
import pandas as pd
from config import GENE_UNIPROT
try:
    from config import APP_AUTHOR
except ImportError:
    APP_AUTHOR = {}
try:
    from config import APP_PUBLIC_URL
except ImportError:
    APP_PUBLIC_URL = ""
try:
    from config import EXAMPLE_VARIANTS
except ImportError:
    EXAMPLE_VARIANTS = []
from predictors import (
    parse_mutation,
    get_alphamissense_prediction,
    get_tissue_interactors,
    get_recommended_pdb,
    estimate_structural_impact,
    resolve_uniprot,
    validate_variant_on_sequence,
)
try:
    from predictors import get_ppi_ddg_predictions
except ImportError:
    get_ppi_ddg_predictions = None
from visualization import fetch_pdb, render_py3dmol_html

# Page config
st.set_page_config(
    page_title="Mutation PPI Pipeline",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
    .banner {
        background: linear-gradient(135deg, #1E88E5 0%, #1565C0 50%, #0D47A1 100%);
        padding: 1.5rem 1.5rem 1.8rem;
        border-radius: 12px;
        margin-bottom: 1rem;
        box-shadow: 0 4px 12px rgba(30, 136, 229, 0.3);
    }
    .banner-title { font-size: 2.2rem; font-weight: 700; color: white; margin: 0; letter-spacing: -0.5px; }
    .banner-subtitle { font-size: 1rem; color: rgba(255,255,255,0.9); margin: 0.4rem 0 0 0; }
    .author-box {
        background: #f8f9fa; border-left: 4px solid #1E88E5;
        padding: 0.8rem 1rem; border-radius: 0 8px 8px 0;
        margin-bottom: 1rem; font-size: 0.95rem;
    }
    .author-name { font-weight: 600; color: #1565C0; }
</style>
""", unsafe_allow_html=True)

# Banner
st.markdown("""
<div class="banner">
    <p class="banner-title">🧬 Mutation PPI Prediction Pipeline</p>
    <p class="banner-subtitle">Variant QC → pathogenicity, structural impact, PPI (interface-aware)</p>
</div>
""", unsafe_allow_html=True)

# Welcome & author details
auth = APP_AUTHOR or {}
if auth.get("name"):
    title_part = f" — {auth['title']}" if auth.get("title") else ""
    aff_part = auth.get("affiliation", "")
    email_part = f' • <a href="mailto:{auth["email"]}">{auth["email"]}</a>' if auth.get("email") else ""
    st.markdown(
        f'<div class="author-box">'
        f'<strong>Welcome to the tool.</strong><br>'
        f'Developed by <span class="author-name">{auth["name"]}</span>{title_part}<br>'
        f'{aff_part}{email_part}'
        f'</div>',
        unsafe_allow_html=True
    )

# Disclaimer — always visible
st.info(
    "**Disclaimer:** For research and educational use only. Not medical advice. Predictions are theoretical models. "
    "Do not use as a substitute for professional medical advice, diagnosis, or treatment."
)

# Data sources, attributions & legal
with st.expander("📋 Data Sources, Attributions & Legal", expanded=False):
    st.markdown("""
    **Disclaimer:** This tool is for research and educational use only. Predictions are theoretical models and do not constitute medical or clinical advice. Do not use as a substitute for professional medical advice, diagnosis, or treatment.

    **Educational use:** This pipeline uses publicly available data from the following sources. Attribution and good scientific practice apply.
    """)
    st.markdown("""
    | Source | Data / Purpose | License / Terms |
    |--------|----------------|-----------------|
    | [UniProt](https://www.uniprot.org/) | Canonical sequences, protein metadata | CC BY-ND; [attribution expected](https://www.uniprot.org/help/license) |
    | [AlphaFold DB](https://alphafold.ebi.ac.uk/) (EMBL-EBI) | Protein structure predictions, pLDDT | [EMBL-EBI terms](https://www.ebi.ac.uk/about/terms-of-use); cite [AlphaFold DB](https://www.ebi.ac.uk/training/online/courses/navigating-alphafold-database/citing-the-database/) |
    | [AlphaMissense](https://alphamissense.hegelab.org/) | Pathogenicity scores | CC BY-NC-SA (non-commercial); [disclaimers apply](https://alphamissense.hegelab.org/help) |
    | [RCSB PDB](https://www.rcsb.org/) | Experimental structures | CC0; [attribution encouraged](https://www.rcsb.org/pages/usage-policy) |
    | [mCSM-PPI2](https://biosig.lab.uq.edu.au/mcsm_ppi2/) | PPI ΔΔG predictions (external link) | [Rodrigues et al., NAR 2019](https://doi.org/10.1093/nar/gkz383) |
    """)
    st.caption("Property and intellectual rights remain with the respective data providers. Users must comply with each source's terms when publishing or redistributing results.")

# Sidebar — Input
with st.sidebar:
    st.header("📥 Input")
    with st.expander("👤 For server admin (you)", expanded=False):
        if APP_PUBLIC_URL:
            st.markdown(f"**Share link:** [{APP_PUBLIC_URL}]({APP_PUBLIC_URL})")
        st.markdown("""
        **Run on your laptop:**
        1. `docker pull ghcr.io/sokrypton/colabfold:1.5.3-cuda12.2.2`
        2. `$env:COLABFOLD_DOCKER="1"; python -m streamlit run app.py --server.port 8501`
        3. `ngrok http 8501`

        See RUN_LOCAL.md for full steps.
        """)
    with st.expander("How to use", expanded=False):
        st.markdown("""
        1. Enter **any gene** (we look up UniProt) and **mutation** (text: `p.R526H` or simple: position + WT + Mut)
        2. Choose **tissue** and **analysis mode**
        3. Click **Run Pipeline**
        4. Interacting proteins come from STRING for any gene (curated for cardiac when available)
        5. If mutant is missing: use **Predict here** (≤400 aa) or **Option B**
        """)
    gene = st.text_input("Gene symbol", value="SCN5A", help="Any gene; we look up UniProt automatically")
    input_mode = st.radio("Mutation input", ["Text format", "Position + AAs"], horizontal=True, help="Text: p.R526H, R526H, Ser1054Ala. Simple: enter position, WT, Mut separately.")
    pos_override, wt_override, mut_override = None, None, None
    mutation_input = ""
    if input_mode == "Text format":
        mutation_input = st.text_input("Mutation", value="c.1577G>A, p.R526H", help="p.R526H, R526H, Ser1054Ala, c.1577G>A, or 526 R→H")
    else:
        c1, c2, c3 = st.columns(3)
        with c1:
            pos_override = st.number_input("Position", min_value=1, value=526, help="Amino acid position in protein")
        aa_list = list("ARNDCEQGHILKMFPSTWYV")
        with c2:
            wt_override = st.selectbox("WT AA", aa_list, index=aa_list.index("R"), help="Wild-type amino acid")
        with c3:
            mut_override = st.selectbox("Mutant AA", aa_list, index=aa_list.index("H"), help="Mutant amino acid")
        mutation_input = f"{pos_override} {wt_override} {mut_override}"
    tissue = st.selectbox(
        "Tissue of interest",
        ["Cardiac myocyte", "Heart", "Skeletal muscle", "Neuron", "Other"],
        index=0,
    )
    analysis_mode = st.radio(
        "Analysis mode",
        ["Fast (no AlphaFold)", "AlphaFold (WT vs MUT)", "AlphaFold-Multimer for PPI"],
        index=0,
    )
    run_button = st.button("🚀 Run Pipeline", type="primary", use_container_width=True)

# Parse mutation (with gene for variant_id)
parsed = parse_mutation(
    mutation_input, gene,
    position=pos_override if input_mode == "Position + AAs" else None,
    wt_aa=wt_override if input_mode == "Position + AAs" else None,
    mut_aa=mut_override if input_mode == "Position + AAs" else None,
)
pos = parsed.get("position") or parsed.get("cds_pos")
wt, mut = parsed.get("wt_aa"), parsed.get("mut_aa")

if not pos and not parsed.get("cds_pos"):
    st.warning("Could not parse mutation. Use p.R526H, R526H, 526 R→H, or Position + AAs.")
    st.stop()

# Resolve UniProt → canonical sequence
uniprot_data = resolve_uniprot(gene)
if isinstance(uniprot_data, dict):
    uniprot_id = uniprot_data.get("uniprot_id")
    canonical_fasta = uniprot_data.get("canonical_fasta", "")
else:
    uniprot_id = uniprot_data or GENE_UNIPROT.get(gene.upper())
    canonical_fasta = ""

if not uniprot_id:
    st.error(f"Could not resolve UniProt ID for gene {gene}")
    st.stop()

# Mandatory Variant QC — hard-stop on WT mismatch
st.header("0️⃣ Variant QC & Mapping")
qc_ok = True
if wt and pos and canonical_fasta:
    qc = validate_variant_on_sequence(canonical_fasta, wt, pos)
    if not qc.get("match"):
        st.error(f"**WT mismatch:** {qc.get('message', '')}")
        st.stop()
else:
    if not wt or not mut:
        st.warning("WT/mutant not fully parsed. Use p.R526H format.")
        qc_ok = False

# PDB coverage (with UniProt->PDB mapping when canonical available; fetches from RCSB for any gene)
pdb_rec = get_recommended_pdb(gene, pos, uniprot_seq=canonical_fasta, uniprot_id=uniprot_id)
has_res = pdb_rec.get("has_residue_coordinates", True)
pdb_ids = pdb_rec.get("pdb_ids", [])
pdb_resseq = pdb_rec.get("pdb_resseq")
mapping_method = pdb_rec.get("mapping_method", "—")

col1, col2, col3, col4, col5 = st.columns(5)
with col1:
    st.metric("UniProt ID", uniprot_id)
with col2:
    st.metric("Canonical length", len(canonical_fasta) or "—")
with col3:
    st.metric("Verified WT at pos", f"{wt or ''}{pos}" if wt and pos else "—")
with col4:
    st.metric("Residue in PDB?", "Yes" if has_res else "No")
with col5:
    st.metric("WT check", "✅ OK" if qc_ok else "⚠️ Warn")
if pdb_resseq is not None or mapping_method not in ("—", "config"):
    st.caption(f"**Mapping:** {mapping_method} · **Mapped PDB resseq:** {pdb_resseq if pdb_resseq is not None else '—'}")

if not (run_button or st.session_state.get("results_ready")):
    st.info("👈 Click **Run Pipeline** to continue.")
    st.stop()

st.session_state["results_ready"] = True

with st.spinner("Running predictions..."):
    # Section 1: Parsed mutation
    st.header("1️⃣ Parsed Mutation")
    st.metric("Variant", parsed.get("variant_id", f"{gene}:{wt}{pos}{mut}"))

    # Section 2: AlphaMissense
    st.header("2️⃣ AlphaMissense Pathogenicity")
    if pos and wt and mut:
        am_result = get_alphamissense_prediction(gene, pos, wt, mut)
        if "error" in am_result:
            st.warning(f"AlphaMissense API: {am_result['error']}")
        else:
            score = am_result.get("pathogenicity")
            if score is not None:
                st.error(f"Pathogenic: {score:.1%}") if score > 0.5 else st.success(f"Pathogenic: {score:.1%}")

    # AlphaFold (if selected)
    af_deltas = None
    if analysis_mode.startswith("AlphaFold") and canonical_fasta and pos and mut:
        try:
            from alphafold_runner import run_wt_mut_alphafold
            from structure_metrics import compute_local_structure_deltas_from_af
            import os
            variant_id = parsed.get("variant_id", f"{gene}_{wt}{pos}{mut}")
            mut_seq = canonical_fasta[: pos - 1] + mut + canonical_fasta[pos:]
            out_dir = os.path.join(os.path.dirname(__file__), "alphafold_out")
            af_result = run_wt_mut_alphafold(canonical_fasta, mut_seq, variant_id, out_dir, uniprot_id=uniprot_id)
            wt_ok = af_result["wt"]["status"] == "cached"
            mut_ok = af_result["mut"]["status"] == "cached"
            if not wt_ok:
                st.warning("AlphaFold cache not found — fetching from AlphaFold DB or manual setup required.")
                st.caption(af_result["wt"].get("how_to_generate_cache", ""))
            elif not mut_ok:
                st.info("WT structure from AlphaFold DB. Mutant not in DB — predict mutant to get ΔpLDDT.")
                plddt_arr = af_result["wt"].get("per_res_plddt", [])
                if pos and plddt_arr and 1 <= pos <= len(plddt_arr):
                    st.metric("WT pLDDT at site", round(plddt_arr[pos - 1], 1), help="Per-residue confidence at mutation site")

                mut_seq = canonical_fasta[: pos - 1] + mut + canonical_fasta[pos:]
                from alphafold_runner import _seq_hash, predict_via_esm_atlas, predict_via_local_colabfold, is_local_colabfold_available
                import os
                out_dir = os.path.join(os.path.dirname(__file__), "alphafold_out")
                mut_hash = _seq_hash(mut_seq)
                cache_dir = os.path.join(out_dir, "cache", mut_hash)
                esm_limit = 400
                colabfold_ok = is_local_colabfold_available()

                with st.expander("📌 How to get mutant structure (ΔpLDDT)", expanded=True):
                    if colabfold_ok:
                        st.markdown("#### **Option A: Predict on this server** *(any length — runs ColabFold on this machine)*")
                        st.caption("Predictions run locally on this laptop/server. Takes several minutes for long proteins.")
                        if st.button("▶ Predict mutant structure now (local ColabFold)", type="primary", key="predict_mut_local"):
                            with st.spinner("Running ColabFold locally (may take 5–30 min for long proteins)…"):
                                res = predict_via_local_colabfold(mut_seq, out_dir)
                            if res:
                                st.success("Done. Reloading…")
                                st.rerun()
                            else:
                                st.error("Local prediction failed. Check ColabFold is installed and GPU available.")
                        st.markdown("---")

                    st.markdown("#### **Option B: Predict here (cloud)** *(≤400 residues only)*")
                    st.caption("Uses ESM Atlas API. No sign-up. ~30–90 seconds.")
                    if len(mut_seq) <= esm_limit:
                        if st.button("▶ Predict via ESM Atlas (cloud)", key="predict_mut_esm"):
                            with st.spinner("Predicting via ESM Atlas (30–90 seconds)…"):
                                res = predict_via_esm_atlas(mut_seq, out_dir)
                            if res:
                                st.success("Done. Reloading…")
                                st.rerun()
                            else:
                                st.error("Prediction failed. Use another option.")
                    else:
                        st.warning(f"Sequence is {len(mut_seq)} residues. Option B supports ≤400. Use Option A (local) or C.")

                    st.markdown("---")
                    st.markdown("#### **Option C: Use an external tool** *(any length)*")
                    st.markdown("**Step 1 — Download mutant sequence**")
                    fasta_content = f">mutant_{gene}_{wt}{pos}{mut}\n{mut_seq}"
                    st.download_button(
                        label="⬇ Download mutant FASTA",
                        data=fasta_content,
                        file_name=f"mutant_{gene}_{wt}{pos}{mut}.fasta",
                        mime="text/plain",
                        key="dl_mut_fasta",
                    )
                    st.markdown("**Step 2 — Run a prediction**")
                    st.markdown("""
                    Open one of these (free, no sign-up), paste the sequence from the FASTA, and run:
                    """)
                    st.markdown("""
                    | Tool | Link |
                    |------|------|
                    | ColabFold | [colab.research.google.com/.../ColabFold](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb) |
                    | AlphaFold Server | [alphafoldserver.com](https://alphafoldserver.com) |
                    | Robetta | [robetta.bakerlab.org](https://robetta.bakerlab.org) |
                    """)
                    st.markdown("**Step 3 — Save outputs locally**")
                    st.markdown("Place these two files in the folder below:")
                    st.code(cache_dir, language=None)
                    st.markdown("""
                    | File | What to save |
                    |------|--------------|
                    | `ranked_0.pdb` | Best model PDB from the tool |
                    | `plddt.json` | `{"plddt": [92.5, 88, ...]}` — one score per residue |
                    """)
                    with st.expander("ℹ If you only have a PDB file"):
                        st.markdown("""
                        Many tools put pLDDT in the PDB B-factor column. Create `plddt.json` by listing those B-factors:
                        `{"plddt": [92.5, 88, 91, ...]}`.
                        """)
                    st.markdown("**Step 4 — Re-run**")
                    st.markdown("Click **Run Pipeline** in the sidebar. The app will load the cached mutant and show ΔpLDDT.")
            else:
                af_deltas = compute_local_structure_deltas_from_af(af_result["wt"], af_result["mut"], pos)
                st.metric("ΔpLDDT (window)", af_deltas.get("delta_mean_plddt_window"), help="AlphaFold local confidence change")
        except Exception as e:
            st.warning(f"AlphaFold: {e}")

    # Section 3: Structural impact (with AF evidence if available)
    st.header("3️⃣ Estimated Structural Impact")
    in_voltage_sensor = gene.upper() == "SCN5A" and pos and 400 < pos < 800
    struct = estimate_structural_impact(wt or "R", mut or "H", pos or 526, in_voltage_sensor, af_deltas=af_deltas)
    impact_color = {"High": "🔴", "Medium": "🟡", "Low": "🟢", "Uncertain": "⚪"}
    st.markdown(f"**Impact:** {impact_color.get(struct['impact'], '')} **{struct['impact']}**")
    for r in struct["reasons"]:
        st.markdown(f"- {r}")

    # Section 4: Tissue-specific PPIs + PPI ΔΔG table
    st.header("4️⃣ Tissue-Specific Protein Interactions & PPI ΔΔG")
    interactors = get_tissue_interactors(gene, tissue, uniprot_id=uniprot_id)
    if interactors and get_ppi_ddg_predictions and wt and mut and pos:
        ppi_rows = get_ppi_ddg_predictions(gene, wt, mut, pos, interactors, canonical_seq=canonical_fasta)
        df = pd.DataFrame(ppi_rows)
        st.dataframe(df, use_container_width=True, hide_index=True)
        st.caption("ΔΔG = NA when no complex model or residue not at interface.")
    elif interactors:
        for ip in interactors:
            label = ip.get("label", ip.get("partner", "?"))
            with st.expander(f"**{label}** — {ip.get('role', '')}"):
                st.write(f"UniProt: {ip.get('uniprot', '')}")
    else:
        st.info(f"No interactors found for {gene} in {tissue}. (STRING API or curated list)")

    # Section 5: 3D structure
    st.header("5️⃣ Interactive 3D Structure")
    pdb_choice = st.selectbox(
        "Select PDB",
        pdb_ids if pdb_ids else ["No structure"],
        index=0,
    )
    if pdb_ids and pdb_choice in pdb_ids:
        pdb_data = fetch_pdb(pdb_choice)
        if pdb_data:
            # Use mapped PDB resseq when available (same PDB as recommended)
            res_pos = pdb_resseq if (pdb_choice == pdb_rec.get("pdb_id") and pdb_resseq is not None) else pos
            col1, col2 = st.columns(2)
            with col1:
                html_wt = render_py3dmol_html(pdb_data, residue_pos=res_pos, highlight_color="royalblue", label=f"WT {wt}{pos}")
                if html_wt:
                    st.components.v1.html(html_wt, height=520, scrolling=False)
            with col2:
                html_mut = render_py3dmol_html(pdb_data, residue_pos=res_pos, highlight_color="red", label=f"Mut {wt}{pos}→{mut}")
                if html_mut:
                    st.components.v1.html(html_mut, height=520, scrolling=False)
        else:
            st.warning("Could not fetch PDB.")
    else:
        st.info("No PDB for this protein. Add to config or use AlphaFold DB.")

    st.header("6️⃣ PPI Binding Affinity (mCSM-PPI2)")
    st.link_button("Open mCSM-PPI2", "https://biosig.lab.uq.edu.au/mcsm_ppi2/submit_prediction")

st.success("✅ Pipeline complete.")
