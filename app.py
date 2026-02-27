"""
Mutation PPI Prediction Pipeline — Interactive Streamlit App
Input: gene, mutation, tissue → Variant QC, pathogenicity, structural impact, PPI (interface-aware)
"""

import streamlit as st
import pandas as pd
from config import GENE_UNIPROT
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
    .main-header { font-size: 2rem; font-weight: 700; color: #1E88E5; margin-bottom: 0.5rem; }
    .sub-header { font-size: 1.1rem; color: #546E7A; margin-bottom: 2rem; }
</style>
""", unsafe_allow_html=True)

st.markdown('<p class="main-header">🧬 Mutation PPI Prediction Pipeline</p>', unsafe_allow_html=True)
st.markdown(
    '<p class="sub-header">Variant QC → pathogenicity, structural impact, PPI (interface-aware)</p>',
    unsafe_allow_html=True,
)

# Sidebar — Input
with st.sidebar:
    st.header("📥 Input")
    gene = st.text_input("Gene symbol", value="SCN5A", help="e.g., SCN5A, MYH7, KCNQ1")
    mutation_input = st.text_input(
        "Mutation",
        value="c.1577G>A, p.R526H",
        help="Formats: c.1577G>A, p.R526H, or R526H",
    )
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
parsed = parse_mutation(mutation_input, gene)
pos = parsed.get("position") or parsed.get("cds_pos")
wt, mut = parsed.get("wt_aa"), parsed.get("mut_aa")

if not pos and not parsed.get("cds_pos"):
    st.warning("Could not parse mutation. Use formats: p.R526H, R526H, or c.1577G>A")
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

# PDB coverage
pdb_rec = get_recommended_pdb(gene, pos)
has_res = pdb_rec.get("has_residue_coordinates", True)
pdb_ids = pdb_rec.get("pdb_ids", [])

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
                st.info("WT structure from AlphaFold DB. Mutant not in DB — ΔpLDDT requires local AlphaFold run.")
                plddt_arr = af_result["wt"].get("per_res_plddt", [])
                if pos and plddt_arr and 1 <= pos <= len(plddt_arr):
                    st.metric("WT pLDDT at site", round(plddt_arr[pos - 1], 1), help="Per-residue confidence at mutation site")
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
    interactors = get_tissue_interactors(gene, tissue)
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
        st.info(f"No predefined interactors for {gene} in {tissue}.")

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
            col1, col2 = st.columns(2)
            with col1:
                html_wt = render_py3dmol_html(pdb_data, residue_pos=pos, highlight_color="royalblue", label=f"WT {wt}{pos}")
                if html_wt:
                    st.components.v1.html(html_wt, height=520, scrolling=False)
            with col2:
                html_mut = render_py3dmol_html(pdb_data, residue_pos=pos, highlight_color="red", label=f"Mut {wt}{pos}→{mut}")
                if html_mut:
                    st.components.v1.html(html_mut, height=520, scrolling=False)
        else:
            st.warning("Could not fetch PDB.")
    else:
        st.info("No PDB for this protein. Add to config or use AlphaFold DB.")

    st.header("6️⃣ PPI Binding Affinity (mCSM-PPI2)")
    st.link_button("Open mCSM-PPI2", "https://biosig.lab.uq.edu.au/mcsm_ppi2/submit_prediction")

st.success("✅ Pipeline complete.")
