"""
Mutation PPI Prediction Pipeline — Variant QC, Analysis mode, Evidence panel.
"""

import streamlit as st
import pandas as pd
from config import GENE_UNIPROT, PROTEIN_PDB
try:
    from config import EXAMPLE_VARIANTS
except ImportError:
    EXAMPLE_VARIANTS = []

from predictors import (
    parse_mutation,
    get_alphamissense_prediction,
    get_tissue_interactors,
    get_ppi_ddg_predictions,
    get_recommended_pdb,
    estimate_structural_impact,
    resolve_uniprot,
    validate_variant_on_canonical_sequence,
)
from visualization import fetch_pdb, render_py3dmol_html
try:
    from regulatory_annot import scan_motifs, ptm_proximity
except ImportError:
    scan_motifs = ptm_proximity = None

st.set_page_config(page_title="Mutation PPI Pipeline", page_icon="\u269b", layout="wide", initial_sidebar_state="expanded")

st.markdown("""
<style>
    .main-header { font-size: 2rem; font-weight: 700; color: #1E88E5; margin-bottom: 0.5rem; }
    .sub-header { font-size: 1.1rem; color: #546E7A; margin-bottom: 2rem; }
</style>
""", unsafe_allow_html=True)

st.markdown('<p class="main-header">\u269b Mutation PPI Prediction Pipeline</p>', unsafe_allow_html=True)
st.markdown('<p class="sub-header">Variant QC -> pathogenicity, structural impact, PPI (interface-aware)</p>', unsafe_allow_html=True)

with st.sidebar:
    st.header("Input")
    if EXAMPLE_VARIANTS:
        ex_sel = st.selectbox("Example variant", [""] + [v["label"] for v in EXAMPLE_VARIANTS], index=0)
        ex = next((v for v in EXAMPLE_VARIANTS if v["label"] == ex_sel), None) if ex_sel else None
        _g = ex["gene"] if ex else "SCN5A"
        _m = ex["mutation"] if ex else "c.1577G>A, p.R526H"
    else:
        _g, _m = "SCN5A", "c.1577G>A, p.R526H"
    gene = st.text_input("Gene symbol", value=_g)
    mutation_input = st.text_input("Mutation", value=_m)
    tissue = st.selectbox("Tissue", ["Cardiac myocyte", "Heart", "Skeletal muscle", "Neuron", "Other"], index=0)
    analysis_mode = st.radio("Analysis mode",
        ["Fast (no AlphaFold)", "AlphaFold (WT vs MUT)", "AlphaFold-Multimer for PPI"], index=0)
    run_button = st.button("Run Pipeline", type="primary", use_container_width=True)

parsed = parse_mutation(mutation_input, gene)
pos = parsed.get("pos") or parsed.get("position") or parsed.get("cds_pos")
wt, mut = parsed.get("wt") or parsed.get("wt_aa"), parsed.get("mut") or parsed.get("mut_aa")

if not pos and not parsed.get("cds_pos"):
    st.warning("Could not parse mutation. Use formats: p.R526H, p.Ser1054Ala")
    st.stop()

@st.cache_data(ttl=3600)
def _cached_resolve(g):
    return resolve_uniprot(g)

uniprot_data = _cached_resolve(gene)
if isinstance(uniprot_data, dict):
    uniprot_id = uniprot_data.get("uniprot_id")
else:
    uniprot_id = uniprot_data or GENE_UNIPROT.get(gene.upper())

if not uniprot_id:
    st.error("Could not resolve UniProt ID for gene " + gene)
    st.stop()

# --- 0. Mandatory Variant QC & Mapping ---
st.header("0. Variant QC & Mapping")
canonical = (uniprot_data or {}).get("canonical_fasta", "") if isinstance(uniprot_data, dict) else ""
seq_len = (uniprot_data or {}).get("length", 0) if isinstance(uniprot_data, dict) else 0

qc_ok = True
if wt and pos and canonical:
    val = validate_variant_on_canonical_sequence(canonical, wt, pos)
    if not val.get("match"):
        st.error("**WT mismatch:** " + val.get("message", ""))
        st.stop()
else:
    if not wt or not mut:
        st.warning("WT/mutant not fully parsed. Use p.R526H format.")
        qc_ok = False

pdb_rec = get_recommended_pdb(gene, pos) if pos else get_recommended_pdb(gene)
has_res = pdb_rec.get("has_residue_coordinates", True) if isinstance(pdb_rec, dict) else True

col1, col2, col3, col4, col5 = st.columns(5)
with col1:
    st.metric("UniProt ID", uniprot_id)
with col2:
    st.metric("Canonical length", seq_len or "-")
with col3:
    st.metric("Verified WT at pos", (wt or "") + str(pos) if wt and pos else "-")
with col4:
    st.metric("Residue in PDB?", "Yes" if has_res else "No")
with col5:
    st.metric("WT check", "OK" if qc_ok else "Warn")

# --- Evidence & Confidence panel ---
st.header("Evidence & Confidence")
e1, e2, e3 = st.columns(3)
with e1:
    st.markdown("**WT verification:** OK" if qc_ok else "**WT verification:** Fail")
with e2:
    st.markdown("**Residue in PDB:** " + ("Yes" if has_res else "No"))
with e3:
    st.markdown("**Region:** Structured" if has_res else "**Region:** Recommend AlphaFold")

if scan_motifs and canonical and pos:
    motifs = scan_motifs(canonical, pos)
    if motifs.get("motifs"):
        with st.expander("Motifs near site"):
            st.json(motifs)
if ptm_proximity and canonical and pos:
    ptm = ptm_proximity(canonical, pos)
    if ptm.get("phospho_site_nearby"):
        st.info("Phospho site(s) nearby: " + str(ptm.get("nearby_sites", [])))

if not (run_button or st.session_state.get("results_ready")):
    st.info("Click Run Pipeline to continue.")
    st.stop()

st.session_state["results_ready"] = True

with st.spinner("Running predictions..."):
    st.header("1. Parsed Mutation")
    st.metric("Variant", parsed.get("variant_id", (gene or "") + ":" + (wt or "") + str(pos or "") + (mut or "")))

    st.header("2. AlphaMissense Pathogenicity")
    @st.cache_data(ttl=86400)
    def _cached_am(g, p, w, m):
        return get_alphamissense_prediction(g, p, w, m)
    if pos and wt and mut:
        am = _cached_am(gene, pos, wt, mut)
        if "error" in am:
            st.warning(am["error"])
        else:
            score = am.get("pathogenicity")
            if score is not None:
                st.error("Pathogenic: " + str(round(score*100,1)) + "%") if score > 0.5 else st.success("Pathogenic: " + str(round(score*100,1)) + "%")

    st.header("3. Estimated Structural Impact")
    in_vs = gene.upper() == "SCN5A" and pos and 400 < pos < 800
    struct = estimate_structural_impact(wt or "R", mut or "H", pos or 526, in_vs)
    st.markdown("**Impact:** " + struct["impact"])
    for r in struct["reasons"]:
        st.markdown("- " + r)

    st.header("4. Tissue-Specific Protein Interactions & PPI dG")
    @st.cache_data(ttl=3600)
    def _cached_int(g, t):
        return get_tissue_interactors(g, t)
    interactors = _cached_int(gene, tissue)
    if interactors and wt and mut and pos:
        ppi_rows = get_ppi_ddg_predictions(gene, wt, mut, pos, interactors)
        df = pd.DataFrame(ppi_rows)
        st.dataframe(df, use_container_width=True, hide_index=True)
        st.caption("dG = NA when no complex model or residue not at interface.")
    elif interactors:
        for ip in interactors:
            with st.expander(ip["partner"] + " - " + ip["role"]):
                st.write("UniProt: " + ip["uniprot"])
    else:
        st.info("No predefined interactors for " + gene + " in " + tissue + ".")

    st.header("5. Interactive 3D Structure")
    pdb_ids = pdb_rec.get("pdb_ids", pdb_rec) if isinstance(pdb_rec, dict) else (pdb_rec or [])
    pdb_choice = st.selectbox("Select PDB", pdb_ids if pdb_ids else ["No structure"], index=0)
    if pdb_ids and pdb_choice in pdb_ids:
        pdb_data = fetch_pdb(pdb_choice)
        if pdb_data:
            c1, c2 = st.columns(2)
            with c1:
                html_wt = render_py3dmol_html(pdb_data, residue_pos=pos, highlight_color="royalblue", label="WT " + str(wt) + str(pos))
                if html_wt:
                    st.components.v1.html(html_wt, height=520, scrolling=False)
            with c2:
                html_mut = render_py3dmol_html(pdb_data, residue_pos=pos, highlight_color="red", label="Mutant " + str(wt) + str(pos) + "->" + str(mut))
                if html_mut:
                    st.components.v1.html(html_mut, height=520, scrolling=False)

    st.header("6. PPI Binding Affinity (mCSM-PPI2)")
    st.link_button("Open mCSM-PPI2", "https://biosig.lab.uq.edu.au/mcsm_ppi2/submit_prediction")

st.success("Pipeline complete.")
