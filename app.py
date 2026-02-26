"""
Mutation PPI Prediction Pipeline — Interactive Streamlit App
Input: gene, mutation, tissue → outputs all predictions + high-res 3D visualization
"""

import streamlit as st
from config import GENE_UNIPROT, PROTEIN_PDB
from predictors import (
    parse_mutation,
    get_alphamissense_prediction,
    get_tissue_interactors,
    get_recommended_pdb,
    estimate_structural_impact,
    resolve_uniprot,
    get_ppi_ddg_predictions,
)
from visualization import fetch_pdb, render_py3dmol_html

# Page config
st.set_page_config(
    page_title="Mutation PPI Pipeline",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Custom CSS for a clean look
st.markdown("""
<style>
    .main-header { font-size: 2rem; font-weight: 700; color: #1E88E5; margin-bottom: 0.5rem; }
    .sub-header { font-size: 1.1rem; color: #546E7A; margin-bottom: 2rem; }
    .metric-card { padding: 1rem; border-radius: 8px; background: #F5F5F5; margin: 0.5rem 0; }
    .result-box { border-left: 4px solid #1E88E5; padding: 1rem; margin: 1rem 0; }
    .warning-box { border-left: 4px solid #FF9800; padding: 1rem; margin: 1rem 0; }
</style>
""", unsafe_allow_html=True)

# Header
st.markdown('<p class="main-header">🧬 Mutation PPI Prediction Pipeline</p>', unsafe_allow_html=True)
st.markdown(
    '<p class="sub-header">Predict pathogenicity, structural impact, and protein-protein interactions for missense mutations in tissue context</p>',
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
    run_button = st.button("🚀 Run Pipeline", type="primary", use_container_width=True)

# Parse mutation
parsed = parse_mutation(mutation_input)
if not parsed.get("position") and not parsed.get("cds_pos"):
    st.warning("Could not parse mutation. Use formats: p.R526H, R526H, or c.1577G>A")
    st.stop()

# Resolve UniProt
uniprot_id = GENE_UNIPROT.get(gene.upper()) or resolve_uniprot(gene)
if not uniprot_id:
    st.error(f"Could not resolve UniProt ID for gene {gene}")

if run_button or st.session_state.get("results_ready"):
    st.session_state["results_ready"] = True

    with st.spinner("Running predictions..."):
        # --- Section 1: Parsed mutation summary ---
        st.header("1️⃣ Parsed Mutation")
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Gene", gene)
        with col2:
            st.metric("Wild-type", parsed.get("wt_aa", "—"))
        with col3:
            st.metric("Position", parsed.get("position", parsed.get("cds_pos", "—")))
        with col4:
            st.metric("Mutant", parsed.get("mut_aa", "—"))
        pos = parsed.get("position") or parsed.get("cds_pos")
        wt, mut = parsed.get("wt_aa"), parsed.get("mut_aa")
        if not wt or not mut:
            st.info("Amino acids inferred from position. Add p.R526H format for full parsing.")

        # --- Section 2: AlphaMissense ---
        st.header("2️⃣ AlphaMissense Pathogenicity")
        if pos and wt and mut:
            am_result = get_alphamissense_prediction(gene, pos, wt, mut)
            if "error" in am_result:
                st.warning(f"AlphaMissense API: {am_result['error']}")
                st.info("You can check manually: https://alphamissense.hegelab.org/")
            else:
                score = am_result.get("pathogenicity")
                if score is not None:
                    if score > 0.5:
                        st.error(f"Pathogenic likelihood: **{score:.2%}**")
                    else:
                        st.success(f"Pathogenic likelihood: **{score:.2%}**")
                st.json(am_result.get("raw", am_result))
        else:
            st.info("Provide mutation in p.R526H format for AlphaMissense.")

        # --- Section 3: Structural impact heuristic ---
        st.header("3️⃣ Estimated Structural Impact")
        in_voltage_sensor = gene.upper() == "SCN5A" and pos and 400 < pos < 800
        struct = estimate_structural_impact(wt or "R", mut or "H", pos or 526, in_voltage_sensor)
        impact_color = {"High": "🔴", "Medium": "🟡", "Low": "🟢"}
        st.markdown(f"**Impact level:** {impact_color.get(struct['impact'], '')} **{struct['impact']}**")
        for r in struct["reasons"]:
            st.markdown(f"- {r}")
        if in_voltage_sensor:
            st.info("Position likely in DII voltage-sensor region — charge changes can severely affect gating.")

        # --- Section 4: Tissue-specific PPIs + PPI ΔΔG (auto-calculated on selection) ---
        st.header("4️⃣ Tissue-Specific Protein Interactions & PPI ΔΔG")
        interactors = get_tissue_interactors(gene, tissue)
        if interactors:
            partner_options = [f"{ip['partner']} — {ip['role']}" for ip in interactors]
            selected_labels = st.multiselect(
                "Select interacting proteins (ΔΔG calculated immediately)",
                options=partner_options,
                default=partner_options,
                key="ppi_select",
            )
            selected_indices = [i for i, l in enumerate(partner_options) if l in selected_labels]
            selected_interactors = [interactors[i] for i in selected_indices]
            if selected_interactors and wt and mut and pos:
                ddg_results = get_ppi_ddg_predictions(gene, selected_interactors, wt, pos, mut)
                for r in ddg_results:
                    ddg = r["mutant_ddg"]
                    ddg_msg = f"ΔΔG ≈ **{ddg} kcal/mol**" + (" (weaker binding)" if ddg < 0 else " (similar binding)")
                    with st.expander(f"**{r['partner']}** — {r['role']} · {ddg_msg}"):
                        st.write(f"UniProt: {r['uniprot']}")
                        st.caption(f"Heuristic method · negative ΔΔG = weaker mutant binding")
            elif selected_interactors:
                for ip in selected_interactors:
                    with st.expander(f"**{ip['partner']}** — {ip['role']}"):
                        st.write(f"UniProt: {ip['uniprot']}")
                st.info("Provide mutation in p.R526H format for ΔΔG calculation.")
        else:
            st.info(f"No predefined interactors for {gene} in {tissue}. Add to config.CARDIAC_MYOCYTE_INTERACTORS.")

        # --- Section 5: High-resolution 3D structure ---
        st.header("5️⃣ Interactive 3D Structure")
        pdb_ids = get_recommended_pdb(gene)
        pdb_choice = st.selectbox(
            "Select PDB structure",
            pdb_ids if pdb_ids else ["No structure in database"],
            index=0,
        )
        if pdb_ids and pdb_choice in pdb_ids:
            pdb_data = fetch_pdb(pdb_choice)
            if pdb_data:
                st.markdown("**Rotate, zoom, pan** — compare wild-type (blue) vs mutant (red) at position")
                wt, mut = parsed.get("wt_aa", "—"), parsed.get("mut_aa", "—")
                col1, col2 = st.columns(2)
                with col1:
                    st.markdown(f"**Wild-type** — {wt}{pos}")
                    html_wt = render_py3dmol_html(pdb_data, residue_pos=pos, highlight_color="royalblue", label=f"Wild-type {wt}{pos}")
                    if html_wt:
                        st.components.v1.html(html_wt, height=520, scrolling=False)
                with col2:
                    st.markdown(f"**Mutant** — {wt}{pos}→{mut}")
                    html_mut = render_py3dmol_html(pdb_data, residue_pos=pos, highlight_color="red", label=f"Mutant {wt}{pos}→{mut}")
                    if html_mut:
                        st.components.v1.html(html_mut, height=520, scrolling=False)
                if not html_wt and not html_mut:
                    st.warning("3D viewer unavailable. Install: pip install py3Dmol")
            else:
                st.warning("Could not fetch PDB structure.")
        else:
            st.info(f"No predefined PDB for {gene}. Add to config.PROTEIN_PDB or use AlphaFold DB.")

        # --- Section 6: Next steps / mCSM-PPI2 ---
        st.header("6️⃣ PPI Binding Affinity (mCSM-PPI2)")
        st.markdown("""
        For interface residues, use **mCSM-PPI2** to predict ΔΔG:
        - [mCSM-PPI2 Web Server](https://biosig.lab.uq.edu.au/mcsm_ppi2/)
        - Upload your PDB complex and mutation (e.g., `A R 526 H`)
        - R526 is in the voltage sensor; PPI tools apply if you identify an interface involving this residue.
        """)
        st.link_button("Open mCSM-PPI2", "https://biosig.lab.uq.edu.au/mcsm_ppi2/submit_prediction")

    st.success("✅ Pipeline complete.")

else:
    st.info("👈 Enter mutation and tissue, then click **Run Pipeline** to start.")
