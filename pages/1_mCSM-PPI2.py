"""
mCSM-PPI2 — Predict ΔΔG for interface residues.
Embedded web tool for protein-protein interaction binding affinity.
"""

import streamlit as st

st.header("🔬 PPI Binding Affinity (mCSM-PPI2)")
st.markdown("""
For **interface residues**, use mCSM-PPI2 to predict the effects of mutations on protein–protein interaction binding affinity (ΔΔG):

- Upload your PDB complex and mutation (e.g., chain `A`, wild-type `R`, position `526`, mutant `H`)
- R526 in SCN5A is in the voltage sensor; PPI tools apply if you identify an interface involving this residue
""")

# Embedded mCSM-PPI2 submission form
MCSM_URL = "https://biosig.lab.uq.edu.au/mcsm_ppi2/submit_prediction"
iframe_html = f"""
<div style="width: 100%; height: 900px;">
    <iframe 
        src="{MCSM_URL}" 
        width="100%" 
        height="900" 
        frameborder="0" 
        allowfullscreen
        style="border: 1px solid #ddd; border-radius: 8px;">
    </iframe>
</div>
"""
st.components.v1.html(iframe_html, height=920, scrolling=False)

st.markdown("---")
st.markdown("**If the tool above does not load**, open it in a new tab:")
st.link_button("Open mCSM-PPI2 in new tab", MCSM_URL)
