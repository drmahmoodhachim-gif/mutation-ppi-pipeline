# Mutation PPI Prediction Pipeline

An **interactive pipeline** that takes a **mutation** (e.g., SCN5A c.1577G>A, p.R526H) and **tissue of interest** (e.g., cardiac myocyte) and automatically runs:

1. **AlphaMissense** — pathogenicity prediction  
2. **Structural impact** — physicochemical/charge change heuristic  
3. **Tissue-specific protein interactions** — curated PPI partners  
4. **High-resolution 3D structure** — interactive PyMOL-like viewer with mutant residue highlighted  
5. **PPI binding affinity** — link to mCSM-PPI2 for interface mutations  

## Quick start (local)

```bash
git clone https://github.com/drmahmoodhachim-gif/mutation-ppi-pipeline.git
cd mutation-ppi-pipeline
pip install -r requirements.txt
streamlit run app.py
```

Then open the URL (usually http://localhost:8501) in your browser.

## Deploy live (Streamlit Community Cloud)

1. **Push to GitHub** (see below)
2. Go to [share.streamlit.io](https://share.streamlit.io)
3. Sign in with GitHub
4. Click **"New app"** → **"Deploy an existing app"**
5. Set:
   - **Repository:** `YOUR_USERNAME/mutation-ppi-pipeline`
   - **Branch:** `main`
   - **Main file path:** `app.py`
6. Click **Deploy**

Your app will be live at `https://YOUR_APP_NAME.streamlit.app`

## Input formats

- **Mutation:** `c.1577G>A`, `p.R526H`, `R526H`, or combined `c.1577G>A, p.R526H`
- **Gene:** SCN5A, MYH7, KCNQ1, etc. (see `config.py` for supported genes)
- **Tissue:** Cardiac myocyte, Heart, Skeletal muscle, Neuron, Other  

## Project structure

```
mutation-ppi-pipeline/
├── app.py           # Streamlit app
├── config.py        # Gene–UniProt, tissue PPIs, PDB IDs
├── predictors.py    # AlphaMissense, structural impact, PPIs
├── visualization.py # 3D structure fetch and viewer
├── requirements.txt
└── README.md
```

## Extending the pipeline

- **Add genes:** Edit `GENE_UNIPROT` in `config.py`
- **Add tissue PPIs:** Edit `CARDIAC_MYOCYTE_INTERACTORS` or add new tissue dicts
- **Add PDB structures:** Edit `PROTEIN_PDB` in `config.py`
- **Local ColabFold:** Uncomment in `requirements.txt` and add a prediction step in `predictors.py`

## Tools used

| Tool | Purpose |
|------|---------|
| AlphaMissense (hegelab API) | Pathogenicity prediction |
| RCSB PDB | High-resolution structures |
| stmol + py3Dmol | Interactive 3D viewer |
| mCSM-PPI2 | PPI ΔΔG (manual link) |

## License

MIT
