# Run the App Locally So Users Submit to Your Laptop

This guide lets you run the Mutation PPI Pipeline on your laptop so **predictions run on your machine** when users submit requests. Users access your app via a public URL (ngrok).

---

## Prerequisites

1. **Python 3.8+** with the project dependencies
2. **ColabFold** installed (for local structure prediction — any length)
3. **ngrok** (free) to expose your laptop to the internet

---

## Step 1: Install ColabFold

ColabFold runs AlphaFold2 locally. Install via Docker (recommended) or conda.

### Option A: Docker (recommended)

```bash
# Pull ColabFold Docker image
docker pull ghcr.io/sokrypton/colabfold:1.5.3-cuda12.2.2

# Test it
docker run --rm ghcr.io/sokrypton/colabfold:1.5.3-cuda12.2.2 colabfold_batch --help
```

**Use Docker mode in the app** — before running Streamlit, set:

```bash
# Windows (PowerShell)
$env:COLABFOLD_DOCKER="1"

# Linux/Mac
export COLABFOLD_DOCKER=1
```

Then run `streamlit run app.py`. Option A (local prediction) will use Docker automatically.

### Option B: Conda

```bash
conda create -n colabfold python=3.10
conda activate colabfold
pip install colabfold[alphafold]
```

---

## Step 2: Run the Streamlit App

```bash
cd C:\Users\chime\mutation-ppi-pipeline
pip install -r requirements.txt

# If using Docker for ColabFold:
# Windows: $env:COLABFOLD_DOCKER="1"
# Linux/Mac: export COLABFOLD_DOCKER=1

python -m streamlit run app.py --server.port 8501
```

The app will be available at **http://localhost:8501** on your machine.

---

## Step 3: Expose with ngrok (so users can reach your laptop)

1. **Sign up** (free): https://dashboard.ngrok.com/signup
2. **Get your authtoken**: https://dashboard.ngrok.com/get-started/your-authtoken
3. **Add authtoken** (one-time):

```bash
ngrok config add-authtoken YOUR_TOKEN
```

4. **Run ngrok**:

```bash
ngrok http 8501
```

You’ll get a URL like `https://xxxx.ngrok-free.dev`. **Share this URL with users.**

**Optional:** Update `APP_PUBLIC_URL` in `config.py` with your ngrok URL — the app will then display the link prominently so users can easily access it.

---

## Step 4: Keep It Running

- Keep your laptop on and connected to the internet
- Leave both `streamlit run app.py` and `ngrok http 8501` running
- Free ngrok URLs change each time you restart ngrok (paid plans offer fixed URLs)

---

## Summary

| Step | Command / Action |
|------|------------------|
| 1 | Install ColabFold (Docker or conda) |
| 2 | `python -m streamlit run app.py --server.port 8501` |
| 3 | `ngrok http 8501` |
| 4 | Share the ngrok URL with users |

Users submit requests → your app receives them → ColabFold runs on your laptop → results are shown to the user.
