# How to Run the Mutation PPI Pipeline

## Step 1: Open terminal in the project folder
```powershell
cd C:\Users\chime\mutation-ppi-pipeline
```

## Step 2: Create a virtual environment (recommended)
```powershell
python -m venv venv
.\venv\Scripts\Activate.ps1
```

## Step 3: Install dependencies
```powershell
pip install -r requirements.txt
```

## Step 4: Run the app
```powershell
python -m streamlit run app.py
```

## Step 5: Open in browser
- The terminal will show: **Local URL: http://localhost:8501** (or 8502, 8503, etc.)
- Click the link or open that URL in your browser

## Step 6: Use the app
1. Gene: **SCN5A** (default)
2. Mutation: **c.1577G>A, p.R526H** (or p.Ser1054Ala, etc.)
3. Tissue: **Cardiac myocyte**
4. Click **Run Pipeline**

---

## Quick run (no venv)
```powershell
cd C:\Users\chime\mutation-ppi-pipeline
pip install streamlit py3Dmol requests pandas
python -m streamlit run app.py
```

## If port is in use
```powershell
python -m streamlit run app.py --server.port 8507
```
