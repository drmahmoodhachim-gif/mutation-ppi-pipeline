"""AlphaFold runner - pluggable (local/remote/cached)."""

import hashlib
import json
import os
import shutil
import subprocess
import tempfile
from typing import Optional

try:
    import requests
except ImportError:
    requests = None

AF_DB_API = "https://alphafold.ebi.ac.uk/api/prediction"
ESM_ATLAS_API = "https://api.esmatlas.com/foldSequence/v1/pdb/"
ESM_ATLAS_MAX_LEN = 400


def _seq_hash(seq: str) -> str:
    return hashlib.sha256(seq.encode()).hexdigest()[:16]


HOW_TO_GENERATE = (
    "Save AlphaFold outputs as: {cache_dir}/ranked_0.pdb and {cache_dir}/plddt.json "
    "(JSON with 'plddt' or 'confidence' array of per-residue scores)."
)


def _parse_plddt_from_ebi(data: dict) -> list:
    """Extract per-residue pLDDT from EBI AlphaFold confidence JSON (various formats)."""
    # Format 1: {"confidenceScore": [92.5, 88, ...], "residueNumber": [1,2,3,...]} (EBI v6)
    scores = data.get("confidenceScore")
    if isinstance(scores, (list, tuple)) and scores:
        return [float(x) for x in scores]
    # Format 2: {"residues": [{"residueNumber": 1, "confidenceScore": 92.5}, ...]}
    residues = data.get("residues", [])
    if residues and isinstance(residues[0], dict):
        out = [0.0] * max(len(residues), max(r.get("residueNumber", i + 1) for i, r in enumerate(residues)))
        for r in residues:
            idx = r.get("residueNumber", 0) - 1
            if 0 <= idx < len(out):
                out[idx] = float(r.get("confidenceScore", r.get("plddt", 0)))
        return [x for x in out if x > 0] if any(x > 0 for x in out) else [float(r.get("confidenceScore", r.get("plddt", 0))) for r in residues]
    # Format 3: {"plddt": [...]} or {"confidence": [...]}
    for key in ("plddt", "confidence"):
        arr = data.get(key)
        if isinstance(arr, (list, tuple)) and arr:
            return [float(x) for x in arr]
    return []


def fetch_from_alphafold_db(uniprot_id: str, sequence: str, out_dir: str) -> Optional[str]:
    """
    Fetch WT structure and pLDDT from AlphaFold Database (EBI), save to cache.
    Returns cache_dir path if successful, None otherwise. Uses canonical UniProt ID (no isoform).
    """
    if not requests or not uniprot_id:
        return None
    uid = uniprot_id.split("-")[0] if "-" in uniprot_id else uniprot_id  # Q14524-4 -> Q14524
    try:
        r = requests.get(f"{AF_DB_API}/{uid}", timeout=30)
        if not r.ok:
            return None
        models = r.json()
        if not models or not isinstance(models, list):
            return None
        # Prefer model matching canonical (Q14524-F1, full length)
        m = next((x for x in models if x.get("uniprotAccession") == uid and x.get("sequenceEnd", 0) >= len(sequence)), models[0])
        pdb_url = m.get("pdbUrl")
        conf_url = m.get("plddtDocUrl") or m.get("confidenceUrl")
        if not pdb_url:
            return None
        cache_key = _seq_hash(sequence)
        cache_dir = os.path.join(out_dir, "cache", cache_key)
        os.makedirs(cache_dir, exist_ok=True)
        pdb_path = os.path.join(cache_dir, "ranked_0.pdb")
        plddt_path = os.path.join(cache_dir, "plddt.json")
        # Download PDB
        rp = requests.get(pdb_url, timeout=60)
        if rp.ok:
            with open(pdb_path, "wb") as f:
                f.write(rp.content)
        else:
            return None
        # Download and convert confidence JSON
        per_res = []
        if conf_url:
            rc = requests.get(conf_url, timeout=30)
            if rc.ok:
                try:
                    raw = rc.json()
                    per_res = _parse_plddt_from_ebi(raw)
                    if not per_res and isinstance(raw, list):
                        per_res = [float(x.get("confidenceScore", x.get("plddt", 0))) for x in raw if isinstance(x, dict)]
                    if not per_res and isinstance(raw, dict):
                        per_res = raw.get("plddt", raw.get("confidence", []))
                except Exception:
                    pass
        # Fallback: extract pLDDT from B-factor in PDB
        if not per_res and os.path.isfile(pdb_path):
            per_res = []
            with open(pdb_path, "r") as f:
                prev = None
                for line in f:
                    if line.startswith(("ATOM", "HETATM")) and len(line) >= 60:
                        try:
                            resi = line[22:26].strip()
                            bfac = float(line[60:66])
                            if (resi, line[21:22]) != prev:
                                prev = (resi, line[21:22])
                                per_res.append(bfac)
                        except (ValueError, IndexError):
                            pass
        with open(plddt_path, "w") as f:
            json.dump({"plddt": per_res}, f)
        return cache_dir
    except Exception:
        return None


def _extract_plddt_from_pdb(pdb_path: str) -> list:
    """Extract per-residue pLDDT from B-factor column in PDB (ESMFold/AlphaFold format)."""
    per_res = []
    prev = None
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
                try:
                    resi = line[22:26].strip()
                    bfac = float(line[60:66])
                    if (resi, line[21:22]) != prev:
                        prev = (resi, line[21:22])
                        per_res.append(bfac)
                except (ValueError, IndexError):
                    pass
    return per_res


def predict_via_esm_atlas(sequence: str, out_dir: str) -> Optional[dict]:
    """
    Predict structure via ESM Atlas API (no auth). Sequences up to 400 residues.
    Saves ranked_0.pdb and plddt.json to cache, returns same dict as run_alphafold_single when cached.
    Returns None if API fails or sequence too long.
    """
    if not requests or len(sequence) > ESM_ATLAS_MAX_LEN:
        return None
    try:
        r = requests.post(ESM_ATLAS_API, data=sequence.encode(), timeout=120)
        if not r.ok:
            return None
        cache_key = _seq_hash(sequence)
        cache_dir = os.path.join(out_dir, "cache", cache_key)
        os.makedirs(cache_dir, exist_ok=True)
        pdb_path = os.path.join(cache_dir, "ranked_0.pdb")
        plddt_path = os.path.join(cache_dir, "plddt.json")
        with open(pdb_path, "wb") as f:
            f.write(r.content)
        per_res = _extract_plddt_from_pdb(pdb_path)
        if not per_res:
            return None
        mean_plddt = float(sum(per_res) / len(per_res))
        with open(plddt_path, "w") as f:
            json.dump({"plddt": per_res}, f)
        return {
            "status": "cached",
            "pdb_path": pdb_path,
            "plddt_path": plddt_path,
            "mean_plddt": mean_plddt,
            "per_res_plddt": per_res,
        }
    except Exception:
        return None


def _find_ranked_pdb(dir_path: str) -> Optional[str]:
    """Find ranked_0.pdb in ColabFold output (may be in subdirs)."""
    for root, _, files in os.walk(dir_path):
        if "ranked_0.pdb" in files:
            return os.path.join(root, "ranked_0.pdb")
    return None


def predict_via_local_colabfold(sequence: str, out_dir: str, colabfold_cmd: str = None) -> Optional[dict]:
    """
    Run ColabFold locally (any length). Requires colabfold_batch or Docker.
    Set env COLABFOLD_CMD for custom command, or COLABFOLD_DOCKER=1 to use Docker.
    Returns result dict on success, None on failure.
    """
    if not sequence:
        return None
    cmd_name = colabfold_cmd or os.environ.get("COLABFOLD_CMD", "colabfold_batch")
    use_docker = os.environ.get("COLABFOLD_DOCKER", "").lower() in ("1", "true", "yes")

    cache_key = _seq_hash(sequence)
    cache_dir = os.path.join(out_dir, "cache", cache_key)
    os.makedirs(cache_dir, exist_ok=True)

    try:
        with tempfile.TemporaryDirectory() as tmp:
            fasta_path = os.path.join(tmp, "query.fasta")
            with open(fasta_path, "w") as f:
                f.write(f">query\n{sequence}\n")

            if use_docker:
                # Docker: mount tmp to /data, run colabfold_batch inside container
                img = os.environ.get("COLABFOLD_IMAGE", "ghcr.io/sokrypton/colabfold:1.5.3-cuda12.2.2")
                tmp_abs = os.path.abspath(tmp)
                # Docker -v format: host_path:container_path
                cmd = [
                    "docker", "run", "--rm",
                    "-v", f"{tmp_abs}:/data", "-w", "/data",
                    img, "colabfold_batch", "/data/query.fasta", "/data",
                ]
            else:
                cmd = [cmd_name, fasta_path, tmp]

            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600, cwd=tmp)
            if proc.returncode != 0:
                return None

            pdb_src = _find_ranked_pdb(tmp)
            if not pdb_src or not os.path.isfile(pdb_src):
                return None

            pdb_dst = os.path.join(cache_dir, "ranked_0.pdb")
            shutil.copy2(pdb_src, pdb_dst)
            per_res = _extract_plddt_from_pdb(pdb_dst)
            if not per_res:
                return None
            mean_plddt = float(sum(per_res) / len(per_res))
            plddt_path = os.path.join(cache_dir, "plddt.json")
            with open(plddt_path, "w") as f:
                json.dump({"plddt": per_res}, f)

            return {
                "status": "cached",
                "pdb_path": pdb_dst,
                "plddt_path": plddt_path,
                "mean_plddt": mean_plddt,
                "per_res_plddt": per_res,
            }
    except (subprocess.TimeoutExpired, FileNotFoundError, Exception):
        return None


def is_local_colabfold_available(colabfold_cmd: str = None) -> bool:
    """Check if ColabFold is available (colabfold_batch or Docker)."""
    cmd = colabfold_cmd or os.environ.get("COLABFOLD_CMD", "colabfold_batch")
    if os.environ.get("COLABFOLD_DOCKER", "").lower() in ("1", "true", "yes"):
        try:
            img = os.environ.get("COLABFOLD_IMAGE", "ghcr.io/sokrypton/colabfold:1.5.3-cuda12.2.2")
            proc = subprocess.run(["docker", "images", "-q", img], capture_output=True, timeout=5)
            return proc.returncode == 0 and len((proc.stdout or b"").strip()) > 0
        except Exception:
            return False
    try:
        proc = subprocess.run([cmd, "--help"], capture_output=True, timeout=5)
        return proc.returncode == 0 or b"colabfold" in (proc.stdout or proc.stderr or b"").lower()
    except Exception:
        return False


def run_alphafold_single(
    sequence: str,
    job_name: str,
    out_dir: str,
    uniprot_id: Optional[str] = None,
) -> dict:
    """
    Run AlphaFold on single sequence. Pluggable: local, remote, or cached.
    Returns: status, pdb_path, plddt_path, pae_path, mean_plddt, per_res_plddt.
    When cached files exist, loads pLDDT from JSON.
    """
    cache_key = _seq_hash(sequence)
    cache_dir = os.path.join(out_dir, "cache", cache_key)
    os.makedirs(cache_dir, exist_ok=True)

    pdb_path = os.path.join(cache_dir, "ranked_0.pdb")
    plddt_path = os.path.join(cache_dir, "plddt.json")
    pae_path = os.path.join(cache_dir, "pae.json")

    if os.path.isfile(pdb_path) and os.path.isfile(plddt_path):
        try:
            with open(plddt_path, "r") as f:
                p = json.load(f)
            per_res = p.get("plddt", p.get("confidence", [])) or []
            if isinstance(per_res, (list, tuple)) and per_res:
                mean_plddt = float(sum(per_res) / len(per_res))
            else:
                mean_plddt = 0.0
                per_res = []
            return {
                "status": "cached",
                "pdb_path": pdb_path,
                "plddt_path": plddt_path,
                "pae_path": pae_path if os.path.isfile(pae_path) else None,
                "mean_plddt": mean_plddt,
                "per_res_plddt": list(per_res),
                "how_to_generate_cache": HOW_TO_GENERATE.format(cache_dir=cache_dir),
            }
        except Exception:
            pass

    # Try AlphaFold DB when cache missing and uniprot_id provided (WT only)
    if uniprot_id and fetch_from_alphafold_db(uniprot_id, sequence, out_dir):
        return run_alphafold_single(sequence, job_name, out_dir, uniprot_id=None)

    return {
        "status": "missing",
        "pdb_path": None,
        "plddt_path": None,
        "pae_path": None,
        "mean_plddt": 0.0,
        "per_res_plddt": [],
        "how_to_generate_cache": HOW_TO_GENERATE.format(cache_dir=cache_dir),
    }


def run_wt_mut_alphafold(
    wt_seq: str,
    mut_seq: str,
    variant_id: str,
    out_dir: str,
    uniprot_id: Optional[str] = None,
) -> dict:
    """Run AF for WT and mutant. Returns wt and mut result dicts. Fetches WT from AlphaFold DB when cache missing."""
    wt = run_alphafold_single(wt_seq, f"wt_{variant_id}", out_dir, uniprot_id=uniprot_id)
    mut = run_alphafold_single(mut_seq, f"mut_{variant_id}", out_dir)
    return {"wt": wt, "mut": mut}
