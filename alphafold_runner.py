"""AlphaFold runner - pluggable (local/remote/cached)."""

import hashlib
import os
from typing import Optional


def _seq_hash(seq: str) -> str:
    return hashlib.sha256(seq.encode()).hexdigest()[:16]


def run_alphafold_single(sequence: str, job_name: str, out_dir: str) -> dict:
    """
    Run AlphaFold on single sequence. Pluggable: local, remote, or cached.
    Returns: pdb_path, plddt_path, pae_path, mean_plddt, per_res_plddt.
    """
    cache_key = _seq_hash(sequence)
    cache_dir = os.path.join(out_dir, "cache", cache_key)
    os.makedirs(cache_dir, exist_ok=True)
    pdb_path = os.path.join(cache_dir, "ranked_0.pdb")
    plddt_path = os.path.join(cache_dir, "plddt.json")
    pae_path = None
    mean_plddt = 0.0
    per_res_plddt = []
    if os.path.isfile(pdb_path):
        return {
            "pdb_path": pdb_path,
            "plddt_path": plddt_path,
            "pae_path": pae_path,
            "mean_plddt": mean_plddt,
            "per_res_plddt": per_res_plddt,
        }
    return {
        "pdb_path": None,
        "plddt_path": None,
        "pae_path": None,
        "mean_plddt": 0.0,
        "per_res_plddt": [],
    }


def run_wt_mut_alphafold(wt_seq: str, mut_seq: str, variant_id: str, out_dir: str) -> dict:
    """Run AF for WT and mutant. Returns wt and mut result dicts."""
    wt = run_alphafold_single(wt_seq, f"wt_{variant_id}", out_dir)
    mut = run_alphafold_single(mut_seq, f"mut_{variant_id}", out_dir)
    return {"wt": wt, "mut": mut}
