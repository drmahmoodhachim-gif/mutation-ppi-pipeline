"""Structure comparison metrics for WT vs mutant."""

import json
import os
from typing import Optional


def _extract_plddt(plddt_json_path: str, pos: int, window: int) -> Optional[list]:
    """Extract pLDDT scores from JSON file path. Returns list for window around pos, or None."""
    if not plddt_json_path or not os.path.isfile(plddt_json_path):
        return None
    try:
        with open(plddt_json_path, "r") as f:
            data = json.load(f)
        scores = data.get("plddt", data.get("confidence", []))
        if not scores or not isinstance(scores, (list, tuple)):
            return None
        start = max(0, pos - 1 - window // 2)
        end = min(len(scores), pos + window // 2)
        return list(scores[start:end])
    except Exception:
        return None


def compute_local_structure_deltas(
    wt_pdb: str,
    mut_pdb: str,
    pos: int,
    window: int = 10,
) -> dict:
    """
    Legacy: accepts PDB paths (unused). Use compute_local_structure_deltas_from_af instead.
    """
    return {
        "delta_mean_plddt_window": None,
        "wt_mean_plddt_window": None,
        "mut_mean_plddt_window": None,
        "local_rmsd_window": None,
        "region_type": "LowConfidence/IDR-like",
    }


def compute_local_structure_deltas_from_af(
    wt_af: dict,
    mut_af: dict,
    pos: int,
    window: int = 10,
) -> dict:
    """
    Compute local structure deltas from AlphaFold result dicts.
    wt_af / mut_af must have plddt_path (file path to pLDDT JSON).
    Returns: delta_mean_plddt_window, wt_mean_plddt_window, mut_mean_plddt_window,
    local_rmsd_window (or None), region_type.
    """
    wt_plddt = _extract_plddt(wt_af.get("plddt_path"), pos, window)
    mut_plddt = _extract_plddt(mut_af.get("plddt_path"), pos, window)
    if wt_plddt is None or mut_plddt is None:
        return {
            "delta_mean_plddt_window": None,
            "wt_mean_plddt_window": None,
            "mut_mean_plddt_window": None,
            "local_rmsd_window": None,
            "region_type": "LowConfidence/IDR-like",
        }
    wt_mean = sum(wt_plddt) / len(wt_plddt) if wt_plddt else 0
    mut_mean = sum(mut_plddt) / len(mut_plddt) if mut_plddt else 0
    delta = mut_mean - wt_mean
    region = "Structured" if (wt_mean >= 60 and min(wt_plddt) >= 50) else "LowConfidence/IDR-like"
    return {
        "delta_mean_plddt_window": round(delta, 2),
        "wt_mean_plddt_window": round(wt_mean, 2),
        "mut_mean_plddt_window": round(mut_mean, 2),
        "local_rmsd_window": None,
        "region_type": region,
    }
