"""Structure comparison metrics for WT vs mutant."""


def compute_local_structure_deltas(
    wt_pdb: str,
    mut_pdb: str,
    pos: int,
    window: int = 10,
) -> dict:
    """
    Compute local structure deltas. Returns:
    delta_mean_plddt_window, wt_mean_plddt_window, mut_mean_plddt_window,
    local_rmsd_window (or None), region_type.
    """
    wt_plddt = _extract_plddt(wt_pdb, pos, window)
    mut_plddt = _extract_plddt(mut_pdb, pos, window)
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


def _extract_plddt(pdb_or_json: str, pos: int, window: int):
    """Extract pLDDT scores from AlphaFold JSON/PAE or return None."""
    import json
    try:
        if pdb_or_json.strip().startswith("{"):
            data = json.loads(pdb_or_json)
            scores = data.get("plddt", data.get("confidence", []))
            if scores and isinstance(scores, (list, tuple)):
                start = max(0, pos - 1 - window // 2)
                end = min(len(scores), pos + window // 2)
                return list(scores[start:end])
    except Exception:
        pass
    return None
