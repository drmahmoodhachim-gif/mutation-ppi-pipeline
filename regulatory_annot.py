"""Motif and PTM proximity analysis for linker/regulatory variants."""

import re


# Common kinase/PTM motifs (simplified patterns)
MOTIF_PATTERNS = [
    ("PKA_site", r"[KR]{2}[^P].[ST]"),
    ("PKC_site", r"[ST].[RK]"),
    ("CK2_site", r"[ST].{0,2}[DE]"),
    ("CaMKII_site", r"[RK]{2}.{0,2}[ST]"),
    ("Proline_rich", r"P{2,}"),
    ("Acidic_patch", r"[DE]{3,}"),
    ("Basic_patch", r"[KR]{2,}"),
]


def scan_motifs(sequence: str, pos: int, window: int = 12) -> dict:
    """
    Scan for motifs in window around position.
    Returns: {
      'window_seq': str,
      'motifs': [{'name': str, 'match': str, 'offset': int}, ...],
    }
    """
    if not sequence or pos is None:
        return {"window_seq": "", "motifs": []}
    one_indexed = max(1, pos)
    start = max(0, one_indexed - 1 - window // 2)
    end = min(len(sequence), one_indexed + window // 2)
    window_seq = sequence[start:end]

    motifs = []
    for name, pattern in MOTIF_PATTERNS:
        for m in re.finditer(pattern, window_seq, re.IGNORECASE):
            offset = start + m.start() + 1
            motifs.append({"name": name, "match": m.group(), "offset": offset})

    return {"window_seq": window_seq, "motifs": motifs}


def ptm_proximity(sequence: str, pos: int) -> dict:
    """
    Check for phosphorylatable S/T nearby.
    Returns: {
      'phospho_site_nearby': bool,
      'nearby_sites': [{'type': str, 'pos': int}, ...],
    }
    """
    if not sequence or pos is None:
        return {"phospho_site_nearby": False, "nearby_sites": []}

    one_indexed = max(1, pos)
    rad = 7
    start = max(0, one_indexed - 1 - rad)
    end = min(len(sequence), one_indexed + rad)
    nearby = []

    for i in range(start, end):
        aa = sequence[i] if i < len(sequence) else ""
        if aa in "ST":
            nearby.append({"type": aa, "pos": i + 1})

    return {
        "phospho_site_nearby": len(nearby) > 0,
        "nearby_sites": nearby,
    }
