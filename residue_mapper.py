"""Map UniProt positions to PDB residue numbers via sequence alignment."""

from typing import Optional

# 3-letter to 1-letter amino acid code
AA_3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def _extract_chain_sequence(pdb_str: str, chain: str) -> list[tuple[int, str]]:
    """Extract (resseq, 1-letter aa) for chain from PDB. Ordered by residue number."""
    seq = []
    seen = set()
    for line in pdb_str.split("\n"):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if len(line) < 55:
            continue
        if line[21:22].strip() != chain:
            continue
        try:
            resseq = int(line[22:26].strip().lstrip("-") or 0)
            resname = line[17:20].strip().upper()
        except (ValueError, IndexError):
            continue
        if (resseq, resname) in seen:
            continue
        seen.add((resseq, resname))
        aa = AA_3_TO_1.get(resname, "X")
        seq.append((resseq, aa))
    seq.sort(key=lambda x: x[0])
    return seq


def map_uniprot_pos_to_pdb_resseq(
    pdb_str: str,
    chain: str,
    uniprot_seq: str,
) -> dict[int, int]:
    """
    Map UniProt position (1-based) -> PDB residue number (resseq).
    Uses sequential alignment: assumes PDB chain order matches UniProt order where aligned.
    Returns {} if alignment fails.
    """
    if not pdb_str or not uniprot_seq:
        return {}
    chain_seq = _extract_chain_sequence(pdb_str, chain)
    if not chain_seq:
        return {}
    pdb_resseqs = [r[0] for r in chain_seq]
    pdb_aas = "".join(r[1] for r in chain_seq)
    # Simple alignment: assume 1:1 where sequences match; handle gaps minimally.
    try:
        from Bio.Align import PairwiseAligner
        from Bio.Seq import Seq
        aligner = PairwiseAligner(mode="global", match_score=1, mismatch_score=-1, open_gap_score=-2, extend_gap_score=-1)
        alns = aligner.align(Seq(uniprot_seq), Seq(pdb_aas))
        if not alns:
            return {}
        a = alns[0]
        mapping = {}
        ui, pi = 0, 0
        for ua, pa in zip(str(a[0]), str(a[1])):
            if ua != "-" and pa != "-":
                ui += 1
                if pi < len(pdb_resseqs):
                    mapping[ui] = pdb_resseqs[pi]
                pi += 1
            elif ua != "-":
                ui += 1
            elif pa != "-":
                pi += 1
        return mapping
    except ImportError:
        # Fallback: assume 1:1 if lengths similar and no obvious gaps
        if abs(len(pdb_aas) - len(uniprot_seq)) < 50 and len(chain_seq) > 0:
            out = {}
            for i in range(min(len(uniprot_seq), len(pdb_resseqs))):
                out[i + 1] = pdb_resseqs[i]
            return out
        return {}
    except Exception:
        return {}
