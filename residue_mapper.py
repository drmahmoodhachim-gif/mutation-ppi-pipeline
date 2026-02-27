"""Map UniProt positions to PDB residue numbers via sequence alignment."""

from typing import Optional

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def extract_chain_residues(pdb_text: str, chain: str) -> list[tuple[int, str, str]]:
    """
    Return list of residues in PDB order: [(resseq, icode, aa1), ...].
    Unique by (resseq, icode); insertion codes preserved.
    """
    if not pdb_text:
        return []
    seen = set()
    residues = []
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        if len(line) < 27:
            continue
        ch = line[21:22].strip()
        if ch != chain:
            continue
        resname = line[17:20].strip().upper()
        resseq_str = line[22:26].strip()
        icode = line[26:27].strip() if len(line) > 26 else ""
        try:
            resseq = int(resseq_str)
        except ValueError:
            continue
        key = (resseq, icode)
        if key in seen:
            continue
        seen.add(key)
        aa1 = AA3_TO_1.get(resname, "X")
        residues.append((resseq, icode, aa1))
    residues.sort(key=lambda x: (x[0], x[1]))
    return residues


def map_uniprot_to_pdb(uniprot_seq: str, pdb_residues: list) -> dict[int, int]:
    """
    Map UniProt position (1-based) -> PDB resseq (int).
    Uses Biopython PairwiseAligner. Only aligned, non-gap positions.
    Returns {} if alignment fails.
    """
    if not uniprot_seq or not pdb_residues:
        return {}
    pdb_seq = "".join(aa for (_, _, aa) in pdb_residues)
    try:
        from Bio.Align import PairwiseAligner
        from Bio.Seq import Seq
        aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1, open_gap_score=-5, extend_gap_score=-0.5)
        alns = aligner.align(Seq(uniprot_seq), Seq(pdb_seq))
        if not alns:
            return {}
        a = alns[0]
        mapping = {}
        ui, pi = 0, 0
        for ua, pa in zip(str(a[0]), str(a[1])):
            if ua != "-":
                ui += 1
            if pa != "-":
                pi += 1
            if ua != "-" and pa != "-" and pi <= len(pdb_residues):
                mapping[ui] = pdb_residues[pi - 1][0]
        return mapping
    except ImportError:
        try:
            from Bio import pairwise2
            alns = pairwise2.align.globalms(uniprot_seq, pdb_seq, 2, -1, -5, -0.5, one_alignment_only=True)
            if not alns:
                return {}
            aln = alns[0]
            u_aln, p_aln = aln.seqA, aln.seqB
            mapping = {}
            u_i, p_i = 0, 0
            for a, b in zip(u_aln, p_aln):
                if a != "-":
                    u_i += 1
                if b != "-":
                    p_i += 1
                if a != "-" and b != "-" and p_i <= len(pdb_residues):
                    mapping[u_i] = pdb_residues[p_i - 1][0]
            return mapping
        except Exception:
            return {}
    except Exception:
        return {}


def map_uniprot_pos_to_pdb_resseq(pdb_str: str, chain: str, uniprot_seq: str) -> dict[int, int]:
    """
    Convenience: extract chain residues and map UniProt -> PDB resseq.
    Returns dict: uniprot_pos (1-based) -> pdb_resseq (int).
    """
    residues = extract_chain_residues(pdb_str, chain)
    return map_uniprot_to_pdb(uniprot_seq, residues)
