"""Prediction modules for mutation effect on protein structure and PPI."""

import re
import requests
from typing import Optional
from config import (
    GENE_UNIPROT,
    CARDIAC_MYOCYTE_INTERACTORS,
    DEFAULT_INTERACTORS,
    PROTEIN_PDB,
    PPI_PDB_COMPLEXES,
    ALPHAMISSENSE_API,
)

AA_3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}


def parse_mutation(mutation_input: str, gene: str = None) -> dict:
    """
    Returns: gene, uniprot_id, wt, pos, mut, variant_id, wt_aa, position, mut_aa, cds_pos, warnings.
    """
    warnings = []
    result = {"gene": gene, "uniprot_id": None, "wt": None, "pos": None, "mut": None,
              "variant_id": None, "wt_aa": None, "position": None, "mut_aa": None, "cds_pos": None, "warnings": warnings}
    if gene:
        result["uniprot_id"] = GENE_UNIPROT.get(gene.upper()) or _resolve_uniprot_id(gene)

    aa3 = re.search(r"[p.]?\s*(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)\s*(\d+)\s*(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)", mutation_input, re.I)
    if aa3:
        result["wt_aa"] = AA_3_TO_1.get(aa3.group(1).capitalize(), aa3.group(1)[0])
        result["position"] = int(aa3.group(2))
        result["mut_aa"] = AA_3_TO_1.get(aa3.group(3).capitalize(), aa3.group(3)[0])
    if not result.get("wt_aa"):
        aa_match = re.search(r"[p.]?\s*([ARNDCEQGHILKMFPSTWYV])\s*(\d+)\s*([ARNDCEQGHILKMFPSTWYV])", mutation_input, re.I)
        if aa_match:
            result["wt_aa"] = aa_match.group(1).upper()
            result["position"] = int(aa_match.group(2))
            result["mut_aa"] = aa_match.group(3).upper()

    result["wt"], result["pos"], result["mut"] = result["wt_aa"], result["position"], result["mut_aa"]
    if result.get("gene") and result.get("wt") and result.get("pos") and result.get("mut"):
        result["variant_id"] = f"{result['gene']}:{result['wt']}{result['pos']}{result['mut']}"

    cds_match = re.search(r"[c.]?\s*(\d+)\s*([ACGT])\s*[>]\s*([ACGT])", mutation_input, re.I)
    if cds_match:
        result["cds_pos"] = int(cds_match.group(1))
    return result


def _resolve_uniprot_id(gene: str) -> Optional[str]:
    uid = GENE_UNIPROT.get(gene.upper())
    if uid:
        return uid
    try:
        r = requests.get("https://rest.uniprot.org/uniprotkb/search",
            params={"query": f"gene_exact:{gene} AND (organism_id:9606)", "format": "json", "size": 1, "fields": "accession"}, timeout=10)
        if r.ok and r.json().get("results"):
            return r.json()["results"][0].get("primaryAccession")
    except Exception:
        pass
    return None


def resolve_uniprot(gene: str) -> Optional[dict]:
    """
    Returns: {uniprot_id, canonical_fasta, length, protein_name, organism}
    """
    uid = GENE_UNIPROT.get(gene.upper()) or _resolve_uniprot_id(gene)
    if not uid:
        return None
    try:
        r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uid}",
            params={"format": "json", "fields": "accession,protein_name,organism_name,sequence"}, timeout=15)
        if not r.ok:
            return {"uniprot_id": uid, "canonical_fasta": "", "length": 0, "protein_name": "", "organism": ""}
        d = r.json()
        seq = d.get("sequence", {}).get("value", "")
        name = d.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "") or ""
        org = d.get("organism", {}).get("scientificName", "")
        return {"uniprot_id": uid, "canonical_fasta": seq, "length": len(seq), "protein_name": name, "organism": org}
    except Exception:
        return {"uniprot_id": uid, "canonical_fasta": "", "length": 0, "protein_name": "", "organism": ""}


def validate_variant_on_canonical_sequence(canonical_fasta: str, wt: str, pos: int) -> dict:
    """Returns: {wt_from_seq, match, message}"""
    if not canonical_fasta or not wt or pos is None:
        return {"wt_from_seq": None, "match": False, "message": "Missing sequence, WT, or position."}
    if pos < 1 or pos > len(canonical_fasta):
        return {"wt_from_seq": None, "match": False, "message": f"Position {pos} out of range (1-{len(canonical_fasta)})."}
    wt_from_seq = canonical_fasta[pos - 1].upper()
    match = wt_from_seq == wt.upper()
    msg = "OK" if match else f"WT mismatch: canonical has {wt_from_seq} at {pos}, variant claims {wt}."
    return {"wt_from_seq": wt_from_seq, "match": match, "message": msg}


def get_alphamissense_prediction(gene: str, position: int, wt_aa: str, mut_aa: str) -> dict:
    uniprot_id = GENE_UNIPROT.get(gene.upper()) or _resolve_uniprot_id(gene)
    if not uniprot_id:
        return {"error": f"Gene {gene} not in database."}
    url = f"{ALPHAMISSENSE_API}?uid={uniprot_id}&resi={position}"
    try:
        r = requests.get(url, timeout=15)
        r.raise_for_status()
        data = r.json()
        if isinstance(data, dict) and "score" in data:
            return {"pathogenicity": data.get("score", 0), "raw": data}
        if isinstance(data, list):
            for item in data:
                if isinstance(item, dict) and str(item.get("aa", "")) == mut_aa:
                    return {"pathogenicity": item.get("score", 0), "raw": item}
        return {"pathogenicity": None, "raw": data, "note": "Mutation not found"}
    except Exception as e:
        return {"error": str(e), "url": url}


def get_tissue_interactors(gene: str, tissue: str) -> list:
    gene_upper = gene.upper()
    t = tissue.lower()
    if any(kw in t for kw in ("cardiac", "heart", "cardiomyocyte", "myocyte")):
        return CARDIAC_MYOCYTE_INTERACTORS.get(gene_upper, [])
    return []


def get_recommended_pdb(gene: str, pos: int = None) -> dict:
    """Returns: {pdb_id, chain, has_residue_coordinates, mapping_method, pdb_ids, notes}"""
    pdb_ids = PROTEIN_PDB.get(gene.upper(), [])
    pdb_id = pdb_ids[0] if pdb_ids else None
    has_res = True
    notes = ""
    if pos and pdb_ids:
        try:
            from visualization import fetch_pdb, residue_in_pdb
            pd = fetch_pdb(pdb_ids[0])
            has_res = bool(pd and residue_in_pdb(pd, pos, "A"))
            if not has_res:
                notes = "Residue not in PDB; recommend AlphaFold."
        except Exception:
            notes = "Could not verify residue."
    return {"pdb_id": pdb_id, "chain": "A", "has_residue_coordinates": has_res,
            "mapping_method": "config", "pdb_ids": pdb_ids, "notes": notes or "OK"}


def _heuristic_ddg(wt_aa: str, mut_aa: str) -> float:
    charge = {"R": 1, "K": 1, "D": -1, "E": -1, "H": 0.5}
    if abs(charge.get(wt_aa, 0) - charge.get(mut_aa, 0)) > 0.5:
        return 2.0
    h = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
         "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
         "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
    return 1.2 if abs(h.get(wt_aa, 0) - h.get(mut_aa, 0)) > 2 else 0.3


def get_ppi_ddg_predictions(gene: str, wt_aa: str, mut_aa: str, position: int, interactors: list) -> list:
    """PPI Delta-Delta-G per partner. NA when no complex or not at interface."""
    rows = []
    for ip in interactors:
        partner, role = ip.get("partner", "?"), ip.get("role", "")
        key = f"{gene.upper()}_{partner}"
        pdb_id = PPI_PDB_COMPLEXES.get(key) if isinstance(PPI_PDB_COMPLEXES, dict) else None
        if isinstance(pdb_id, dict):
            pdb_id = pdb_id.get("pdb_id")

        has_c = bool(pdb_id and isinstance(pdb_id, str) and len(pdb_id) == 4)
        complex_model = None
        interface_info = None
        ppi_ddg = None
        method = "Heuristic"
        conclusion = "NA"
        notes = "No interface model; residue not confirmed at interface."
        got_pdb = False

        if has_c:
            try:
                from interface import is_residue_at_interface
                from visualization import fetch_pdb
                pd = fetch_pdb(str(pdb_id))
                got_pdb = bool(pd)
                if pd:
                    ir = is_residue_at_interface(pd, "A", "B", position, 8.0)
                    interface_info = ir
                    if ir.get("is_interface"):
                        ppi_ddg = round(_heuristic_ddg(wt_aa, mut_aa), 2)
                        method = "Heuristic"
                        conclusion = "LikelyDisruptive" if ppi_ddg > 1.0 else ("LikelyStabilizing" if ppi_ddg < -0.5 else "Neutral/Unclear")
                        notes = "OK"
                    else:
                        notes = f"Not interface (min_dist={ir.get('min_distance', 0):.1f} A)"
                else:
                    notes = "Complex model unavailable."
            except Exception:
                notes = "Complex unavailable."

        rows.append({
            "Interacting protein": partner,
            "Role": role,
            "Tissue evidence": "Present",
            "Complex model available?": "Y" if (has_c and got_pdb) else "N",
            "Interface residue?": "Y" if (interface_info and interface_info.get("is_interface")) else ("N" if interface_info else "—"),
            "PPI Delta-Delta-G (kcal/mol)": ppi_ddg if ppi_ddg is not None else "NA",
            "Conclusion": conclusion,
            "Notes": notes,
        })
    return rows


def estimate_structural_impact(wt_aa: str, mut_aa: str, position: int, in_voltage_sensor: bool = False) -> dict:
    charge_change = {"R": 1, "K": 1, "D": -1, "E": -1, "H": 0.5}
    hydrophobicity = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
                     "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
                     "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
                     "Y": -1.3, "V": 4.2}
    c_wt, c_mut = charge_change.get(wt_aa, 0), charge_change.get(mut_aa, 0)
    h_wt, h_mut = hydrophobicity.get(wt_aa, 0), hydrophobicity.get(mut_aa, 0)
    impact = "Low"
    reasons = []
    if abs(c_wt - c_mut) > 0.5:
        impact = "High" if in_voltage_sensor else "Medium"
        reasons.append(f"Charge change: {wt_aa}({c_wt}) -> {mut_aa}({c_mut})")
    if abs(h_wt - h_mut) > 2:
        if impact != "High":
            impact = "Medium"
        reasons.append(f"Hydrophobicity change: Delta={abs(h_wt - h_mut):.1f}")
    if not reasons:
        reasons.append("Conservative substitution")
    return {"impact": impact, "reasons": reasons}
