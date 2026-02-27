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

# Single-letter amino acid code
AA_3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}


def parse_mutation(mutation_input: str, gene: str = None) -> dict:
    """Parse mutation. Returns gene, wt, pos, mut, variant_id, wt_aa, position, mut_aa, cds_pos."""
    result = {"gene": gene, "wt_aa": None, "position": None, "mut_aa": None, "cds_pos": None, "wt": None, "pos": None, "mut": None, "variant_id": None}

    aa3_match = re.search(r"[p.]?\s*(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)\s*(\d+)\s*(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)", mutation_input, re.I)
    if aa3_match:
        result["wt_aa"] = AA_3_TO_1.get(aa3_match.group(1).capitalize(), aa3_match.group(1)[0])
        result["position"] = int(aa3_match.group(2))
        result["mut_aa"] = AA_3_TO_1.get(aa3_match.group(3).capitalize(), aa3_match.group(3)[0])
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


def get_alphamissense_prediction(gene: str, position: int, wt_aa: str, mut_aa: str) -> dict:
    """Fetch AlphaMissense pathogenicity prediction via REST API."""
    uniprot_id = GENE_UNIPROT.get(gene.upper())
    if not uniprot_id:
        return {"error": f"Gene {gene} not in database. Add to config.GENE_UNIPROT."}

    url = f"{ALPHAMISSENSE_API}?uid={uniprot_id}&resi={position}"
    try:
        r = requests.get(url, timeout=15)
        r.raise_for_status()
        data = r.json()
        # API returns all substitutions; find our mutation
        if isinstance(data, dict) and "score" in data:
            return {"pathogenicity": data.get("score", 0), "raw": data}
        if isinstance(data, list):
            for item in data:
                if isinstance(item, dict) and str(item.get("aa", "")) == mut_aa:
                    return {"pathogenicity": item.get("score", 0), "raw": item}
        return {"pathogenicity": None, "raw": data, "note": "Mutation not found in response"}
    except Exception as e:
        return {"error": str(e), "url": url}


def resolve_uniprot(gene: str) -> Optional[str]:
    """Resolve gene symbol to UniProt ID via config or REST."""
    uid = GENE_UNIPROT.get(gene.upper())
    if uid:
        return uid
    try:
        r = requests.get(
            "https://rest.uniprot.org/uniprotkb/search",
            params={
                "query": f"gene_exact:{gene} AND (organism_id:9606)",
                "format": "json",
                "size": 1,
                "fields": "accession",
            },
            timeout=10,
        )
        if r.ok:
            data = r.json()
            results = data.get("results", [])
            if results:
                return results[0].get("primaryAccession")
    except Exception:
        pass
    return None



def validate_variant_on_canonical_sequence(canonical_fasta: str, wt: str, pos: int) -> dict:
    if not canonical_fasta or not wt or pos is None: return {'wt_from_seq': None, 'match': False, 'message': 'Missing'}
    if pos < 1 or pos > len(canonical_fasta): return {'wt_from_seq': None, 'match': False, 'message': 'Out of range'}
    wt_from_seq = canonical_fasta[pos - 1].upper()
    match = wt_from_seq == wt.upper()
    return {'wt_from_seq': wt_from_seq, 'match': match, 'message': 'OK' if match else 'Mismatch'}

def get_tissue_interactors(gene: str, tissue: str) -> list:
    """Get protein interaction partners for tissue of interest."""
    gene_upper = gene.upper()
    t = tissue.lower()
    if any(kw in t for kw in ("cardiac", "heart", "cardiomyocyte", "myocyte")):
        return CARDIAC_MYOCYTE_INTERACTORS.get(gene_upper, [])
    return []


def get_recommended_pdb(gene: str, pos: int = None):
    pdb_ids = PROTEIN_PDB.get(gene.upper(), [])
    rec = {"pdb_ids": pdb_ids, "has_residue_coordinates": True}
    if pos and pdb_ids:
        try:
            from visualization import fetch_pdb, residue_in_pdb
            pd = fetch_pdb(pdb_ids[0])
            rec["has_residue_coordinates"] = bool(pd and residue_in_pdb(pd,pos,"A"))
        except: pass
    return rec


    """Get recommended PDB IDs for structure visualization."""
    return PROTEIN_PDB.get(gene.upper(), [])



def _heuristic_ddg(wt_aa, mut_aa):
    charge = {"R":1,"K":1,"D":-1,"E":-1,"H":0.5}
    if abs(charge.get(wt_aa,0)-charge.get(mut_aa,0))>0.5: return 2.0
    h = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,"G":-0.4,"H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,"P":-1.6,"S":-0.8,"T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}
    return 1.2 if abs(h.get(wt_aa,0)-h.get(mut_aa,0))>2 else 0.3

def get_ppi_ddg_predictions(gene, wt_aa, mut_aa, position, interactors):
    rows = []
    for ip in interactors:
        pdb_id = PPI_PDB_COMPLEXES.get(gene.upper()+"_"+ip.get("partner","")) if isinstance(PPI_PDB_COMPLEXES,dict) else None
        has_c = bool(pdb_id and len(str(pdb_id))==4)
        ddg, conc, notes = None, "NA", "No interface model; residue not confirmed at interface."
        if has_c:
            try:
                from interface import is_residue_at_interface
                from visualization import fetch_pdb
                pd = fetch_pdb(str(pdb_id))
                if pd:
                    ir = is_residue_at_interface(pd,"A","B",position,8.0)
                    if ir.get("is_interface"):
                        ddg = round(_heuristic_ddg(wt_aa,mut_aa),2)
                        conc = "LikelyDisruptive" if ddg>1.0 else ("LikelyStabilizing" if ddg<-0.5 else "Neutral/Unclear")
                        notes = "OK"
                    else: notes = "Not interface"
            except: notes = "Complex unavailable"
        rows.append({"Interacting protein":ip.get("partner","?"),"Role":ip.get("role",""),"Tissue evidence":"Present","Complex model available?":"Y" if has_c else "N","Interface residue?":"-","PPI dG (kcal/mol)":ddg if ddg else "NA","Conclusion":conc,"Notes":notes})
    return rows

def estimate_structural_impact(wt_aa: str, mut_aa: str, position: int, in_voltage_sensor: bool = False) -> dict:
    """
    Heuristic for structural impact based on physicochemical properties.
    """
    charge_change = {"R": 1, "K": 1, "D": -1, "E": -1, "H": 0.5}
    hydrophobicity = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
                     "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
                     "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
                     "Y": -1.3, "V": 4.2}
    c_wt = charge_change.get(wt_aa, 0)
    c_mut = charge_change.get(mut_aa, 0)
    h_wt = hydrophobicity.get(wt_aa, 0)
    h_mut = hydrophobicity.get(mut_aa, 0)

    impact = "Low"
    reasons = []
    if abs(c_wt - c_mut) > 0.5:
        impact = "High" if in_voltage_sensor else "Medium"
        reasons.append(f"Charge change: {wt_aa}({c_wt}) → {mut_aa}({c_mut})")
    if abs(h_wt - h_mut) > 2:
        if impact != "High":
            impact = "Medium"
        reasons.append(f"Hydrophobicity change: Δ={abs(h_wt - h_mut):.1f}")
    if not reasons:
        reasons.append("Conservative substitution")

    return {"impact": impact, "reasons": reasons}
