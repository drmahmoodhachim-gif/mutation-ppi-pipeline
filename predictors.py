"""Prediction modules for mutation effect on protein structure and PPI."""

import re
import requests
from typing import Optional
from config import (
    GENE_UNIPROT,
    CARDIAC_MYOCYTE_INTERACTORS,
    DEFAULT_INTERACTORS,
    PROTEIN_PDB,
    ALPHAMISSENSE_API,
)

# Single-letter amino acid code
AA_3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}


def parse_mutation(mutation_input: str) -> dict:
    """
    Parse mutation input in formats:
    - c.1577G>A
    - p.R526H
    - R526H
    Returns: {gene, wt_aa, position, mut_aa, cds_pos}
    """
    result = {"gene": None, "wt_aa": None, "position": None, "mut_aa": None, "cds_pos": None}

    # p.R526H or R526H
    aa_match = re.search(r"[p.]?\s*([ARNDCEQGHILKMFPSTWYV])\s*(\d+)\s*([ARNDCEQGHILKMFPSTWYV])", mutation_input, re.I)
    if aa_match:
        result["wt_aa"] = aa_match.group(1).upper()
        result["position"] = int(aa_match.group(2))
        result["mut_aa"] = aa_match.group(3).upper()

    # c.1577G>A
    cds_match = re.search(r"[c.]?\s*(\d+)\s*([ACGT])\s*[>]\s*([ACGT])", mutation_input, re.I)
    if cds_match:
        result["cds_pos"] = int(cds_match.group(1))

    return result


def _parse_alphamissense_categories(raw: str) -> set:
    """Parse '6:C,G,H,L,P,S' -> {'C','G','H','L','P','S'}."""
    if not raw or not isinstance(raw, str):
        return set()
    part = raw.split(":", 1)[-1].strip() if ":" in raw else raw
    return {a.strip().upper() for a in part.split(",") if a.strip()}


def get_alphamissense_prediction(gene: str, position: int, wt_aa: str, mut_aa: str) -> dict:
    """Fetch AlphaMissense pathogenicity prediction via REST API."""
    uniprot_id = GENE_UNIPROT.get(gene.upper())
    if not uniprot_id:
        return {"error": f"Gene {gene} not in database. Add to config.GENE_UNIPROT."}

    url = f"{ALPHAMISSENSE_API}?uid={uniprot_id}&resi={position}"
    try:
        r = requests.get(url, timeout=20)
        r.raise_for_status()
        data = r.json()
        if not isinstance(data, dict):
            return {"pathogenicity": None, "raw": data, "error": "Unexpected API response"}

        mut = mut_aa.upper()
        # API returns benign/pathogenic/ambiguous per residue (e.g. "6:C,G,H,L,P,S")
        benign = _parse_alphamissense_categories(data.get("benign", ""))
        pathogenic = _parse_alphamissense_categories(data.get("pathogenic", ""))
        ambiguous = _parse_alphamissense_categories(data.get("ambiguous", ""))

        if mut in pathogenic:
            score = 0.9
        elif mut in ambiguous:
            score = 0.5
        elif mut in benign:
            score = 0.15
        else:
            mean = data.get("mean_all") or data.get("mean")
            score = float(mean) if mean is not None else 0.5

        return {"pathogenicity": score, "raw": data}
    except requests.RequestException as e:
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


def get_tissue_interactors(gene: str, tissue: str) -> list:
    """Get protein interaction partners for tissue of interest."""
    gene_upper = gene.upper()
    t = tissue.lower()
    if any(kw in t for kw in ("cardiac", "heart", "cardiomyocyte", "myocyte")):
        return CARDIAC_MYOCYTE_INTERACTORS.get(gene_upper, [])
    return []


def get_recommended_pdb(gene: str) -> list:
    """Get recommended PDB IDs for structure visualization."""
    return PROTEIN_PDB.get(gene.upper(), [])


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
