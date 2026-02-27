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


def parse_mutation(mutation_input: str) -> dict:
    """
    Parse mutation input in formats:
    - c.1577G>A
    - p.R526H
    - R526H
    Returns: {gene, wt_aa, position, mut_aa, cds_pos}
    """
    result = {"gene": None, "wt_aa": None, "position": None, "mut_aa": None, "cds_pos": None}

    # p.Ser1054Ala — 3-letter codes first
    aa3_match = re.search(
        r"[p.]?\s*(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)\s*(\d+)\s*(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)",
        mutation_input, re.I,
    )
    if aa3_match:
        result["wt_aa"] = AA_3_TO_1.get(aa3_match.group(1).capitalize(), aa3_match.group(1)[0])
        result["position"] = int(aa3_match.group(2))
        result["mut_aa"] = AA_3_TO_1.get(aa3_match.group(3).capitalize(), aa3_match.group(3)[0])
    # p.R526H or R526H — single-letter
    if not result.get("wt_aa"):
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


def _heuristic_ddg(wt_aa: str, mut_aa: str) -> float:
    """Heuristic ΔΔG (kcal/mol). Positive = destabilizing."""
    charge = {"R": 1, "K": 1, "D": -1, "E": -1, "H": 0.5}
    c_wt, c_mut = charge.get(wt_aa, 0), charge.get(mut_aa, 0)
    if abs(c_wt - c_mut) > 0.5:
        return 2.0
    hydrophobicity = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
                      "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
                      "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
                      "Y": -1.3, "V": 4.2}
    h_wt, h_mut = hydrophobicity.get(wt_aa, 0), hydrophobicity.get(mut_aa, 0)
    return 1.2 if abs(h_wt - h_mut) > 2 else 0.3

def get_ppi_ddg_predictions(gene: str, wt_aa: str, mut_aa: str, position: int, interactors: list) -> list:
    """Build PPI ΔΔG prediction rows for table."""
    ddg = _heuristic_ddg(wt_aa, mut_aa)
    rows = []
    for ip in interactors:
        conclusion = "Destabilizing" if ddg > 1.0 else ("Stabilizing" if ddg < -0.5 else "Neutral / unclear")
        pathway = "May reduce binding affinity; downstream pathway may be impaired." if ddg > 1.0 else (
            "May enhance binding; verify functional consequences." if ddg < -0.5 else
            "Small effect; experimental validation recommended.")
        rows.append({
            "Interacting protein": ip.get("partner", "?"),
            "Role": ip.get("role", ""),
            "Wild-type ref": f"{wt_aa}{position}",
            "Mutant ΔΔG (kcal/mol)": round(ddg, 2),
            "Conclusion": conclusion,
            "Pathway effect downstream": pathway,
        })
    return rows

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
