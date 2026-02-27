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


def parse_mutation(mutation_input: str, gene: str = None, position: int = None, wt_aa: str = None, mut_aa: str = None) -> dict:
    """
    Returns: gene, uniprot_id, wt, pos, mut, variant_id, wt_aa, position, mut_aa, cds_pos, warnings.
    Supports: p.R526H, R526H, Ser1054Ala, c.1577G>A; also position+wt+mut (526 R H, 526 R→H).
    If position, wt_aa, mut_aa are provided, they override parsing.
    """
    warnings = []
    result = {"gene": gene, "uniprot_id": None, "wt": None, "pos": None, "mut": None,
              "variant_id": None, "wt_aa": None, "position": None, "mut_aa": None, "cds_pos": None, "warnings": warnings}
    if gene:
        result["uniprot_id"] = GENE_UNIPROT.get(gene.upper()) or _resolve_uniprot_id(gene)

    # Direct override when position, wt_aa, mut_aa provided
    if position is not None and wt_aa and mut_aa:
        result["position"] = int(position)
        result["wt_aa"] = wt_aa.upper() if len(wt_aa) == 1 else AA_3_TO_1.get(wt_aa.capitalize()[:3], wt_aa[0])
        result["mut_aa"] = mut_aa.upper() if len(mut_aa) == 1 else AA_3_TO_1.get(mut_aa.capitalize()[:3], mut_aa[0])
    else:
        # Parse from text
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
        if not result.get("wt_aa"):
            # Simple format: position wt mut (526 R H, 526 R→H, 526 R->H, 526 R to H)
            simple = re.search(r"(?:^|position\s*)?(\d+)\s*([ARNDCEQGHILKMFPSTWYV])\s*(?:->|→|to)\s*([ARNDCEQGHILKMFPSTWYV])", mutation_input, re.I)
            if not simple:
                simple = re.search(r"(?:^|\s)(\d+)\s+([ARNDCEQGHILKMFPSTWYV])\s+([ARNDCEQGHILKMFPSTWYV])\s*$", mutation_input, re.I)
            if simple:
                result["position"] = int(simple.group(1))
                result["wt_aa"] = simple.group(2).upper()
                result["mut_aa"] = simple.group(3).upper()

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
            return {"uniprot_id": uid, "canonical_fasta": "", "length": 0, "protein_name": "", "organism": "", "warnings": ["UniProt returned error"]}
        d = r.json()
        seq = d.get("sequence", {}).get("value", "")
        name = d.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "") or ""
        org = d.get("organism", {}).get("scientificName", "")
        return {"uniprot_id": uid, "canonical_fasta": seq, "length": len(seq), "protein_name": name, "organism": org, "warnings": []}
    except Exception:
        return {"uniprot_id": uid, "canonical_fasta": "", "length": 0, "protein_name": "", "organism": "", "warnings": ["UniProt fetch failed"]}


def validate_variant_on_canonical_sequence(canonical_fasta: str, wt: str, pos: int) -> dict:
    """Returns: {wt_from_seq, match, message}. Alias: validate_variant_on_sequence."""
    return validate_variant_on_sequence(canonical_fasta, wt, pos)


def validate_variant_on_sequence(seq: str, wt_aa: str, pos: int) -> dict:
    """Returns: {wt_from_seq, match, message}. Hard-stop on mismatch recommended."""
    if not seq or not wt_aa or pos is None:
        return {"wt_from_seq": None, "match": False, "message": "Missing sequence, WT, or position."}
    if pos < 1 or pos > len(seq):
        return {"wt_from_seq": None, "match": False, "message": f"Position {pos} out of range (1-{len(seq)})."}
    wt_from_seq = seq[pos - 1].upper()
    match = wt_from_seq == wt_aa.upper()
    return {
        "wt_from_seq": wt_from_seq,
        "match": match,
        "message": "OK" if match else f"WT mismatch: input {wt_aa}{pos} but sequence has {wt_from_seq}{pos}",
    }


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


STRING_API = "https://string-db.org/api"
STRING_SPECIES = 9606  # Human


def _fetch_string_interactors(uniprot_id: str, limit: int = 15) -> list:
    """Fetch interaction partners from STRING API for any UniProt ID. Returns list of {partner, label, uniprot, role}."""
    if not uniprot_id or not requests:
        return []
    try:
        r = requests.post(
            f"{STRING_API}/tsv/interaction_partners",
            data={"identifiers": uniprot_id, "species": STRING_SPECIES, "limit": limit, "caller_identity": "mutation-ppi-pipeline"},
            timeout=15,
        )
        if not r.ok or not r.text.strip():
            return []
        lines = r.text.strip().split("\n")
        if len(lines) < 2:
            return []
        header = [h.strip().lower().replace("-", "_") for h in lines[0].split("\t")]
        # interaction_partners: preferredName_B is partner (or column 1/2)
        idx_name = next((i for i, h in enumerate(header) if "preferredname_b" in h or "preferred_name_b" in h or ("preferred" in h and "b" in h)), None)
        if idx_name is None:
            idx_name = next((i for i, h in enumerate(header) if "name" in h and "b" in h), 1)
        idx_score = next((i for i, h in enumerate(header) if "score" in h), -1)
        result = []
        seen = set()
        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) <= max(idx_name, 0):
                continue
            name = (parts[idx_name] or "?").strip()
            if not name or name in seen:
                continue
            seen.add(name)
            score = float(parts[idx_score]) if 0 <= idx_score < len(parts) and parts[idx_score] else 0
            result.append({
                "partner": name,
                "label": name,
                "uniprot": "",
                "role": f"STRING score {score:.0%}" if score else "Interaction partner",
            })
        return result
    except Exception:
        return []


def get_tissue_interactors(gene: str, tissue: str, uniprot_id: str = None) -> list:
    """Get interactors: curated when available, else STRING API for any gene."""
    gene_upper = gene.upper()
    t = tissue.lower()
    curated = []
    if any(kw in t for kw in ("cardiac", "heart", "cardiomyocyte", "myocyte")):
        curated = CARDIAC_MYOCYTE_INTERACTORS.get(gene_upper, [])
    if curated:
        return curated
    # Fallback: fetch from STRING for any gene
    uid = uniprot_id or GENE_UNIPROT.get(gene_upper) or _resolve_uniprot_id(gene)
    if uid:
        return _fetch_string_interactors(uid)
    return []


def _fetch_pdb_by_uniprot(uniprot_id: str, limit: int = 10) -> list:
    """Fetch PDB IDs from RCSB Search API by UniProt ID. Returns list of PDB IDs."""
    if not uniprot_id or not requests:
        return []
    try:
        import json
        query = {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {"type": "terminal", "service": "text", "parameters": {"operator": "exact_match", "value": uniprot_id, "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession"}},
                    {"type": "terminal", "service": "text", "parameters": {"operator": "exact_match", "value": "UniProt", "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name"}},
                ],
            },
            "return_type": "entry",
            "request_options": {"paginate": {"start": 0, "rows": limit}},
        }
        r = requests.post("https://search.rcsb.org/rcsbsearch/v2/query", json=query, timeout=15)
        if not r.ok:
            return []
        data = r.json()
        ids = data.get("result_set", []) or []
        if not isinstance(ids, list):
            return []
        # RCSB returns [{"identifier": "1ABC", ...}, ...]
        return [x.get("identifier", x) if isinstance(x, dict) else str(x) for x in ids[:limit]]
    except Exception:
        return []


def get_recommended_pdb(
    gene: str,
    pos: int = None,
    uniprot_seq: str = None,
    uniprot_id: str = None,
) -> dict:
    """
    Returns: {pdb_id, chain, has_residue_coordinates, mapping_method, pdb_ids, pdb_resseq, notes}.
    Uses PROTEIN_PDB when available; else fetches from RCSB by UniProt ID for any gene.
    """
    pdb_ids = PROTEIN_PDB.get(gene.upper(), [])
    if not pdb_ids and (uniprot_id or GENE_UNIPROT.get(gene.upper()) or _resolve_uniprot_id(gene)):
        uid = uniprot_id or GENE_UNIPROT.get(gene.upper()) or _resolve_uniprot_id(gene)
        pdb_ids = _fetch_pdb_by_uniprot(uid)
    chosen = pdb_ids[0] if pdb_ids else None
    has_res = True
    notes = "OK"
    mapping_method = "config"
    pdb_resseq = None

    if pos and pdb_ids:
        try:
            from visualization import fetch_pdb, residue_in_pdb
            from residue_mapper import extract_chain_residues, map_uniprot_to_pdb

            found = False
            for pid in pdb_ids:
                pd = fetch_pdb(pid)
                if not pd:
                    continue
                # Map UniProt pos -> PDB resseq when sequence available
                if uniprot_seq:
                    residues = extract_chain_residues(pd, "A")
                    mapping = map_uniprot_to_pdb(uniprot_seq, residues)
                    pdb_resseq = mapping.get(pos)
                    if pdb_resseq is not None and residue_in_pdb(pd, pdb_resseq, "A"):
                        chosen = pid
                        has_res = True
                        found = True
                        mapping_method = "alignment"
                        notes = "OK"
                        break
                else:
                    # Fallback: assume UniProt pos == PDB resseq
                    if residue_in_pdb(pd, pos, "A"):
                        chosen = pid
                        has_res = True
                        found = True
                        pdb_resseq = pos
                        mapping_method = "assumed_1:1"
                        break

            if not found:
                chosen = pdb_ids[0]
                has_res = False
                if uniprot_seq and pdb_resseq is None:
                    notes = "No UniProt↔PDB mapping; residue not found. Recommend AlphaFold."
                else:
                    notes = "Residue not present/resolved in available PDBs; recommend AlphaFold."
        except ImportError:
            # residue_mapper not available — fall back to raw position
            from visualization import fetch_pdb, residue_in_pdb
            for pid in pdb_ids:
                pd = fetch_pdb(pid)
                if pd and residue_in_pdb(pd, pos, "A"):
                    chosen, has_res, found = pid, True, True
                    break
            else:
                has_res, found = False, False
                notes = "Could not verify residue (mapping unavailable); recommend AlphaFold."
        except Exception:
            has_res = False
            notes = "Could not verify residue in PDBs; recommend AlphaFold or check mapping."

    return {
        "pdb_id": chosen,
        "chain": "A",
        "has_residue_coordinates": has_res if pos else True,
        "mapping_method": mapping_method,
        "pdb_ids": pdb_ids,
        "pdb_resseq": pdb_resseq,
        "notes": notes,
    }


def _heuristic_ddg(wt_aa: str, mut_aa: str) -> float:
    charge = {"R": 1, "K": 1, "D": -1, "E": -1, "H": 0.5}
    if abs(charge.get(wt_aa, 0) - charge.get(mut_aa, 0)) > 0.5:
        return 2.0
    h = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
         "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
         "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
    return 1.2 if abs(h.get(wt_aa, 0) - h.get(mut_aa, 0)) > 2 else 0.3


def get_ppi_ddg_predictions(
    gene: str,
    wt_aa: str,
    mut_aa: str,
    position: int,
    interactors: list,
    canonical_seq: Optional[str] = None,
) -> list:
    """PPI Delta-Delta-G per partner. NA when no complex or not at interface. Uses residue mapping if canonical_seq provided."""
    rows = []
    for ip in interactors:
        partner = ip.get("partner", "?")
        label = ip.get("label", partner)
        role = ip.get("role", "")
        key = f"{gene.upper()}_{partner}"
        pdb_id = PPI_PDB_COMPLEXES.get(key) if isinstance(PPI_PDB_COMPLEXES, dict) else None
        if isinstance(pdb_id, dict):
            pdb_id = pdb_id.get("pdb_id")

        has_c = bool(pdb_id and isinstance(pdb_id, str) and len(pdb_id) == 4)
        interface_info = None
        ppi_ddg = None
        method = "NA_no_complex"
        conclusion = "Neutral/Unclear"
        notes = "No complex model available; interface-specific ΔΔG not computed."
        confidence = "—"
        got_pdb = False
        pdb_resseq = None

        if has_c:
            try:
                from interface import is_residue_at_interface
                from visualization import fetch_pdb
                pd = fetch_pdb(str(pdb_id))
                got_pdb = bool(pd)
                if pd:
                    # Residue mapping: UniProt pos -> PDB resseq
                    if canonical_seq:
                        try:
                            from residue_mapper import map_uniprot_pos_to_pdb_resseq
                            mapping = map_uniprot_pos_to_pdb_resseq(pd, "A", canonical_seq)
                            pdb_resseq = mapping.get(position)
                            if pdb_resseq is None:
                                method = "NA_no_mapping"
                                notes = "No reliable residue mapping (UniProt↔PDB)."
                            else:
                                ir = is_residue_at_interface(pd, "A", "B", int(pdb_resseq), 8.0)
                                interface_info = ir
                                if ir.get("is_interface"):
                                    ppi_ddg = round(_heuristic_ddg(wt_aa, mut_aa), 2)
                                    method = "Heuristic_interface_proxy"
                                    confidence = "Low"
                                    conclusion = "LikelyDisruptive" if ppi_ddg > 1.0 else ("LikelyStabilizing" if ppi_ddg < -0.5 else "Neutral/Unclear")
                                    notes = "OK"
                                else:
                                    notes = f"Not interface (min_dist={ir.get('min_distance', 0):.1f} A)"
                        except ImportError:
                            ir = is_residue_at_interface(pd, "A", "B", position, 8.0)
                            interface_info = ir
                            if ir.get("is_interface"):
                                ppi_ddg = round(_heuristic_ddg(wt_aa, mut_aa), 2)
                                method = "Heuristic_interface_proxy"
                                confidence = "Low"
                                conclusion = "LikelyDisruptive" if ppi_ddg > 1.0 else ("LikelyStabilizing" if ppi_ddg < -0.5 else "Neutral/Unclear")
                                notes = "OK (no mapping; used UniProt pos)"
                            else:
                                notes = f"Not interface (min_dist={ir.get('min_distance', 0):.1f} A)"
                    else:
                        ir = is_residue_at_interface(pd, "A", "B", position, 8.0)
                        interface_info = ir
                        if ir.get("is_interface"):
                            ppi_ddg = round(_heuristic_ddg(wt_aa, mut_aa), 2)
                            method = "Heuristic_interface_proxy"
                            confidence = "Low"
                            conclusion = "LikelyDisruptive" if ppi_ddg > 1.0 else ("LikelyStabilizing" if ppi_ddg < -0.5 else "Neutral/Unclear")
                            notes = "OK (no canonical seq for mapping)"
                        else:
                            notes = f"Not interface (min_dist={ir.get('min_distance', 0):.1f} A)"
                else:
                    notes = "Complex model unavailable."
            except Exception:
                notes = "Complex unavailable."

        mapping_conf = "—"
        if has_c and got_pdb and canonical_seq:
            mapping_conf = "mapped" if pdb_resseq is not None else "not mapped"

        rows.append({
            "Interacting protein": label,
            "Role": role,
            "Tissue evidence": "Present",
            "Complex model available?": "Y" if (has_c and got_pdb) else "N",
            "Mapped PDB resseq": pdb_resseq if (has_c and pdb_resseq is not None) else "—",
            "Mapping confidence": mapping_conf,
            "Interface residue?": "Y" if (interface_info and interface_info.get("is_interface")) else ("N" if interface_info else "—"),
            "PPI Delta-Delta-G (kcal/mol)": ppi_ddg if ppi_ddg is not None else "NA",
            "Method": method,
            "Confidence": confidence,
            "Conclusion": conclusion,
            "Notes": notes,
        })
    return rows


def estimate_structural_impact(
    wt_aa: str,
    mut_aa: str,
    position: int,
    in_voltage_sensor: bool = False,
    af_deltas: Optional[dict] = None,
) -> dict:
    charge_change = {"R": 1, "K": 1, "D": -1, "E": -1, "H": 0.5}
    hydrophobicity = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
                     "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
                     "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
                     "Y": -1.3, "V": 4.2}
    c_wt, c_mut = charge_change.get(wt_aa, 0), charge_change.get(mut_aa, 0)
    h_wt, h_mut = hydrophobicity.get(wt_aa, 0), hydrophobicity.get(mut_aa, 0)
    impact = "Low"
    reasons = []

    # AlphaFold evidence augments impact
    if af_deltas:
        region = af_deltas.get("region_type", "")
        if region == "LowConfidence/IDR-like":
            impact = "Uncertain"
            reasons.append("Region low-confidence (IDR-like); structural impact uncertain.")
            return {"impact": impact, "reasons": reasons}
        delta_plddt = af_deltas.get("delta_mean_plddt_window")
        if delta_plddt is not None:
            if delta_plddt <= -10:
                impact = "High" if impact != "High" else impact
                reasons.append(f"AlphaFold pLDDT drop: Δ{delta_plddt:.1f}")
            elif delta_plddt <= -5:
                if impact == "Low":
                    impact = "Medium"
                reasons.append(f"AlphaFold pLDDT drop: Δ{delta_plddt:.1f}")

    # Physicochemical (never claim High from hydrophobicity alone; require charge or AF)
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
