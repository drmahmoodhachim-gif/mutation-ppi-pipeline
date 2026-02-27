"""Interface detection for protein complexes."""


def is_residue_at_interface(
    complex_pdb: str,
    chain_a: str,
    chain_b: str,
    resi_a: int,
    cutoff_angstrom: float = 8.0,
) -> dict:
    """
    Check if residue resi_a in chain_a is at interface with chain_b.
    Returns: {
      'is_interface': bool,
      'min_distance': float,
      'n_contacts': int,
    }
    """
    if not complex_pdb or resi_a is None:
        return {"is_interface": False, "min_distance": 999.0, "n_contacts": 0}

    def parse_pdb_coords(pdb_text, chain_id, resi):
        coords = []
        for line in pdb_text.split("\n"):
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if len(line) < 54:
                continue
            ch = line[21:22].strip()
            try:
                res = int(line[22:26].strip().lstrip("-") or 0)
            except ValueError:
                continue
            if ch != chain_id or res != resi:
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append((x, y, z))
            except ValueError:
                continue
        return coords

    def parse_chain_coords(pdb_text, chain_id):
        coords = []
        for line in pdb_text.split("\n"):
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if len(line) < 54:
                continue
            if line[21:22].strip() != chain_id:
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append((x, y, z))
            except ValueError:
                continue
        return coords

    coords_a = parse_pdb_coords(complex_pdb, chain_a, resi_a)
    coords_b = parse_chain_coords(complex_pdb, chain_b)

    if not coords_a or not coords_b:
        return {"is_interface": False, "min_distance": 999.0, "n_contacts": 0}

    min_dist = 999.0
    n_contacts = 0
    for (xa, ya, za) in coords_a:
        for (xb, yb, zb) in coords_b:
            d = ((xa - xb) ** 2 + (ya - yb) ** 2 + (za - zb) ** 2) ** 0.5
            if d < min_dist:
                min_dist = d
            if d <= cutoff_angstrom:
                n_contacts += 1

    is_interface = min_dist <= cutoff_angstrom
    return {
        "is_interface": is_interface,
        "min_distance": round(min_dist, 2),
        "n_contacts": n_contacts,
    }
