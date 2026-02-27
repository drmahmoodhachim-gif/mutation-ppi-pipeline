"""Interactive 3D protein structure visualization."""

import requests
import streamlit as st
from typing import Optional

try:
    import py3Dmol
    HAS_PY3DMOL = True
except ImportError:
    HAS_PY3DMOL = False


def fetch_pdb(pdb_id: str) -> Optional[str]:
    """Fetch PDB structure from RCSB."""
    url = f"https://files.rcsb.org/view/{pdb_id}.pdb"
    try:
        r = requests.get(url, timeout=30)
        if r.ok:
            return r.text
    except Exception:
        pass
    return None


def residue_in_pdb(pdb_data: str, residue_pos: int, chain: str = "A") -> bool:
    """Check if residue position (PDB resseq) has coordinates in PDB. Accepts signed residue numbers."""
    if not pdb_data:
        return False
    for line in pdb_data.split("\n"):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if len(line) < 26 or line[21:22].strip() != chain:
            continue
        res_str = line[22:26].strip()
        try:
            res = int(res_str)
        except ValueError:
            continue
        if res == residue_pos:
            return True
    return False


def _extract_chain(pdb_data: str, chain: str) -> str:
    """Extract only ATOM/HETATM for chain — lighter 3D (Nav1.5 has 4 chains)."""
    out = [l for l in pdb_data.split("\n") if l.startswith(("ATOM", "HETATM")) and len(l) >= 22 and l[21:22].strip() == chain]
    return "\n".join(out) + "\nEND\n" if out else pdb_data


def create_mol_viewer(
    pdb_data: str,
    chain: str = "A",
    residue_pos: Optional[int] = None,
    highlight_color: str = "red",
    style: str = "cartoon",
) -> Optional[object]:
    """Create interactive py3Dmol viewer with optional residue highlight."""
    if not HAS_PY3DMOL or not pdb_data:
        return None

    view = py3Dmol.view(width=800, height=600)
    view.addModel(pdb_data, "pdb")
    view.setStyle({style: {"colorscheme": "spectrum"}})

    if residue_pos is not None:
        selector = {"resi": str(residue_pos), "chain": chain}
        view.setStyle(selector, {"stick": {"colorscheme": "whiteCarbon"}, "cartoon": {"color": highlight_color}})
        view.addLabel(
            f"Mutant position {residue_pos}",
            {"fontColor": "black", "fontSize": 14, "backgroundColor": "white"},
            selector,
        )
        view.zoomTo(selector)

    view.zoomTo()
    return view


def render_stmol(pdb_data: str, residue_pos: Optional[int] = None, height: int = 500) -> None:
    """Render structure in Streamlit using stmol (fallback)."""
    try:
        import stmol
        xyzview = create_mol_viewer(pdb_data, residue_pos=residue_pos)
        if xyzview:
            xyz_html = xyzview.write_html()
            stmol.showmol(xyzview, height=height)
    except Exception as e:
        st.warning(f"3D viewer unavailable: {e}. Install: pip install stmol py3Dmol")


def render_py3dmol_html(
    pdb_data: str,
    residue_pos: Optional[int] = None,
    highlight_color: str = "red",
    label: Optional[str] = None,
    width: int = 600,
    height: int = 500,
) -> str:
    """Generate embeddable HTML for py3Dmol 3D structure viewer."""
    if not HAS_PY3DMOL or not pdb_data:
        return ""
    pdb_a = _extract_chain(pdb_data, "A")
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_a, "pdb")
    view.setStyle({"cartoon": {"colorscheme": "spectrum", "opacity": 0.85}})
    if residue_pos is not None:
        sel = {"resi": str(residue_pos), "chain": "A"}
        view.setStyle(sel, {
            "stick": {"color": highlight_color, "radius": 0.3},
            "sphere": {"scale": 0.6, "color": highlight_color},
            "cartoon": {"color": highlight_color},
        })
        lbl = label or f"Position {residue_pos}"
        view.addLabel(lbl, {"fontColor": "white", "fontSize": 16, "backgroundColor": highlight_color}, sel)
    view.zoomTo()
    return view.write_html()
