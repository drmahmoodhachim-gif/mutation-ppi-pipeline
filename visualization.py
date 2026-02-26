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
        r = requests.get(url, timeout=10)
        if r.ok:
            return r.text
    except Exception:
        pass
    return None


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


def render_py3dmol_html(pdb_data: str, residue_pos: Optional[int] = None, width: int = 800, height: int = 600) -> str:
    """Generate embeddable HTML for py3Dmol 3D structure viewer."""
    if not HAS_PY3DMOL or not pdb_data:
        return ""
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_data, "pdb")
    view.setStyle({"cartoon": {"colorscheme": "spectrum"}})
    if residue_pos is not None:
        view.setStyle({"resi": str(residue_pos)}, {"stick": {"colorscheme": "whiteCarbon"}, "cartoon": {"color": "red"}})
        view.addLabel(f"Mutant {residue_pos}", {"fontColor": "black", "fontSize": 12, "backgroundColor": "white"}, {"resi": str(residue_pos)})
    view.zoomTo()
    return view.write_html()
