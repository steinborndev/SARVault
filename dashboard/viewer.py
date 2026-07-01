"""Embedded 3D structure viewer (PDBe Mol*) for co-crystal PDB entries.

Builds a self-contained HTML fragment that loads the official PDBe Mol* web
component from a CDN and renders a single PDB entry. The structure data is
fetched client-side by the viewer from PDBe/RCSB, so this adds no server-side
dependency — it is a browser embed dropped into ``st.components.v1.html``.
"""

_PDBE_MOLSTAR_VERSION = "3.1.3"
_CDN = f"https://cdn.jsdelivr.net/npm/pdbe-molstar@{_PDBE_MOLSTAR_VERSION}/build"


def pdbe_molstar_html(pdb_id: str, height: int = 460) -> str:
    """Return an HTML embed of the PDBe Mol* viewer for one PDB entry.

    The background is set to match the dark dashboard theme (#0e1117) and water
    molecules are hidden so the bound ligand reads clearly against the protein.
    """
    pid = str(pdb_id).lower()
    return f"""
<link rel="stylesheet" href="{_CDN}/pdbe-molstar.css">
<script src="{_CDN}/pdbe-molstar-component.js"></script>
<div style="position:relative; width:100%; height:{height}px; border-radius:8px; overflow:hidden;">
  <pdbe-molstar
    molecule-id="{pid}"
    hide-water
    hide-expand-icon
    bg-color-r="14" bg-color-g="17" bg-color-b="23"
    style="position:absolute; inset:0; width:100%; height:100%;">
  </pdbe-molstar>
</div>
"""
