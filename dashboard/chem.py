"""Chemical structure rendering (SMILES -> 2D SVG) via RDKit."""


def smiles_to_svg(smiles: str, width: int = 420, height: int = 340) -> str | None:
    """Render a SMILES string to a 2D structure SVG, or None if it cannot be parsed."""
    if not smiles:
        return None
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()
