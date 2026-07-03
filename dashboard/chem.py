"""Chemical structure rendering (SMILES -> 2D SVG) and descriptors via RDKit."""


def heavy_atom_count(smiles: str) -> int | None:
    """Number of non-hydrogen atoms in a molecule, or None if it cannot be parsed."""
    if not smiles:
        return None
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    return mol.GetNumHeavyAtoms() if mol is not None else None


def is_valid_smarts(smarts: str) -> bool:
    """True if ``smarts`` parses as a valid RDKit SMARTS query pattern."""
    if not smarts:
        return False
    from rdkit import Chem

    return Chem.MolFromSmarts(smarts) is not None


def has_substructure(smiles: str, smarts: str) -> bool | None:
    """Whether ``smiles`` contains the ``smarts`` substructure.

    Returns None when either the molecule or the query cannot be parsed, so callers
    can distinguish "no match" (False) from "unusable input" (None).
    """
    if not smiles or not smarts:
        return None
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    query = Chem.MolFromSmarts(smarts)
    if mol is None or query is None:
        return None
    return mol.HasSubstructMatch(query)


def smiles_to_svg(
    smiles: str, width: int = 420, height: int = 340, highlight_smarts: str | None = None
) -> str | None:
    """Render a SMILES to a 2D SVG, optionally highlighting a SMARTS substructure.

    Kept below the plain renderer so the highlight-aware version is the one imported
    elsewhere. Falls back to no highlight if the SMARTS is invalid or unmatched.
    """
    if not smiles:
        return None
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    highlight_atoms = []
    if highlight_smarts:
        query = Chem.MolFromSmarts(highlight_smarts)
        if query is not None:
            for match in mol.GetSubstructMatches(query):
                highlight_atoms.extend(match)

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms or None)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()
