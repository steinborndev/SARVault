"""Chemical structure rendering (SMILES -> 2D SVG) and descriptors via RDKit."""

# Pale-green wash used to subtly mark the shared scaffold in a series depiction.
_SCAFFOLD_HL = (0.16, 0.77, 0.42, 0.22)

# Template mols (scaffold SMILES -> mol with fixed 2D coords), cached per process so a
# whole series aligns to one reference frame without recomputing coordinates per render.
_TEMPLATE_CACHE: dict = {}


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


def _template_mol(scaffold_smiles: str):
    """Scaffold mol with fixed 2D coords, cached per SMILES (None if unparseable)."""
    if scaffold_smiles in _TEMPLATE_CACHE:
        return _TEMPLATE_CACHE[scaffold_smiles]
    from rdkit import Chem
    from rdkit.Chem import rdDepictor

    templ = Chem.MolFromSmiles(scaffold_smiles)
    if templ is not None:
        rdDepictor.Compute2DCoords(templ)
    _TEMPLATE_CACHE[scaffold_smiles] = templ
    return templ


def smiles_to_svg(
    smiles: str,
    width: int = 420,
    height: int = 340,
    highlight_smarts: str | None = None,
    scaffold_smiles: str | None = None,
    align_to_scaffold: bool = False,
    highlight_scaffold: bool = False,
) -> str | None:
    """Render a SMILES to a 2D SVG.

    Optionally (a) highlights a SMARTS substructure (compound-library search), and/or
    (b) orients the depiction so a shared Bemis-Murcko ``scaffold_smiles`` sits in a
    fixed frame — so that clicking through a series keeps the common core still and only
    the substituents move — with the core optionally given a subtle highlight.

    Falls back to a plain depiction when the SMARTS is invalid/unmatched or the scaffold
    does not substructure-match (rare aromaticity/kekulisation cases), so a single odd
    member never breaks the page.
    """
    if not smiles:
        return None
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    template = _template_mol(scaffold_smiles) if scaffold_smiles else None

    # Constrain the depiction to the scaffold's coordinates (shared reference frame).
    if template is not None and align_to_scaffold:
        try:
            rdDepictor.GenerateDepictionMatching2DStructure(mol, template)
        except (ValueError, RuntimeError):
            rdDepictor.Compute2DCoords(mol)  # member doesn't match — plain layout

    highlight_atoms: list[int] = []
    highlight_bonds: list[int] = []
    atom_colors: dict = {}
    bond_colors: dict = {}

    # Subtle wash over the shared scaffold so substituent changes stand out.
    if template is not None and highlight_scaffold:
        match = mol.GetSubstructMatch(template)
        if match:
            match_set = set(match)
            for a in match:
                highlight_atoms.append(a)
                atom_colors[a] = _SCAFFOLD_HL
            for bond in mol.GetBonds():
                if bond.GetBeginAtomIdx() in match_set and bond.GetEndAtomIdx() in match_set:
                    highlight_bonds.append(bond.GetIdx())
                    bond_colors[bond.GetIdx()] = _SCAFFOLD_HL

    # SMARTS substructure highlight (compound-library search), kept independent.
    if highlight_smarts:
        query = Chem.MolFromSmarts(highlight_smarts)
        if query is not None:
            for m in mol.GetSubstructMatches(query):
                highlight_atoms.extend(m)

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(
        mol,
        highlightAtoms=highlight_atoms or None,
        highlightBonds=highlight_bonds or None,
        highlightAtomColors=atom_colors or None,
        highlightBondColors=bond_colors or None,
    )
    drawer.FinishDrawing()
    return drawer.GetDrawingText()
