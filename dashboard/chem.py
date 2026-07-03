"""Chemical structure rendering (SMILES -> 2D SVG) and descriptors via RDKit."""

import math
import statistics

# Pale-green wash used to subtly mark the shared scaffold in a series depiction.
_SCAFFOLD_HL = (0.16, 0.77, 0.42, 0.22)

# Template mols (scaffold SMILES -> mol with fixed 2D coords), cached per process so a
# whole series aligns to one reference frame without recomputing coordinates per render.
_TEMPLATE_CACHE: dict = {}

# Small breathing margin (Angstrom) added around a molecule's drawing window.
_FRAME_PAD = 0.6

_COORDGEN_ENABLED = False


def _ensure_coordgen() -> None:
    """Prefer the CoordGen layout engine (idempotent, process-global).

    CoordGen produces cleaner 2D layouts for fused and heavily-substituted systems and,
    under a template, re-orients pendant groups to avoid the overlaps the classic
    depictor leaves behind - while still placing the matched core on the template.
    """
    global _COORDGEN_ENABLED
    if not _COORDGEN_ENABLED:
        from rdkit.Chem import rdDepictor

        rdDepictor.SetPreferCoordGen(True)
        _COORDGEN_ENABLED = True


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

    _ensure_coordgen()
    templ = Chem.MolFromSmiles(scaffold_smiles)
    if templ is not None:
        rdDepictor.Compute2DCoords(templ)
    _TEMPLATE_CACHE[scaffold_smiles] = templ
    return templ


def _template_centroid(template) -> tuple:
    conf = template.GetConformer()
    n = template.GetNumAtoms()
    cx = sum(conf.GetAtomPosition(i).x for i in range(n)) / n
    cy = sum(conf.GetAtomPosition(i).y for i in range(n)) / n
    return cx, cy


def _radius_from(mol, cx: float, cy: float) -> float:
    """Largest distance from (cx, cy) to any atom of ``mol`` (its conformer)."""
    conf = mol.GetConformer()
    return max(
        math.dist((conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y), (cx, cy))
        for i in range(mol.GetNumAtoms())
    )


def smiles_to_svg(
    smiles: str,
    width: int = 420,
    height: int = 340,
    highlight_smarts: str | None = None,
    scaffold_smiles: str | None = None,
    align_to_scaffold: bool = False,
    highlight_scaffold: bool = False,
    frame: tuple | None = None,
) -> str | None:
    """Render a SMILES to a 2D SVG.

    Optionally (a) highlights a SMARTS substructure (compound-library search), and/or
    (b) orients the depiction so a shared Bemis-Murcko ``scaffold_smiles`` sits in a
    fixed frame, so that clicking through a series keeps the common core still and only
    the substituents move, with the core optionally given a subtle highlight.

    When ``frame`` (from ``scaffold_frame``) is given and the member aligns, the drawing
    is centred on the shared scaffold centroid so the core sits at the same spot for
    every member. Typical members share one scale (the core is the same size too); an
    unusually large "outlier" member is scaled down to fit around that same centre
    instead of being clipped, so only its size changes, not the core's position.

    Falls back to a plain depiction when the SMARTS is invalid/unmatched or the scaffold
    does not substructure-match (rare aromaticity/kekulisation cases), so a single odd
    member never breaks the page.
    """
    if not smiles:
        return None
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Geometry import Point2D

    _ensure_coordgen()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    template = _template_mol(scaffold_smiles) if scaffold_smiles else None

    # Constrain the depiction to the scaffold's coordinates (shared reference frame).
    aligned = False
    if template is not None and align_to_scaffold:
        try:
            rdDepictor.GenerateDepictionMatching2DStructure(mol, template)
            aligned = True
        except (ValueError, RuntimeError):
            rdDepictor.Compute2DCoords(mol)  # member doesn't match, plain layout

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
    # Centre on the shared scaffold centroid so the core is position-stable across the
    # series. Typical members share one radius (stable core size); an outlier gets its
    # own radius so it fits fully instead of clipping - only its size differs.
    if aligned and frame is not None:
        cx, cy, frame_radius, fence = frame
        r_member = _radius_from(mol, cx, cy)
        r = r_member if (fence is not None and r_member > fence) else frame_radius
        r = max(r, r_member) + _FRAME_PAD
        drawer.SetScale(width, height, Point2D(cx - r, cy - r), Point2D(cx + r, cy + r))
    drawer.DrawMolecule(
        mol,
        highlightAtoms=highlight_atoms or None,
        highlightBonds=highlight_bonds or None,
        highlightAtomColors=atom_colors or None,
        highlightBondColors=bond_colors or None,
    )
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def scaffold_frame(scaffold_smiles: str, member_smiles) -> tuple | None:
    """Shared centred drawing window for a scaffold's members.

    Returns ``(cx, cy, frame_radius, fence)`` where ``(cx, cy)`` is the scaffold centroid
    (identical across aligned members), ``frame_radius`` is the window half-size that
    fits the typical members, and ``fence`` is the Tukey outlier threshold (Q3 + 1.5*IQR)
    on member radii - members beyond it are scaled to fit individually rather than
    shrinking the whole series. With fewer than four members no fence is set. Returns None
    when the scaffold is unusable or no member yields coordinates.
    """
    if not scaffold_smiles:
        return None
    from rdkit import Chem
    from rdkit.Chem import rdDepictor

    template = _template_mol(scaffold_smiles)
    if template is None:
        return None
    cx, cy = _template_centroid(template)

    radii: list[float] = []
    for smi in member_smiles:
        if not smi:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        try:
            rdDepictor.GenerateDepictionMatching2DStructure(mol, template)
        except (ValueError, RuntimeError):
            continue
        radii.append(_radius_from(mol, cx, cy))

    if not radii:
        return None

    fence = None
    if len(radii) >= 4:
        q1, _, q3 = statistics.quantiles(radii, n=4)
        fence = q3 + 1.5 * (q3 - q1)

    typical = [r for r in radii if fence is None or r <= fence]
    frame_radius = max(typical) if typical else max(radii)
    return (cx, cy, frame_radius, fence)
