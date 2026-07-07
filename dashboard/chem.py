"""Chemical structure rendering (SMILES -> 2D SVG) and descriptors via RDKit."""

import math
import statistics

# Pale-green wash used to subtly mark the shared scaffold in a series depiction.
_SCAFFOLD_HL = (0.16, 0.77, 0.42, 0.22)

# Template mols (scaffold SMILES -> mol with fixed 2D coords), cached per process so a
# whole series aligns to one reference frame without recomputing coordinates per render.
_TEMPLATE_CACHE: dict = {}

# Fractional breathing margin added around a molecule's drawing window.
_FRAME_PAD = 0.08

# Wall-clock ceiling (seconds) for the maximum-common-substructure search used to
# align an activity-cliff pair. Cliff partners are similar by construction, so the
# search is fast; the cap only guards against a pathological pair stalling a render.
_MCS_TIMEOUT = 5

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


def _extents_from(mol, cx: float, cy: float) -> tuple:
    """Return (radius, half_width, half_height) of ``mol`` measured from (cx, cy).

    ``radius`` is the farthest atom (used for outlier detection); ``half_width`` and
    ``half_height`` are the largest |dx| and |dy| from the centre, used to size a window
    that actually fills the canvas rather than a square sized to the diagonal.
    """
    conf = mol.GetConformer()
    dxs = [abs(conf.GetAtomPosition(i).x - cx) for i in range(mol.GetNumAtoms())]
    dys = [abs(conf.GetAtomPosition(i).y - cy) for i in range(mol.GetNumAtoms())]
    radius = max(math.hypot(dx, dy) for dx, dy in zip(dxs, dys))
    return radius, max(dxs), max(dys)


def _largest_fragment(mol):
    """Largest fragment of a possibly multi-component (salt/solvent) mol.

    For depiction only: a cliff pair is compared on its parent structure, so a
    counter-ion or solvent sitting off to one side would otherwise inflate the shared
    drawing window and shrink the core we actually want to compare. Returns the input
    unchanged when it is already a single fragment.
    """
    from rdkit import Chem

    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) <= 1:
        return mol
    return max(frags, key=lambda m: m.GetNumHeavyAtoms())


def _apply_frame(drawer, mol, frame, width: int, height: int) -> None:
    """Centre and scale ``drawer`` on a shared drawing window (see ``smiles_to_svg``).

    Centres on the shared core centroid so the core is position-stable, sizes the window
    to the canvas aspect so the structure fills it, and - for a series - lets an outlier
    beyond ``fence`` use its own extents so it fits fully instead of clipping. A pair
    passes ``fence=None`` (no outlier to fence off), so both members share one window.
    """
    from rdkit.Geometry import Point2D

    cx, cy, frame_hx, frame_hy, fence = frame
    r_member, mhx, mhy = _extents_from(mol, cx, cy)
    base_hx, base_hy = (
        (mhx, mhy) if (fence is not None and r_member > fence) else (frame_hx, frame_hy)
    )
    base_hx, base_hy = max(base_hx, 1e-6), max(base_hy, 1e-6)
    # Grow the shorter axis to the canvas aspect so nothing is letterboxed.
    if base_hx / base_hy >= width / height:
        hx, hy = base_hx, base_hx * height / width
    else:
        hx, hy = base_hy * width / height, base_hy
    hx *= 1.0 + _FRAME_PAD
    hy *= 1.0 + _FRAME_PAD
    drawer.SetScale(width, height, Point2D(cx - hx, cy - hy), Point2D(cx + hx, cy + hy))


def _draw_svg(
    mol,
    width: int,
    height: int,
    frame: tuple | None,
    highlight_atoms: list,
    highlight_bonds: list,
    atom_colors: dict,
    bond_colors: dict,
) -> str:
    """Draw ``mol`` to an SVG, applying a shared frame and optional highlights."""
    from rdkit.Chem.Draw import rdMolDraw2D

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    if frame is not None:
        _apply_frame(drawer, mol, frame, width, height)
    drawer.DrawMolecule(
        mol,
        highlightAtoms=highlight_atoms or None,
        highlightBonds=highlight_bonds or None,
        highlightAtomColors=atom_colors or None,
        highlightBondColors=bond_colors or None,
    )
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def _core_highlight(mol, core) -> tuple:
    """Atoms/bonds of ``mol`` matching ``core``, tinted with the scaffold wash.

    Returns ``(atoms, bonds, atom_colors, bond_colors)`` for ``_draw_svg``; empty when
    the core does not match (so the depiction is simply left untinted).
    """
    atoms: list[int] = []
    bonds: list[int] = []
    atom_colors: dict = {}
    bond_colors: dict = {}
    match = mol.GetSubstructMatch(core)
    if match:
        match_set = set(match)
        for a in match:
            atoms.append(a)
            atom_colors[a] = _SCAFFOLD_HL
        for bond in mol.GetBonds():
            if bond.GetBeginAtomIdx() in match_set and bond.GetEndAtomIdx() in match_set:
                bonds.append(bond.GetIdx())
                bond_colors[bond.GetIdx()] = _SCAFFOLD_HL
    return atoms, bonds, atom_colors, bond_colors


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

    # Centre on the shared scaffold centroid so the core is position-stable across the
    # series, sizing the window to the canvas aspect so the structure fills it. Only
    # applied when the member actually aligned to the scaffold.
    return _draw_svg(
        mol,
        width,
        height,
        frame if (aligned and frame is not None) else None,
        highlight_atoms,
        highlight_bonds,
        atom_colors,
        bond_colors,
    )


def scaffold_frame(scaffold_smiles: str, member_smiles) -> tuple | None:
    """Shared centred drawing window for a scaffold's members.

    Returns ``(cx, cy, frame_hx, frame_hy, fence)`` where ``(cx, cy)`` is the scaffold
    centroid (identical across aligned members), ``frame_hx``/``frame_hy`` are the window
    half-extents that fit the typical members, and ``fence`` is the Tukey outlier
    threshold (Q3 + 1.5*IQR) on member radii - members beyond it are scaled to fit
    individually rather than shrinking the whole series. With fewer than four members no
    fence is set. Returns None when the scaffold is unusable or no member yields
    coordinates.
    """
    if not scaffold_smiles:
        return None
    from rdkit import Chem
    from rdkit.Chem import rdDepictor

    template = _template_mol(scaffold_smiles)
    if template is None:
        return None
    cx, cy = _template_centroid(template)

    extents: list[tuple] = []  # (radius, half_width, half_height) per member
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
        extents.append(_extents_from(mol, cx, cy))

    if not extents:
        return None

    radii = [r for r, _, _ in extents]
    fence = None
    if len(radii) >= 4:
        q1, _, q3 = statistics.quantiles(radii, n=4)
        fence = q3 + 1.5 * (q3 - q1)

    typical = [(hx, hy) for (r, hx, hy) in extents if fence is None or r <= fence]
    if not typical:
        typical = [(hx, hy) for (_, hx, hy) in extents]
    frame_hx = max(hx for hx, _ in typical)
    frame_hy = max(hy for _, hy in typical)
    return (cx, cy, frame_hx, frame_hy, fence)


def pair_core(smiles_a: str, smiles_b: str, strip_salts: bool = True):
    """Align two activity-cliff partners onto their maximum common substructure.

    Unlike a chemical series - where one Bemis-Murcko scaffold is a substructure of
    every member - the two compounds in a cliff pair need not share a scaffold at all
    (this is exactly what the mart's ``same_scaffold`` flag distinguishes). So the shared
    reference is their maximum common substructure (MCS), computed per pair, rather than
    a fixed scaffold. Cliff partners are similar by construction, so the MCS is large and
    captures precisely the common core; aligning both to it makes the differing
    substituents line up for a direct read.

    Molecule A keeps its own fresh 2D layout and becomes the reference frame; molecule B
    is redrawn onto A over the shared core, so both place the core at identical
    coordinates. Salts and solvents are stripped to the largest fragment first (depiction
    only). ``core`` is returned as an RDKit query mol.

    Returns ``(mol_a, mol_b, core)``, or None when either SMILES is unparseable, the MCS
    search is cancelled, or no common core is found - so the caller can fall back to two
    independent depictions and never break the page on an odd pair.
    """
    if not smiles_a or not smiles_b:
        return None
    from rdkit import Chem
    from rdkit.Chem import rdDepictor, rdFMCS

    _ensure_coordgen()
    mol_a = Chem.MolFromSmiles(smiles_a)
    mol_b = Chem.MolFromSmiles(smiles_b)
    if mol_a is None or mol_b is None:
        return None
    if strip_salts:
        mol_a = _largest_fragment(mol_a)
        mol_b = _largest_fragment(mol_b)

    # completeRingsOnly + ringMatchesRingOnly keep the core chemically whole (no half
    # rings), which is what makes the aligned depiction read cleanly.
    mcs = rdFMCS.FindMCS(
        [mol_a, mol_b],
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
        bondCompare=rdFMCS.BondCompare.CompareOrderExact,
        timeout=_MCS_TIMEOUT,
    )
    if mcs.canceled or mcs.numAtoms == 0:
        return None
    core = Chem.MolFromSmarts(mcs.smartsString)
    if core is None:
        return None

    rdDepictor.Compute2DCoords(mol_a)
    try:
        # Hard constraint: B's core is pinned to A's coordinates, the rest built around it.
        rdDepictor.GenerateDepictionMatching2DStructure(mol_b, mol_a, refPatt=core)
    except (ValueError, RuntimeError):
        return None
    return mol_a, mol_b, core


def pair_frame(mol_a, mol_b, core) -> tuple | None:
    """Shared, centred drawing window for an aligned cliff pair.

    Centres on the shared-core centroid (identical in both mols after ``pair_core``) and
    sizes the window to the larger of the two molecules' extents, so both panels use one
    scale and the core sits at the same place and size in each - the differing
    substituents then compare directly. Returns the same ``(cx, cy, hx, hy, fence)`` tuple
    ``_draw_svg`` consumes, with ``fence=None`` (a two-member pair has no outlier to fence
    off). Returns None if the core no longer matches molecule A.
    """
    match_a = mol_a.GetSubstructMatch(core)
    if not match_a:
        return None
    conf_a = mol_a.GetConformer()
    cx = sum(conf_a.GetAtomPosition(i).x for i in match_a) / len(match_a)
    cy = sum(conf_a.GetAtomPosition(i).y for i in match_a) / len(match_a)
    hx = hy = 0.0
    for mol in (mol_a, mol_b):
        _, mhx, mhy = _extents_from(mol, cx, cy)
        hx, hy = max(hx, mhx), max(hy, mhy)
    return (cx, cy, hx, hy, None)


def render_pair(
    smiles_a: str,
    smiles_b: str,
    width: int = 360,
    height: int = 300,
    align: bool = True,
    highlight_core: bool = False,
    strip_salts: bool = True,
) -> tuple:
    """Render an activity-cliff pair as two SVGs, aligned on their shared core.

    With ``align`` set, both structures are oriented so their maximum common substructure
    sits in one fixed frame at the same scale, so the substituent changes line up side by
    side; ``highlight_core`` washes that shared core. Falls back to two independent
    depictions when ``align`` is off or the pair cannot be aligned (no common core,
    unparseable input). Returns ``(svg_a, svg_b)``; an element is None if that SMILES is
    unusable.
    """
    plain = (
        smiles_to_svg(smiles_a, width=width, height=height),
        smiles_to_svg(smiles_b, width=width, height=height),
    )
    if not align:
        return plain
    result = pair_core(smiles_a, smiles_b, strip_salts=strip_salts)
    if result is None:
        return plain
    mol_a, mol_b, core = result
    frame = pair_frame(mol_a, mol_b, core)
    svgs = []
    for mol in (mol_a, mol_b):
        atoms, bonds, atom_colors, bond_colors = (
            _core_highlight(mol, core) if highlight_core else ([], [], {}, {})
        )
        svgs.append(_draw_svg(mol, width, height, frame, atoms, bonds, atom_colors, bond_colors))
    return tuple(svgs)
