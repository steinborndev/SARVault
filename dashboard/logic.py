"""Pure DataFrame helpers for scope filtering and landing-page metrics."""


def _target_keys(target_sar, targets):
    if not targets:
        return set(target_sar["compound_key"])
    return set(target_sar.loc[target_sar["target_pref_name"].isin(targets), "compound_key"])


def resolve_scope_keys(target_sar, catalog, scope):
    """Compound keys passing the scope's target / approval / min-potency facets."""
    scope = scope or {}
    keys = _target_keys(target_sar, scope.get("targets"))
    cat = catalog[catalog["compound_key"].isin(keys)]
    approval = scope.get("approval", "all")
    if approval == "approved":
        cat = cat[cat["is_approved_drug"]]
    elif approval == "research":
        cat = cat[~cat["is_approved_drug"]]
    min_p = scope.get("min_pchembl") or 0
    if min_p > 0:
        cat = cat[cat["best_pchembl"].fillna(-1) >= min_p]
    if scope.get("structure_only"):
        cat = cat[cat["has_pdb"].fillna(False)]
    return set(cat["compound_key"])


def scoped_target_sar(target_sar, scope, keys):
    """SAR pairs limited to in-scope compounds and (if set) selected targets."""
    scope = scope or {}
    df = target_sar[target_sar["compound_key"].isin(keys)]
    targets = scope.get("targets")
    if targets:
        df = df[df["target_pref_name"].isin(targets)]
    return df


_RO5_RULES = (
    ("MW ≤ 500", "mw_freebase", 500),
    ("logP ≤ 5", "alogp", 5),
    ("HBD ≤ 5", "hbd", 5),
    ("HBA ≤ 10", "hba", 10),
)


def _missing(value) -> bool:
    return value is None or (isinstance(value, float) and value != value)


def ro5_breakdown(row):
    """Per-criterion Lipinski Ro5 pass/fail plus the computed violation count.

    A criterion with a missing descriptor is reported as unknown (pass=None) and
    is not counted as a violation.
    """
    items = []
    violations = 0
    for label, col, threshold in _RO5_RULES:
        value = row.get(col)
        if _missing(value):
            passed = None
        else:
            passed = float(value) <= threshold
            if not passed:
                violations += 1
        items.append({"label": label, "value": value, "pass": passed})
    return {"items": items, "violations": violations}


def ligand_efficiency(pchembl, heavy_atoms):
    """Hopkins ligand efficiency (kcal/mol per heavy atom): 1.37 * pChEMBL / HAC.

    Returns None when potency or heavy-atom count is missing or non-positive.
    """
    if _missing(pchembl) or _missing(heavy_atoms) or not heavy_atoms:
        return None
    return 1.37 * float(pchembl) / float(heavy_atoms)


def lipophilic_efficiency(pchembl, logp):
    """Lipophilic ligand efficiency: pChEMBL - logP (None if either is missing)."""
    if _missing(pchembl) or _missing(logp):
        return None
    return float(pchembl) - float(logp)


def add_efficiency(df, hac_by_smiles):
    """Attach heavy_atoms, ligand_efficiency and lipophilic_efficiency columns.

    ``hac_by_smiles`` maps canonical_smiles -> heavy-atom count (None if
    unparseable). Requires columns: canonical_smiles, best_pchembl, alogp.
    """
    out = df.copy()
    out["heavy_atoms"] = out["canonical_smiles"].map(hac_by_smiles)
    out["ligand_efficiency"] = [
        ligand_efficiency(p, h) for p, h in zip(out["best_pchembl"], out["heavy_atoms"])
    ]
    out["lipophilic_efficiency"] = [
        lipophilic_efficiency(p, lp) for p, lp in zip(out["best_pchembl"], out["alogp"])
    ]
    return out


def overview_metrics(target_sar, catalog, scope):
    """Headline metrics for the landing page, restricted to the current scope."""
    keys = resolve_scope_keys(target_sar, catalog, scope)
    sar = scoped_target_sar(target_sar, scope, keys)
    cat = catalog[catalog["compound_key"].isin(keys)]
    return {
        "compounds": int(len(keys)),
        "activities": int(sar["n_measurements"].sum()),
        "targets": int(sar["target_pref_name"].nunique()),
        "pairs": int(len(sar)),
        "multi_target": int((cat["n_targets"] >= 2).sum()),
        "approved": int(cat["is_approved_drug"].sum()),
    }
