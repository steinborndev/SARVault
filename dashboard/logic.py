"""Pure DataFrame helpers for scope filtering and landing-page metrics."""


def scoped_target_sar(target_sar, scope):
    """Restrict the SAR mart to the selected target names (None/empty = all)."""
    if not scope:
        return target_sar
    return target_sar[target_sar["target_pref_name"].isin(scope)]


def scope_compound_keys(target_sar, scope):
    """Compound keys measured against any target in scope."""
    return set(scoped_target_sar(target_sar, scope)["compound_key"])


def overview_metrics(target_sar, selectivity, chemical_space, scope):
    """Headline metrics for the landing page, restricted to the current scope."""
    sar = scoped_target_sar(target_sar, scope)
    keys = set(sar["compound_key"])
    sel = selectivity[selectivity["compound_key"].isin(keys)]
    chem = chemical_space[chemical_space["compound_key"].isin(keys)]
    return {
        "compounds": int(chem["compound_key"].nunique()),
        "activities": int(sar["n_measurements"].sum()),
        "targets": int(sar["target_pref_name"].nunique()),
        "pairs": int(len(sar)),
        "multi_target": int((sel["n_targets"] >= 2).sum()),
        "approved": int(chem["is_approved_drug"].sum()),
    }
