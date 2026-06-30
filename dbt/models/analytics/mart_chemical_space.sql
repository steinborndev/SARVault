-- Chemical-space mart (grain: compound): physicochemical profile + potency
-- summary + approved-drug flag, for distribution and approved-vs-research views.
with potency as (
    select
        compound_key,
        max(median_pchembl) as best_pchembl,
        count(*)            as n_targets
    from {{ ref('mart_target_sar') }}
    group by compound_key
)

select
    c.compound_key,
    c.molecule_chembl_id,
    c.pref_name,
    c.mw_freebase,
    c.alogp,
    c.hba,
    c.hbd,
    c.psa,
    c.rotatable_bonds,
    c.num_ro5_violations,
    c.aromatic_rings,
    c.qed_weighted,
    c.is_approved_drug,
    p.best_pchembl,
    coalesce(p.n_targets, 0) as n_targets
from {{ ref('dim_compound') }} c
left join potency p on c.compound_key = p.compound_key
