-- Chemical-space mart (grain: compound): physicochemical profile + potency
-- summary + approved-drug flag, plus the 2-D UMAP embedding, best target and
-- Murcko scaffold, for the distribution, approved-vs-research and structural-
-- embedding views.
with potency as (
    select
        compound_key,
        max(median_pchembl) as best_pchembl,
        count(*)            as n_targets
    from {{ ref('mart_target_sar') }}
    group by compound_key
),

selectivity as (
    select compound_key, best_target from {{ ref('mart_compound_selectivity') }}
),

fingerprint as (
    select compound_key, murcko_scaffold_smiles from {{ ref('mart_compound_fingerprint') }}
),

embedding as (
    select molecule_chembl_id, umap_x, umap_y from {{ ref('stg_compound_embedding') }}
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
    coalesce(p.n_targets, 0) as n_targets,
    s.best_target,
    f.murcko_scaffold_smiles,
    e.umap_x,
    e.umap_y
from {{ ref('dim_compound') }} c
left join potency p     on c.compound_key = p.compound_key
left join selectivity s on c.compound_key = s.compound_key
left join fingerprint f on c.compound_key = f.compound_key
left join embedding e   on c.molecule_chembl_id = e.molecule_chembl_id
