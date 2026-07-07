-- Compound catalog (grain: compound). One row per measured compound with
-- identity, structure, physicochemical properties and potency summary.
with selectivity as (
    select * from {{ ref('mart_compound_selectivity') }}
),
compound as (
    select * from {{ ref('dim_compound') }}
),
pdb as (
    -- Co-crystal coverage per compound; LEFT-joined so compounds without a
    -- resolved PDB entry default to zero (keeps the catalog grain = compound).
    select compound_key, count(distinct pdb_id) as n_pdb_entries
    from {{ ref('mart_compound_pdb') }}
    group by compound_key
)

select
    s.compound_key,
    c.molecule_chembl_id,
    c.pref_name,
    c.canonical_smiles,
    c.inchi_key,
    c.mw_freebase,
    c.alogp,
    c.hba,
    c.hbd,
    c.psa,
    c.rotatable_bonds,
    c.num_ro5_violations,
    c.ro3_pass,
    c.aromatic_rings,
    c.qed_weighted,
    c.max_phase,
    c.is_approved_drug,
    c.molecule_type,
    s.n_targets,
    s.best_pchembl,
    s.best_target,
    t.payload_class,
    s.selectivity_index,
    coalesce(p.n_pdb_entries, 0)     as n_pdb_entries,
    coalesce(p.n_pdb_entries, 0) > 0 as has_pdb
from selectivity s
join compound c on s.compound_key = c.compound_key
left join pdb p on s.compound_key = p.compound_key
-- payload_class of the compound's most-potent (best) target; the catalog grain
-- stays one row per compound. Multi-class compounds are deferred to a bridge (F3.x).
left join {{ ref('dim_target') }} t on s.best_target = t.target_chembl_id
