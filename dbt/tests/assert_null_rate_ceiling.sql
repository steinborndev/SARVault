-- Anomaly guard: null-rate ceiling on a key column. A spike in missing canonical
-- SMILES means the structure-dependent layer (fingerprints, scaffolds, similarity)
-- would silently degrade, so a null rate above the tolerated ceiling fails loudly.
with rates as (
    select
        'dim_compound.canonical_smiles' as column_name,
        avg(case when canonical_smiles is null then 1.0 else 0.0 end) as null_rate
    from {{ ref('dim_compound') }}
)
select column_name, null_rate
from rates
where null_rate > {{ var('anomaly_max_null_rate', 0.20) }}
