-- SAR mart (grain: compound x target): potency summary per compound-target pair.
with fact as (
    select compound_key, target_key, assay_key, pchembl_value
    from {{ ref('fact_activity') }}
)

select
    f.compound_key,
    f.target_key,
    c.molecule_chembl_id,
    t.target_chembl_id,
    t.pref_name                       as target_pref_name,
    t.payload_class,
    median(f.pchembl_value)           as median_pchembl,
    max(f.pchembl_value)              as max_pchembl,
    count(*)                          as n_measurements,
    count(distinct f.assay_key)       as n_assays
from fact f
join {{ ref('dim_compound') }} c on f.compound_key = c.compound_key
join {{ ref('dim_target') }}   t on f.target_key   = t.target_key
group by f.compound_key, f.target_key, c.molecule_chembl_id, t.target_chembl_id, t.pref_name, t.payload_class
