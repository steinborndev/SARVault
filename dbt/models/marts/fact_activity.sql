-- Activity fact: one row per measured activity.
-- Inner joins to the dimensions; joining to dim_assay (confidence >= floor)
-- drops activities whose assay was filtered out in staging.
with activities as (
    select * from {{ ref('stg_activities') }}
),
compounds as (select compound_key, molecule_chembl_id from {{ ref('dim_compound') }}),
targets   as (select target_key,   target_chembl_id   from {{ ref('dim_target') }}),
assays    as (select assay_key,     assay_chembl_id    from {{ ref('dim_assay') }})

select
    a.activity_id,
    c.compound_key,
    t.target_key,
    s.assay_key,
    a.standard_type,
    a.standard_relation,
    a.standard_value,
    a.standard_units,
    a.pchembl_value,
    a.document_chembl_id
from activities a
join compounds c on a.molecule_chembl_id = c.molecule_chembl_id
join targets   t on a.target_chembl_id   = t.target_chembl_id
join assays    s on a.assay_chembl_id    = s.assay_chembl_id
