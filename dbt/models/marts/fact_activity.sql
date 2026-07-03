-- Activity fact: one row per measured activity.
-- Inner joins to the dimensions; joining to dim_assay (confidence >= floor)
-- drops activities whose assay was filtered out in staging.
--
-- Materialised INCREMENTALLY on activity_id: the ChEMBL slice grows monotonically
-- across releases (activities are append-only, keyed by a stable activity_id), so a
-- refresh only needs to load activities not already present. On the first run the
-- relation is built in full; on later runs the is_incremental() branch restricts the
-- scan to new activity_ids, making a no-delta refresh a true no-op. delete+insert on
-- the unique_key keeps the load idempotent if an activity is ever re-emitted.
{{ config(
    materialized='incremental',
    unique_key='activity_id',
    incremental_strategy='delete+insert',
    on_schema_change='append_new_columns'
) }}

with activities as (
    select * from {{ ref('stg_activities') }}

    {% if is_incremental() %}
    -- Only activities not yet loaded (append-only source, stable activity_id).
    where activity_id not in (select activity_id from {{ this }})
    {% endif %}
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
