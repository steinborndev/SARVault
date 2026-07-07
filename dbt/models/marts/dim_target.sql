-- Target dimension: one row per target (human, in scope).
-- payload_class labels each target with its ADC-payload mechanism class. The
-- mapping mirrors config/target_set.yml as the payload_classes dbt var (see
-- dbt_project.yml); tests/test_payload_class.py guards the two against drift.
-- Sourced as an inline var (not a seed) so the dimension builds under a plain
-- `dbt run` (which does not load seeds), e.g. the incremental fact_activity path.
with targets as (
    select * from {{ ref('stg_targets') }}
),

payload_class as (
    select * from (
        values
        {%- for m in var('payload_classes') %}
            ('{{ m.target_chembl_id }}', '{{ m.payload_class }}'){{ "," if not loop.last }}
        {%- endfor %}
    ) as v(target_chembl_id, payload_class)
)

select
    row_number() over (order by t.target_chembl_id) as target_key,
    t.target_chembl_id,
    t.pref_name,
    t.target_type,
    t.organism,
    pc.payload_class
from targets t
left join payload_class pc on t.target_chembl_id = pc.target_chembl_id
