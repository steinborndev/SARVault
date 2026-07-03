-- Anomaly guard: per-layer row-count floor. A near-empty layer signals a broken or
-- partial extract/build (e.g. an API that silently returned nothing) rather than a
-- genuinely tiny slice. Floors are dbt vars so production can raise them well above
-- the CI-fixture values.
with counts as (
    select 'stg_activities' as model, count(*) as n from {{ ref('stg_activities') }}
    union all
    select 'dim_compound',   count(*) from {{ ref('dim_compound') }}
    union all
    select 'fact_activity',  count(*) from {{ ref('fact_activity') }}
)
select model, n
from counts
where n < {{ var('anomaly_min_rows', 3) }}
