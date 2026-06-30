-- Target dimension: one row per target (human, in scope).
with targets as (
    select * from {{ ref('stg_targets') }}
)

select
    row_number() over (order by target_chembl_id) as target_key,
    target_chembl_id,
    pref_name,
    target_type,
    organism
from targets
