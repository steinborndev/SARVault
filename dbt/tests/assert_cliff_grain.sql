-- mart_activity_cliff grain: one row per (target, compound_a, compound_b), a<b.
-- Fails if any triple is duplicated, or any row violates the a<b ordering.
with dups as (
    select target_key, compound_key_a, compound_key_b
    from {{ ref('mart_activity_cliff') }}
    group by target_key, compound_key_a, compound_key_b
    having count(*) > 1
),
misordered as (
    select target_key, compound_key_a, compound_key_b
    from {{ ref('mart_activity_cliff') }}
    where compound_key_a >= compound_key_b
)
select * from dups
union all
select * from misordered
