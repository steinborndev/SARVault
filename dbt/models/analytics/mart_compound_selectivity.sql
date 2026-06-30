-- Selectivity mart (grain: compound). Pairs are summarized by median_pchembl;
-- the selectivity index is best minus second-best target median (log10 fold).
with pair as (
    select compound_key, molecule_chembl_id, target_chembl_id, median_pchembl
    from {{ ref('mart_target_sar') }}
),

agg as (
    select
        compound_key,
        molecule_chembl_id,
        count(*)            as n_targets,
        max(median_pchembl) as best_pchembl,
        min(median_pchembl) as worst_pchembl
    from pair
    group by compound_key, molecule_chembl_id
),

ranked as (
    select
        compound_key,
        target_chembl_id,
        median_pchembl,
        row_number() over (
            partition by compound_key
            order by median_pchembl desc, target_chembl_id
        ) as rn
    from pair
),

best as (
    select compound_key, target_chembl_id as best_target
    from ranked where rn = 1
),

second as (
    select compound_key, median_pchembl as second_best_pchembl
    from ranked where rn = 2
)

select
    a.compound_key,
    a.molecule_chembl_id,
    a.n_targets,
    a.best_pchembl,
    b.best_target,
    s.second_best_pchembl,
    case when a.n_targets >= 2 then a.best_pchembl - s.second_best_pchembl end as selectivity_index,
    a.best_pchembl - a.worst_pchembl as pchembl_spread
from agg a
join best b      on a.compound_key = b.compound_key
left join second s on a.compound_key = s.compound_key
