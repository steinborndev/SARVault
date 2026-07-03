-- Every cliff must satisfy the configured floor (fails => a pair slipped the filter).
select *
from {{ ref('mart_activity_cliff') }}
where tanimoto < {{ var('cliff_min_tanimoto') }}
   or delta_pchembl < {{ var('cliff_min_delta_pchembl') }}
