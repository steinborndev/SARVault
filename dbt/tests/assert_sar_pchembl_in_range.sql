-- Fails if any aggregated potency falls outside the plausible 0-14 range.
select compound_key, target_key, median_pchembl, max_pchembl
from {{ ref('mart_target_sar') }}
where median_pchembl < 0 or median_pchembl > 14
   or max_pchembl < 0 or max_pchembl > 14
