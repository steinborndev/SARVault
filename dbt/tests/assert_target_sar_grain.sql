-- Fails if the compound x target grain is not unique.
select compound_key, target_key, count(*) as n
from {{ ref('mart_target_sar') }}
group by compound_key, target_key
having count(*) > 1
