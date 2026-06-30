-- Fails if the selectivity index is negative (best should be >= second-best).
select compound_key, selectivity_index
from {{ ref('mart_compound_selectivity') }}
where selectivity_index < 0
