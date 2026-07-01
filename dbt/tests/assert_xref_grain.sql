-- Fails if the compound x source x xref_id grain is not unique.
select compound_key, source, xref_id, count(*) as n
from {{ ref('mart_compound_xref') }}
group by compound_key, source, xref_id
having count(*) > 1
