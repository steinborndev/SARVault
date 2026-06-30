-- Fails if any staged activity has a pchembl_value outside the plausible 0-14 range.
select
    activity_id,
    pchembl_value
from {{ ref('stg_activities') }}
where pchembl_value < 0
   or pchembl_value > 14
