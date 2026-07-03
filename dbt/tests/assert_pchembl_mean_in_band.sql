-- Anomaly guard: distribution shift. pChEMBL values are -log10(molar potency) and for
-- this bioactivity slice cluster in a plausible band; a mean outside it flags a unit
-- error, a bad merge, or a corrupted extract rather than real chemistry.
with stats as (
    select avg(pchembl_value) as mean_pchembl from {{ ref('fact_activity') }}
)
select mean_pchembl
from stats
where mean_pchembl < {{ var('anomaly_pchembl_mean_low', 3.0) }}
   or mean_pchembl > {{ var('anomaly_pchembl_mean_high', 10.0) }}
