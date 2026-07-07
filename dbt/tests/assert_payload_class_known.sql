-- Every in-scope target must carry a known ADC-payload mechanism class.
-- Guards against a target landing in the warehouse without a payload_class
-- (seed drift, or a new target added to config without a class).
select
    target_chembl_id,
    payload_class
from {{ ref('dim_target') }}
where payload_class is null
   or payload_class not in ('tubulin_inhibitor', 'topo1_inhibitor', 'topo2_inhibitor')
