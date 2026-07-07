-- Computed cellular potency must sit in a plausible pChEMBL-like band (roughly
-- 1 mM to 1 pM). Values outside flag a unit-conversion or bad-data problem.
select molecule_chembl_id, cell_line_chembl_id, median_p_cyto
from {{ ref('mart_compound_cytotoxicity') }}
where median_p_cyto < 2 or median_p_cyto > 13
