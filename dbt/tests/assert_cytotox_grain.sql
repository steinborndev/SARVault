-- One row per (reference payload, cell line) in the cytotoxicity mart.
select molecule_chembl_id, cell_line_chembl_id, count(*) as n
from {{ ref('mart_compound_cytotoxicity') }}
group by molecule_chembl_id, cell_line_chembl_id
having count(*) > 1
