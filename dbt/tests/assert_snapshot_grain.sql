-- Snapshot grain: at most one window may open per (molecule, dbt_valid_from). A
-- duplicate (molecule_chembl_id, dbt_valid_from) pair would mean two windows opened
-- for the same compound at the same instant.
select molecule_chembl_id, dbt_valid_from, count(*) as n
from {{ ref('compound_status') }}
group by molecule_chembl_id, dbt_valid_from
having count(*) > 1
