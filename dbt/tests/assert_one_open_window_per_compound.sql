-- SCD2 integrity: at any point in time each molecule has exactly one OPEN validity
-- window (dbt_valid_to is null). More than one open window per molecule_chembl_id
-- means a snapshot update failed to close the prior record.
select molecule_chembl_id, count(*) as open_windows
from {{ ref('compound_status') }}
where dbt_valid_to is null
group by molecule_chembl_id
having count(*) > 1
