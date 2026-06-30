-- One row per target, human only.
with source as (
    select * from {{ source('raw', 'targets') }}
)

select
    target_chembl_id,
    pref_name,
    organism,
    target_type,
    cast(tax_id as bigint) as tax_id
from source
where organism = 'Homo sapiens'
qualify row_number() over (partition by target_chembl_id order by target_chembl_id) = 1
