-- Explode ChEMBL molecule cross_references (JSON) into one row per reference.
with molecules as (
    select molecule_chembl_id, cross_references
    from {{ source('raw', 'molecules') }}
    where cross_references is not null and cross_references != 'null'
),

exploded as (
    select
        molecule_chembl_id,
        unnest(json_extract(cross_references, '$[*]')) as ref
    from molecules
)

select
    molecule_chembl_id,
    lower(json_extract_string(ref, '$.xref_src')) as source,
    json_extract_string(ref, '$.xref_id')          as xref_id,
    json_extract_string(ref, '$.xref_name')        as xref_name
from exploded
where json_extract_string(ref, '$.xref_src') is not null
  and json_extract_string(ref, '$.xref_id') is not null
