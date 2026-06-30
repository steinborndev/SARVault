-- One row per assay, keeping only assays at or above the confidence floor.
with source as (
    select * from {{ source('raw', 'assays') }}
)

select
    assay_chembl_id,
    description,
    assay_type,
    cast(confidence_score as integer) as confidence_score,
    target_chembl_id
from source
where cast(confidence_score as integer) >= {{ var('min_confidence_score') }}
qualify row_number() over (partition by assay_chembl_id order by assay_chembl_id) = 1
