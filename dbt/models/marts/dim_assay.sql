-- Assay dimension: one row per assay at/above the confidence floor.
with assays as (
    select * from {{ ref('stg_assays') }}
)

select
    row_number() over (order by assay_chembl_id) as assay_key,
    assay_chembl_id,
    description,
    assay_type,
    confidence_score
from assays
