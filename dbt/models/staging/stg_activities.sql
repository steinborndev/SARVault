-- Fact grain: one row per measured activity.
-- Filters: pchembl present, exact relation, whitelisted endpoint, no validity flag.
with source as (
    select * from {{ source('raw', 'activities') }}
),

filtered as (
    select
        cast(activity_id as bigint)      as activity_id,
        molecule_chembl_id,
        target_chembl_id,
        assay_chembl_id,
        standard_type,
        standard_relation,
        cast(standard_value as double)   as standard_value,
        standard_units,
        cast(pchembl_value as double)    as pchembl_value,
        document_chembl_id
    from source
    where pchembl_value is not null
      and standard_relation = '='
      and data_validity_comment is null
      and standard_type in (
        {%- for st in var('standard_types') %}'{{ st }}'{{ ", " if not loop.last }}{%- endfor %}
      )
)

select *
from filtered
qualify row_number() over (partition by activity_id order by activity_id) = 1
