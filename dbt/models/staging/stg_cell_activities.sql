-- Cellular cytotoxicity (grain: one row per cellular activity of a reference payload).
-- Cellular GI50/IC50 carry a raw standard_value in a concentration unit but no
-- pchembl_value, so potency (p_cyto = -log10(molar)) is computed here. Only exact
-- ('=') concentration readouts in known units are kept.
with source as (
    select * from {{ source('raw', 'cell_activities') }}
),

typed as (
    select
        cast(activity_id as bigint)      as activity_id,
        molecule_chembl_id,
        molecule_pref_name,
        canonical_smiles,
        payload_class,
        reference_name,
        target_chembl_id                 as cell_line_chembl_id,
        target_pref_name                 as cell_line,
        assay_chembl_id,
        standard_type,
        standard_relation,
        cast(standard_value as double)   as standard_value,
        standard_units
    from source
    where standard_relation = '='
      and standard_type in ('GI50', 'IC50')
      and standard_units in ('nM', 'uM', 'M', 'pM')
      and try_cast(standard_value as double) > 0
),

with_potency as (
    select
        *,
        -log10(
            standard_value * case standard_units
                when 'M'  then 1.0
                when 'uM' then 1e-6
                when 'nM' then 1e-9
                when 'pM' then 1e-12
            end
        ) as p_cyto
    from typed
)

select *
from with_potency
qualify row_number() over (partition by activity_id order by activity_id) = 1
