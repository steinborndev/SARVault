-- Cytotoxicity mart (grain: reference payload x cell line). Per-compound-per-cell-line
-- cellular potency (p_cyto) for the reference ADC payloads, from the cellular pull
-- (F3.2). Self-contained: keyed on molecule_chembl_id with structure/name carried in
-- the source records, so it does not depend on dim_compound / dim_target.
with cell as (
    select * from {{ ref('stg_cell_activities') }}
)

select
    molecule_chembl_id,
    any_value(molecule_pref_name)  as molecule_pref_name,
    any_value(reference_name)      as reference_name,
    any_value(payload_class)       as payload_class,
    cell_line_chembl_id,
    any_value(cell_line)           as cell_line,
    count(*)                       as n_measurements,
    round(median(p_cyto), 3)       as median_p_cyto,
    round(max(p_cyto), 3)          as max_p_cyto
from cell
group by molecule_chembl_id, cell_line_chembl_id
