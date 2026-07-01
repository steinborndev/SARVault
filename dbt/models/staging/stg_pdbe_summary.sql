-- PDBe entry metadata: one row per PDB id with title, method, year, resolution.
-- Sourced from the PDBe summary + experiment endpoints during extract stage 3.
select
    lower(pdb_id) as pdb_id,
    title,
    method,
    try_cast(year as integer)        as year,
    try_cast(resolution as double)   as resolution
from {{ source('raw', 'pdbe_summary') }}
where pdb_id is not null
