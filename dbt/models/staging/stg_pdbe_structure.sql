-- PDBe in_pdb resolution: one row per (ligand code, PDB entry).
-- The raw layer lands the ligand (het) code as xref_id; rename to ligand_code
-- so downstream joins read clearly against the compound cross-references.
select
    xref_id  as ligand_code,
    pdb_id
from {{ source('raw', 'pdbe_structures') }}
where xref_id is not null
  and pdb_id is not null
