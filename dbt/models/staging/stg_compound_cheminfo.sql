-- One row per compound with its ECFP4 (Morgan r2, 2048-bit) fingerprint and
-- Bemis-Murcko scaffold, computed by the RDKit cheminfo stage (extract/cheminfo.py).
with source as (
    select * from {{ source('raw', 'compound_cheminfo') }}
)

select
    molecule_chembl_id,
    canonical_smiles,
    ecfp4_hex,
    cast(n_onbits as integer)         as n_onbits,
    cast(heavy_atom_count as integer) as heavy_atom_count,
    murcko_scaffold_smiles,
    murcko_generic_smiles
from source
qualify row_number() over (partition by molecule_chembl_id order by molecule_chembl_id) = 1
