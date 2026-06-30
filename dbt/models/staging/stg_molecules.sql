-- One row per molecule, with flattened structures and physicochemical properties.
with source as (
    select * from {{ source('raw', 'molecules') }}
)

select
    molecule_chembl_id,
    pref_name,
    cast(max_phase as double)                                  as max_phase,
    coalesce(cast(max_phase as double) >= 4, false)            as is_approved_drug,
    molecule_type,
    "molecule_structures.canonical_smiles"                     as canonical_smiles,
    "molecule_structures.standard_inchi_key"                   as inchi_key,
    cast("molecule_properties.mw_freebase" as double)          as mw_freebase,
    cast("molecule_properties.alogp" as double)                as alogp,
    cast("molecule_properties.hba" as integer)                 as hba,
    cast("molecule_properties.hbd" as integer)                 as hbd,
    cast("molecule_properties.psa" as double)                  as psa,
    cast("molecule_properties.rtb" as integer)                 as rotatable_bonds,
    cast("molecule_properties.num_ro5_violations" as integer)  as num_ro5_violations,
    "molecule_properties.ro3_pass"                             as ro3_pass,
    cast("molecule_properties.aromatic_rings" as integer)      as aromatic_rings,
    cast("molecule_properties.qed_weighted" as double)         as qed_weighted
from source
qualify row_number() over (partition by molecule_chembl_id order by molecule_chembl_id) = 1
