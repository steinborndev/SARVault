-- Compound dimension: one row per molecule.
with molecules as (
    select * from {{ ref('stg_molecules') }}
)

select
    row_number() over (order by molecule_chembl_id) as compound_key,
    molecule_chembl_id,
    pref_name,
    canonical_smiles,
    inchi_key,
    mw_freebase,
    alogp,
    hba,
    hbd,
    psa,
    rotatable_bonds,
    num_ro5_violations,
    ro3_pass,
    aromatic_rings,
    qed_weighted,
    max_phase,
    is_approved_drug,
    molecule_type
from molecules
