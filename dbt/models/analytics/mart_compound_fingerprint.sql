-- Fingerprint mart (grain: compound). ECFP4 + Murcko scaffold per compound, keyed
-- to dim_compound and (where the compound is cyclic) dim_scaffold. This is the
-- substrate for similarity search, substructure filtering, scaffold series and
-- activity-cliff analytics. Compounds with an unparseable/absent SMILES have no
-- row here (they were dropped by the RDKit stage).
with cheminfo as (
    select * from {{ ref('stg_compound_cheminfo') }}
),
compound as (
    select compound_key, molecule_chembl_id from {{ ref('dim_compound') }}
),
scaffold as (
    select scaffold_key, murcko_scaffold_smiles from {{ ref('dim_scaffold') }}
)

select
    c.compound_key,
    ci.molecule_chembl_id,
    ci.ecfp4_hex,
    ci.n_onbits,
    ci.heavy_atom_count,
    ci.murcko_scaffold_smiles,
    s.scaffold_key
from cheminfo ci
join compound c on ci.molecule_chembl_id = c.molecule_chembl_id
left join scaffold s on ci.murcko_scaffold_smiles = s.murcko_scaffold_smiles
