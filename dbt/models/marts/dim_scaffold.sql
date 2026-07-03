-- Scaffold dimension: one row per distinct Bemis-Murcko scaffold in the set.
-- Acyclic compounds have no scaffold and are excluded here; their fingerprint
-- rows carry a null scaffold_key. n_compounds is a convenience rollup of members.
with cheminfo as (
    select molecule_chembl_id, murcko_scaffold_smiles, murcko_generic_smiles
    from {{ ref('stg_compound_cheminfo') }}
    where murcko_scaffold_smiles is not null
),

scaffolds as (
    select
        murcko_scaffold_smiles,
        any_value(murcko_generic_smiles)   as murcko_generic_smiles,
        count(distinct molecule_chembl_id) as n_compounds
    from cheminfo
    group by murcko_scaffold_smiles
)

select
    row_number() over (order by murcko_scaffold_smiles) as scaffold_key,
    murcko_scaffold_smiles,
    murcko_generic_smiles,
    n_compounds
from scaffolds
