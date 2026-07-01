-- Co-crystal PDB entries per compound (grain: compound x ligand_code x pdb_id).
-- Each compound's PDBe cross-reference is a chemical-component (het) code; the
-- PDBe in_pdb extract (stage 3) resolves that code to the PDB entries containing
-- it, so this mart exposes the concrete structures for the dashboard's viewer.
with pdbe_xref as (
    select molecule_chembl_id, xref_id as ligand_code
    from {{ ref('stg_compound_xref') }}
    where source = 'pdbe'
    union
    select molecule_chembl_id, xref_id as ligand_code
    from {{ ref('stg_compound_xref_unichem') }}
    where source = 'pdbe'
),
structures as (
    select ligand_code, pdb_id from {{ ref('stg_pdbe_structure') }}
),
compound as (
    select compound_key, molecule_chembl_id from {{ ref('dim_compound') }}
)

select distinct
    c.compound_key,
    c.molecule_chembl_id,
    x.ligand_code,
    s.pdb_id
from pdbe_xref x
join compound c on x.molecule_chembl_id = c.molecule_chembl_id
join structures s on x.ligand_code = s.ligand_code
