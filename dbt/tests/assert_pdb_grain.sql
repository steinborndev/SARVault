-- Fails if the compound x ligand_code x pdb_id grain is not unique.
select compound_key, ligand_code, pdb_id, count(*) as n
from {{ ref('mart_compound_pdb') }}
group by compound_key, ligand_code, pdb_id
having count(*) > 1
