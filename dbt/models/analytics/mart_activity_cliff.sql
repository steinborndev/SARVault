-- Activity-cliff mart (grain: target x compound_a x compound_b, with a < b).
--
-- An activity cliff is a pair of structurally similar compounds with a large
-- potency difference *on the same target* — the sharpest, most information-rich
-- signal in structure-activity analysis. This model self-joins the per-target
-- potency summary on the fingerprint mart, computes ECFP4 Tanimoto via on-bit
-- list intersection, and keeps pairs above the configured floor (vars:
-- cliff_min_tanimoto / cliff_min_delta_pchembl). The floor is deliberately
-- permissive; the dashboard filters and ranks (by SALI) up from here.
--
-- SALI = |dpChEMBL| / (1 - Tanimoto), the cliff steepness. For identical 2D
-- fingerprints (Tanimoto = 1) SALI is undefined (division by zero); such pairs
-- are flagged is_identical_fp and almost always reflect stereochemistry,
-- tautomers or replicate-measurement variance rather than a 2D-structural change.
with sar as (
    select compound_key, target_key, median_pchembl
    from {{ ref('mart_target_sar') }}
),

fp as (
    select compound_key, molecule_chembl_id, ecfp4_onbits, n_onbits,
           murcko_scaffold_smiles
    from {{ ref('mart_compound_fingerprint') }}
),

-- Candidate pairs within a target (a < b), pre-filtered by the cheap potency
-- delta before the fingerprint intersection is computed.
pairs as (
    select
        a.target_key,
        a.compound_key                       as compound_key_a,
        b.compound_key                       as compound_key_b,
        a.median_pchembl                      as pchembl_a,
        b.median_pchembl                      as pchembl_b,
        abs(a.median_pchembl - b.median_pchembl) as delta_pchembl
    from sar a
    join sar b
      on a.target_key = b.target_key
     and a.compound_key < b.compound_key
    where abs(a.median_pchembl - b.median_pchembl) >= {{ var('cliff_min_delta_pchembl') }}
),

scored as (
    select
        p.target_key,
        p.compound_key_a,
        p.compound_key_b,
        fa.molecule_chembl_id                as molecule_chembl_id_a,
        fb.molecule_chembl_id                as molecule_chembl_id_b,
        p.pchembl_a,
        p.pchembl_b,
        p.delta_pchembl,
        fa.murcko_scaffold_smiles            as scaffold_a,
        fb.murcko_scaffold_smiles            as scaffold_b,
        len(list_intersect(fa.ecfp4_onbits, fb.ecfp4_onbits))::double
            / (fa.n_onbits + fb.n_onbits
               - len(list_intersect(fa.ecfp4_onbits, fb.ecfp4_onbits))) as tanimoto
    from pairs p
    join fp fa on p.compound_key_a = fa.compound_key
    join fp fb on p.compound_key_b = fb.compound_key
)

select
    s.target_key,
    t.pref_name                              as target_pref_name,
    s.compound_key_a,
    s.compound_key_b,
    s.molecule_chembl_id_a,
    s.molecule_chembl_id_b,
    round(s.pchembl_a, 2)                    as pchembl_a,
    round(s.pchembl_b, 2)                    as pchembl_b,
    round(s.delta_pchembl, 2)                as delta_pchembl,
    round(s.tanimoto, 3)                     as tanimoto,
    s.tanimoto >= 1.0                        as is_identical_fp,
    s.scaffold_a = s.scaffold_b              as same_scaffold,
    case
        when s.tanimoto < 1.0
        then round(s.delta_pchembl / (1.0 - s.tanimoto), 2)
    end                                      as sali
from scored s
join {{ ref('dim_target') }} t on s.target_key = t.target_key
where s.tanimoto >= {{ var('cliff_min_tanimoto') }}
