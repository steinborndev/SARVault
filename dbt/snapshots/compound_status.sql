-- Slowly-changing (Type 2) history of a compound's development status.
--
-- ChEMBL revises a molecule's max_phase as it progresses through the clinic, so the
-- same molecule_chembl_id can flip from a research compound to an approved drug across
-- releases. This snapshot records WHEN that happened: the `check` strategy opens a new
-- validity window whenever max_phase, is_approved_drug or pref_name changes, giving the
-- warehouse an auditable "research -> approved" timeline per payload.
--
-- Keyed on molecule_chembl_id (the ChEMBL natural key) rather than dim_compound's
-- surrogate compound_key, which is a row_number() and is NOT stable across releases.
{% snapshot compound_status %}
{{ config(
    unique_key='molecule_chembl_id',
    strategy='check',
    check_cols=['max_phase', 'is_approved_drug', 'pref_name'],
    invalidate_hard_deletes=True
) }}

select
    molecule_chembl_id,
    pref_name,
    max_phase,
    is_approved_drug
from {{ ref('stg_molecules') }}

{% endsnapshot %}
