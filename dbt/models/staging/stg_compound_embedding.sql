-- One row per compound with its 2-D UMAP coordinates (extract/embedding.py).
with source as (
    select * from {{ source('raw', 'compound_embedding') }}
)

select
    molecule_chembl_id,
    cast(umap_x as double) as umap_x,
    cast(umap_y as double) as umap_y
from source
qualify row_number() over (partition by molecule_chembl_id order by molecule_chembl_id) = 1
