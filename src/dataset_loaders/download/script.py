import anndata as ad
import scanpy as sc
import pandas as pd

## VIASH START
par = {
  "input": "resources/datasets/SCP2167/raw_data/SCP2167",
  "output": "resources/datasets/SCP2167/raw_dataset.h5ad"
}
## VIASH END

expression = f"{par['input']}/expression/64265d4084cbeaef62ef36a9"

def read_typed_csv(path: str) -> pd.DataFrame:
  col_names = pd.read_csv(path, nrows=1).columns.tolist()
  col_types = pd.read_csv(path, skiprows=0, nrows=1).iloc[0].tolist()

  # Convert type names to appropriate Pandas dtypes
  dtype_mapping = {
      "group": "category",
      "category": "category",
      "TYPE": "string",
      "numeric": "float64",
  }
  col_types = [dtype_mapping.get(t, t) for t in col_types]  

  # Read the rest of the CSV, applying the column names and types
  df = pd.read_csv(
    path,
    skiprows=2,
    names=col_names,
    dtype=dict(zip(col_names, col_types)),
    index_col="NAME"
  )

  return df


adata = sc.read_10x_mtx(expression)
adata
# AnnData object with n_obs × n_vars = 14165 × 36601
#     var: 'gene_ids', 'feature_types'

cluster = read_typed_csv(f"{par['input']}/cluster/humancortex_cluster.csv")
cluster
#                             X          Y    cell_type
# NAME                                                 
# CTACATTCAGCTTTGA-1 -16.746614  -1.999025   Excitatory
# AACCTTTCACTGGATT-1 -19.254899 -17.433621   Excitatory
# CCATCACGTTAGTCGT-1 -16.703370  -1.978855   Excitatory
# CTCCTTTCAGACCATT-1 -16.329322  -1.227991   Excitatory
# GCAACCGCACCAAATC-1 -20.231160  -2.528652   Excitatory
# ...                       ...        ...          ...
# AGTTCGACAGACGCTC-1   4.417137 -22.180945  Endothelial
# CCTCATGTCCTTATCA-1   4.948904 -23.954339  Endothelial
# TACGCTCTCCATCTCG-1   4.965150 -24.211097  Endothelial
# ATTCCTAGTGGAAATT-1  16.381543   9.752836        Oligo
# GCATCGGAGCACCGTC-1  12.886445  11.456965        Oligo

# [4067 rows x 3 columns]

spatial = read_typed_csv(f"{par['input']}/cluster/humancortex_spatial.csv")
spatial
#                               X            Y    cell_type
# NAME                                                     
# CTACATTCAGCTTTGA-1  4982.182919  2715.304054   Excitatory
# AACCTTTCACTGGATT-1  3877.271881  1197.231379   Excitatory
# CCATCACGTTAGTCGT-1  2977.881263  4087.932046   Excitatory
# CTCCTTTCAGACCATT-1  4738.346031  2468.137565   Excitatory
# GCAACCGCACCAAATC-1  3010.716289  2978.643333   Excitatory
# ...                         ...          ...          ...
# AGTTCGACAGACGCTC-1  1274.691325  1424.097687  Endothelial
# CCTCATGTCCTTATCA-1  3541.376000  3306.541636  Endothelial
# TACGCTCTCCATCTCG-1  3164.842000  5547.069250  Endothelial
# ATTCCTAGTGGAAATT-1  2257.890000  1391.706675        Oligo
# GCATCGGAGCACCGTC-1  2657.236500  2220.441000        Oligo

# [4067 rows x 3 columns]


metadata = read_typed_csv(f"{par['input']}/metadata/humancortex_metadata.csv")
metadata
#                    biosample_id       donor_id         species  ... library_preparation_protocol__ontology_label     sex      cluster
# NAME                                                            ...                                                                  
# CTACATTCAGCTTTGA-1   cortex_rna  human_cortex1  NCBITaxon_9606  ...                                    10x 3' v3  female   Excitatory
# AACCTTTCACTGGATT-1   cortex_rna  human_cortex1  NCBITaxon_9606  ...                                    10x 3' v3  female   Excitatory
# CCATCACGTTAGTCGT-1   cortex_rna  human_cortex1  NCBITaxon_9606  ...                                    10x 3' v3  female   Excitatory
# CTCCTTTCAGACCATT-1   cortex_rna  human_cortex1  NCBITaxon_9606  ...                                    10x 3' v3  female   Excitatory
# GCAACCGCACCAAATC-1   cortex_rna  human_cortex1  NCBITaxon_9606  ...                                    10x 3' v3  female   Excitatory
# ...                         ...            ...             ...  ...                                          ...     ...          ...
# AGTTCGACAGACGCTC-1   cortex_rna  human_cortex1  NCBITaxon_9606  ...                                    10x 3' v3  female  Endothelial
# CCTCATGTCCTTATCA-1   cortex_rna  human_cortex1  NCBITaxon_9606  ...                                    10x 3' v3  female  Endothelial
# TACGCTCTCCATCTCG-1   cortex_rna  human_cortex1  NCBITaxon_9606  ...                                    10x 3' v3  female  Endothelial
# ATTCCTAGTGGAAATT-1   cortex_rna  human_cortex1  NCBITaxon_9606  ...                                    10x 3' v3  female        Oligo
# GCATCGGAGCACCGTC-1   cortex_rna  human_cortex1  NCBITaxon_9606  ...                                    10x 3' v3  female        Oligo

# [4067 rows x 12 columns]

# >>> metadata.columns
# Index(['biosample_id', 'donor_id', 'species', 'species__ontology_label',
#        'disease', 'disease__ontology_label', 'organ', 'organ__ontology_label',
#        'library_preparation_protocol',
#        'library_preparation_protocol__ontology_label', 'sex', 'cluster'],
#       dtype='object')

# find intersect of indices
index_intersect = adata.obs.index\
  .intersection(metadata.index)\
  .intersection(cluster.index)\
  .intersection(spatial.index)\
  .astype(str)

# copy data to 
output = adata[index_intersect].copy()
for metadata_col in metadata.columns:
  output.obs[metadata_col] = metadata.loc[index_intersect, metadata_col]

# output.obs = metadata.loc[index_intersect, metadata.columns].copy()
output.obs["cell_type"] = list(cluster.loc[index_intersect, "cell_type"].values)
output.obsm["X_umap"] = cluster.loc[index_intersect, ["X", "Y"]].values
output.obsm["spatial"] = spatial.loc[index_intersect, ["X", "Y"]].values

# AnnData object with n_obs × n_vars = 4065 × 36601
#     obs: 'biosample_id', 'donor_id', 'species', 'species__ontology_label', 'disease', 'disease__ontology_label', 'organ', 'organ__ontology_label', 'library_preparation_protocol', 'library_preparation_protocol__ontology_label', 'sex', 'cluster', 'cell_type'
#     var: 'gene_ids', 'feature_types'
#     obsm: 'X_umap', 'spatial'

output.write_h5ad(par["output"], compression="gzip")
