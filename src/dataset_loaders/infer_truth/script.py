import anndata as ad
import pandas as pd

## VIASH START
par = {
  "input": "resources/datasets/SCP2167/raw_dataset.h5ad",
  "output": "resources/datasets/SCP2167/dataset_with_inferred_truth.h5ad"
}
## VIASH END

# read the dataset
adata = ad.read_h5ad(par["input"])

# todo: add information
assumed_truth = pd.DataFrame({
  "source_cell_type": ...,
  "target_cell_type": ...,
  "ligand": ...,
  "receptor": ...,
  "colocalized": ...
})

adata.uns["assumed_truth"] = assumed_truth

# save the new dataset
adata.write_h5ad(par["output"], compression="gzip")