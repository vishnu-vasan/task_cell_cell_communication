#!/bin/bash

# download data
common/sync_resources/sync_resources \
  --input s3://openproblems-data/resources_test/cell_cell_communication \
  --output resources

# aws s3 cp resources/datasets/SCP2167/raw_dataset.h5ad \
#   s3://openproblems-data/resources_test/cell_cell_communication/datasets/SCP2167/raw_dataset.h5ad