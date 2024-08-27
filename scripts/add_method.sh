#!/bin/bash

task_name="cell_cell_communication"
component_name="cellphonedb_max"
component_lang="python" # change this to "r" if need be

common/create_component/create_component \
  --task $task_name \
  --language "$component_lang" \
  --name "$component_name" \
  --api_file src/api/comp_method.yaml \
  --output "src/methods/$component_name"