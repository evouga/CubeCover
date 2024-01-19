#!/bin/bash

# Check if exactly three arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <mesh_path> <fra_path> <subd_mesh_path>"
    exit 1
fi

# Assign arguments to variables
mesh_path="$1"
fra_path="$2"
subd_mesh_path="$3"

# Call the MATLAB function with these arguments
matlab -nodisplay -nodesktop -r "addpath('./io/', './subdivide/', './mesh/'); subdivide_mesh('${mesh_path}', '${fra_path}', '${subd_mesh_path}'); exit"