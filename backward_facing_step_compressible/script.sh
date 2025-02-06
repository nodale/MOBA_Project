#!/bin/bash
eval "$(conda shell.bash hook)"
conda env list
source /opt/anaconda3/etc/profile.d/conda.sh
echo "---------------------------------------------------------"
for i in $(seq 100 100 1000); do
    echo "Starting process for Re = $i"    
    conda deactivate
    conda activate fenicsproject
    python export_matrices.py $i  #Step 1
    conda deactivate
    conda activate petsc
    python eigenvalue_solver.py $i #Step 2
    python resolvent_analysis.py $i #Step 3
    conda deactivate
    conda activate fenicsproject
    python visualize_eigenmodes.py $i #Step 4
    conda deactivate
    echo "Done for Re = $i"
done