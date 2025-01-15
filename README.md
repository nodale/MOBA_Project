# MOBA_Project
The final Project of "Computational Methods for Operator-Based Analysis" at TUM winter 2025.
The focus is on applying the learned methods, analysing a base flow on a backward-facing step.

## Terminal Comands

Upload the archive and unpack *.tar* file: 
> *tar -xvzf name.tar*

Activate anaconda: 
> *source /opt/anaconda/bin/activate* 

(must show (base)-environment)

*export_matrices.py* and *visualize_optimal_forcing.py* are run from the **fenicsproject** 
environment: 
> *conda activate fenicsproject*,

*resolvent_analysis.py* is run from the **petsc** environment: 
> *conda activate petsc*,

To switch environment you should deactivate the current one: 
> *conda deactivate*

Generate mesh from my_mesh.geo
> *gmsh -2 -format msh2 my_mesh.geo*

> *conda activate fenicsproject*

> *dolfin-convert my_mesh.msh my_mesh.xml*

## Project Report

The project report will be done with latex. I (tobi) will share the document with the team!