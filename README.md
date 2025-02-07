# MOBA_Project
The final Project of "Computational Methods for Operator-Based Analysis" at TUM winter 2025.
The focus is on applying the learned methods, analysing a base flow on a backward-facing step.

## Terminal Comands

### How to upload git repo to MOBA server
Use the script *archive.sh* to compress the repository into a *.tar* file:

> *./archive.sh*

Upload the archive and unpack *.tar* file on server: 
> *tar -xvf archive.tar*

### Other commands

Unzipping a *.zip* file is simple:
> *unzip name.zip*

If this command is not installed on the VM, run it on a local WSL and compressit into a *.tar.gz* file:
> *tar -czvf newTarName.tar.gz folderName*

to *cd /mnt/c/Users/User/Downloads*.

Useful when downloading the git respository to run it on the chair's maschine.

Activate anaconda: 
> *source /opt/anaconda/bin/activate* 

(must show (base)-environment)

*export_matrices.py* and *visualize_optimal_forcing.py* are run from the **fenicsproject** 
environment: 
> *conda activate fenicsproject*,

*resolvent_analysis.py* is run from the **petsc** environment: 
> *conda activate petsc*,

Generate mesh from my_mesh.geo
> *gmsh -2 -format msh2 my_mesh.geo*

> *conda activate fenicsproject*

> *dolfin-convert my_mesh.msh my_mesh.xml*

## Project Report

The project report will be done with latex. I (tobi) will share the document with the team!
