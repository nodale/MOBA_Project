################################################################################################
### Computational Methods for Operator-Based Analysis
#Technical University of Munich, Professorship for Thermofluid Dynamics | Pr. Polifke, Ph.D
#
#Created: 06/2024 | P. Brokof, G. Varillon
################################################################################################

# Include path to flow solver.
import os, sys
sys.path.append("/opt/dynX")

# Import packages.
import numpy as np
from equations.compressible.compressibleNSE import CompressibleNSE
from settings                 import Settings
from geometry                 import Geometry
from utilities.io             import Io
from dolfin                   import *

# Check anaconda environment. This script needs FEniCS loaded. Scalar type should be float64.
from petsc4py import PETSc
print("Scalar type: " + str(PETSc.ScalarType))

# Get re from system arguments.
if(len(sys.argv) == 2):
    case_path = str(sys.argv[1])
    print("Using case path: " + case_path)
else:
    case_path = "100"
    print("No case path provided. Using default case path: " + case_path)



settings = Settings(mesh_folder     = "mesh_ref",
                    baseflow_folder = case_path,
                    mesh_name       = "backward_facing_step",
                    result_folder   = case_path,
                    result_name     = "eigenmodes",
                   )

# Load coefficients of base flow.
coefficients = np.load(settings.baseflow_folder + '/coefficients.npz')
settings.coefficients = {"Re":    coefficients["Re"],
                         "Pe":    coefficients["Pe"],
                         "gamma": coefficients["gamma"],
                         "Ma":    coefficients["Ma"]}

# Make sure this flag is set if the base flow was computed with it!!!
parameters["reorder_dofs_serial"] = False

# Create discretization on mesh.
geometry = Geometry(settings)

# Setup equation.
equation  = CompressibleNSE(settings=settings, geometry=geometry)

EWs = np.load(case_path + '/spectrum.npz')['spectrum']

# Load eigenmodes.
for i in range(len(EWs)):
    evR_np = np.load(case_path + '/ev_' + str(i) + '.npz')['evR']
    evL_np = np.load(case_path + '/ev_' + str(i) + '.npz')['evL']

    equation.q_real_list[equation.dof["u"]].vector().set_local(np.real(evR_np[equation.VMixed.sub(equation.dof["u"]).dofmap().dofs()]))
    equation.q_imag_list[equation.dof["u"]].vector().set_local(np.imag(evR_np[equation.VMixed.sub(equation.dof["u"]).dofmap().dofs()]))

    equation.f_real_list[equation.dof["u"]].vector().set_local(np.real(evL_np[equation.VMixed.sub(equation.dof["u"]).dofmap().dofs()]))
    equation.f_imag_list[equation.dof["u"]].vector().set_local(np.imag(evL_np[equation.VMixed.sub(equation.dof["u"]).dofmap().dofs()]))

    fields_to_write = {}
    fields_to_write["evR_real"] = equation.q_real_list[equation.dof["u"]]
    fields_to_write["evR_imag"] = equation.q_imag_list[equation.dof["u"]]
    fields_to_write["evL_real"] = equation.f_real_list[equation.dof["u"]]
    fields_to_write["evL_imag"] = equation.f_imag_list[equation.dof["u"]]

    # Write eigenmodes.
    io = Io()
    io.write_paraview(geometry, settings, "eigenmode_" + str(i), fields_to_write)

print("Done!")
