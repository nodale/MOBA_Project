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
if(len(sys.argv) == 3):
    case_path = str(sys.argv[1]) + "/"
    St = float(sys.argv[2])
    print("Using case path: " + case_path)
    print("Using Strouhal number: " + str(St))
elif (len(sys.argv) > 1):
    print("ERROR: Enter case path and Strouhal number!")
    print("python visualize_optimal_forcing_response.py case_path Strouhal_number")
    exit()
else:
    case_path = "40/"
    # Strouhal number
    St = 0.1
    print("No case path provided. Using default case path: " + case_path)
    print("Using default Strouhal number: " + str(St))



settings = Settings(mesh_folder     = "mesh_full",
                    baseflow_folder = case_path,
                    mesh_name       = "cylinder_full",
                    result_folder   = case_path,
                    result_name     = "result",
                   )

# Load coefficients of base flow.
coefficients = np.load(settings.baseflow_folder + 'coefficients.npz')
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

# Load optimal forcing and response.
f_np = np.load(case_path + 'fq_' + str("{:4f}".format(St)) + '.npz')["f"]  # Optimal forcing.
q_np = np.load(case_path + 'fq_' + str("{:4f}".format(St)) + '.npz')["q"]  # Optimal response.

equation.q_real_list[equation.dof["u"]].vector().set_local(np.real(q_np[equation.VMixed.sub(equation.dof["u"]).dofmap().dofs()]))
equation.q_real_list[equation.dof["p"]].vector().set_local(np.real(q_np[equation.VMixed.sub(equation.dof["p"]).dofmap().dofs()]))
equation.q_real_list[equation.dof["T"]].vector().set_local(np.real(q_np[equation.VMixed.sub(equation.dof["T"]).dofmap().dofs()]))

equation.q_imag_list[equation.dof["u"]].vector().set_local(np.imag(q_np[equation.VMixed.sub(equation.dof["u"]).dofmap().dofs()]))
equation.q_imag_list[equation.dof["p"]].vector().set_local(np.imag(q_np[equation.VMixed.sub(equation.dof["p"]).dofmap().dofs()]))
equation.q_imag_list[equation.dof["T"]].vector().set_local(np.imag(q_np[equation.VMixed.sub(equation.dof["T"]).dofmap().dofs()]))

equation.f_real_list[equation.dof["u"]].vector().set_local(np.real(f_np[equation.VMixed.sub(equation.dof["u"]).dofmap().dofs()]))
equation.f_real_list[equation.dof["p"]].vector().set_local(np.real(f_np[equation.VMixed.sub(equation.dof["p"]).dofmap().dofs()]))
equation.f_real_list[equation.dof["T"]].vector().set_local(np.real(f_np[equation.VMixed.sub(equation.dof["T"]).dofmap().dofs()]))

equation.f_imag_list[equation.dof["u"]].vector().set_local(np.imag(f_np[equation.VMixed.sub(equation.dof["u"]).dofmap().dofs()]))

fields_to_write = {}
fields_to_write["u_real"]   = equation.q_real_list[equation.dof["u"]]
fields_to_write["u_imag"]   = equation.q_imag_list[equation.dof["u"]]
fields_to_write["p_real"]  = equation.q_real_list[equation.dof["p"]]
fields_to_write["p_imag"]  = equation.q_imag_list[equation.dof["p"]]
fields_to_write["T_real"]   = equation.q_real_list[equation.dof["T"]]
fields_to_write["T_imag"]   = equation.q_imag_list[equation.dof["T"]]
fields_to_write["uf_real"]  = equation.f_real_list[equation.dof["u"]]
fields_to_write["uf_imag"]  = equation.f_imag_list[equation.dof["u"]]

# Write response.
io = Io()
io.write_paraview(geometry, settings, "fq_" + str(St), fields_to_write)

print("Done!")