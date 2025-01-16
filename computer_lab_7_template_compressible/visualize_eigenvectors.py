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

# Check anaconda environment. This script needs FEniCsS loaded. Scalar type should be float64.
from petsc4py import PETSc
print("Scalar type: " + str(PETSc.ScalarType))

# Set here the eigenpair to post-process.
iEv = 0

settings = Settings(mesh_folder     = "mesh_full",
                    baseflow_folder = "baseflow",
                    mesh_name       = "cylinder_full",
                    result_folder   = "result",
                    result_name     = "result",
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


# Load direct eigenvector.
evR_np = np.load(settings.result_folder + '/ev_' + str(iEv) + '.npz')["evR"]
evL_np = np.load(settings.result_folder + '/ev_' + str(iEv) + '.npz')["evL"]

equation.q_real_list[equation.dof["u"]].vector().set_local(np.real(evR_np[equation.VMixed.sub(equation.dof["u"]).dofmap().dofs()]))
equation.q_real_list[equation.dof["p"]].vector().set_local(np.real(evR_np[equation.VMixed.sub(equation.dof["p"]).dofmap().dofs()]))
equation.q_real_list[equation.dof["T"]].vector().set_local(np.real(evR_np[equation.VMixed.sub(equation.dof["T"]).dofmap().dofs()]))

equation.q_imag_list[equation.dof["u"]].vector().set_local(np.imag(evR_np[equation.VMixed.sub(equation.dof["u"]).dofmap().dofs()]))
equation.q_imag_list[equation.dof["p"]].vector().set_local(np.imag(evR_np[equation.VMixed.sub(equation.dof["p"]).dofmap().dofs()]))
equation.q_imag_list[equation.dof["T"]].vector().set_local(np.imag(evR_np[equation.VMixed.sub(equation.dof["T"]).dofmap().dofs()]))

# Load adjoint eigenvector.
equation.qH_real_list[equation.dof["u"]].vector().set_local(np.real(evL_np[equation.VMixed.sub(equation.dof["u"]).dofmap().dofs()]))
equation.qH_real_list[equation.dof["p"]].vector().set_local(np.real(evL_np[equation.VMixed.sub(equation.dof["p"]).dofmap().dofs()]))
equation.qH_real_list[equation.dof["T"]].vector().set_local(np.real(evL_np[equation.VMixed.sub(equation.dof["T"]).dofmap().dofs()]))

equation.qH_imag_list[equation.dof["u"]].vector().set_local(np.imag(evL_np[equation.VMixed.sub(equation.dof["u"]).dofmap().dofs()]))
equation.qH_imag_list[equation.dof["p"]].vector().set_local(np.imag(evL_np[equation.VMixed.sub(equation.dof["p"]).dofmap().dofs()]))
equation.qH_imag_list[equation.dof["T"]].vector().set_local(np.imag(evL_np[equation.VMixed.sub(equation.dof["T"]).dofmap().dofs()]))

fields_to_write = {}
fields_to_write["u_real"]   = equation.q_real_list[equation.dof["u"]]
fields_to_write["u_imag"]   = equation.q_imag_list[equation.dof["u"]]
fields_to_write["p_real"]   = equation.q_real_list[equation.dof["p"]]
fields_to_write["p_imag"]   = equation.q_imag_list[equation.dof["p"]]
fields_to_write["T_real"]   = equation.q_real_list[equation.dof["T"]]
fields_to_write["T_imag"]   = equation.q_imag_list[equation.dof["T"]]
fields_to_write["uH_real"]  = equation.qH_real_list[equation.dof["u"]]
fields_to_write["uH_imag"]  = equation.qH_imag_list[equation.dof["u"]]
fields_to_write["pH_real"]  = equation.qH_real_list[equation.dof["p"]]
fields_to_write["pH_imag"]  = equation.qH_imag_list[equation.dof["p"]]
fields_to_write["TH_real"]  = equation.qH_real_list[equation.dof["T"]]
fields_to_write["TH_imag"]  = equation.qH_imag_list[equation.dof["T"]]

# Write response.
io = Io()
io.write_paraview(geometry, settings, "mode=" + str(iEv), fields_to_write)

print("Done.")