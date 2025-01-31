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
from equations.compressible.compressibleNSE     import CompressibleNSE
from settings                                   import Settings
from fixpoint                                   import Fixpoint_solver
from geometry                                   import Geometry
from utilities.io                               import Io
from utilities.convert2csr                      import convert_2_csrmat

from dolfin import *
import numpy as np
import scipy
import os

# Check for correct anaconda environment. This scrip needs FEniCS loaded. Scalar type should be float64.
from petsc4py import PETSc
print("Scalar type: " + str(PETSc.ScalarType))

# Get re from system arguments.
if(len(sys.argv) > 1):
    Re = int(sys.argv[1])
    if Re != 100:
        restart_name = "baseflow"
        #check if Re - 100 dir exists
        if not os.path.exists(str(Re - 100)):
            print("ERROR: Directory " + str(Re - 100) + " does not exist. Please compute baseflow for Re = " + str(Re - 100) + " first.")
            exit()
    else:
        restart_name = "no_restart"
else:
    restart_name = "no_restart"
    Re = 100 # default

print("Using Re: " + str(Re))
    
try:
    os.mkdir(os.path.join(os.getcwd(), str(Re)))
except:
    pass

# Set boundary conditions.
class BackwardFacingStepNSE(CompressibleNSE):

    def set_boundary_conditions_nonlinear_opreator(self, geometry, settings):
        BCmethod = "geometric"
        # Gamma must be the same ratio as in the .geo file.
        U_expr = Expression("(x[1]-Gamma*2)*(2-x[1])/((1-Gamma)*(1-Gamma))", Gamma=0.5, degree=2)
        self.bcs_nonlinear_operator = [
            # inlet
            DirichletBC(self.VMixed.sub(self.dof["u"]).sub(0),  U_expr, geometry.boundaries, settings.boundary["inlet"],    BCmethod),
            DirichletBC(self.VMixed.sub(self.dof["u"]).sub(1),  0.0,    geometry.boundaries, settings.boundary["inlet"],    BCmethod),
            DirichletBC(self.VMixed.sub(self.dof["T"]),         1.0,    geometry.boundaries, settings.boundary["inlet"],    BCmethod),
            # walls (adiabatic & no-slip)
            DirichletBC(self.VMixed.sub(self.dof["u"]).sub(0),  0.0,    geometry.boundaries, settings.boundary["walls"], BCmethod),
            DirichletBC(self.VMixed.sub(self.dof["u"]).sub(1),  0.0,    geometry.boundaries, settings.boundary["walls"], BCmethod),
            # outlet
            # Stress free outflow imposed via Neumann boundary.
        ]
        self.set_backflow_boundaries(geometry.ds(settings.boundary["outlet"]))

    def set_boundary_conditions_linear_operator_homogeneous(self, geometry, settings):
        BCmethod = "geometric"
        self.bcs_linear_operator_homogeneous = [
            # inlet
            DirichletBC(self.VMixed.sub(self.dof["u"]).sub(0),  0.0,    geometry.boundaries, settings.boundary["inlet"],    BCmethod),
            DirichletBC(self.VMixed.sub(self.dof["u"]).sub(1),  0.0,    geometry.boundaries, settings.boundary["inlet"],    BCmethod),
            DirichletBC(self.VMixed.sub(self.dof["T"]),         0.0,    geometry.boundaries, settings.boundary["inlet"],    BCmethod),
            # walls (adiabatic & no-slip)
            DirichletBC(self.VMixed.sub(self.dof["u"]).sub(0),  0.0,    geometry.boundaries, settings.boundary["walls"], BCmethod),
            DirichletBC(self.VMixed.sub(self.dof["u"]).sub(1),  0.0,    geometry.boundaries, settings.boundary["walls"], BCmethod),
            # outlet
            # Stress free outflow imposed via Neumann boundary.
            ]        
        self.set_backflow_boundaries(geometry.ds(settings.boundary["outlet"]))

# Compute the base flow.
parameters["refinement_algorithm"] = "plaza_with_parent_facets"
parameters["reorder_dofs_serial"]  = False
parameters["allow_extrapolation"]  = True

settings = Settings(mesh_folder        = "mesh_ref",
                    mesh_base_name     = "backward_facing_step", 
                    mesh_name          = "backward_facing_step",
                    mesh_out_name      = "backward_facing_step",
                    result_folder      = str(Re), 
                    result_name        = "baseflow",
                    restart_name       = restart_name, #baseflow or no_restart
                    boundary           = {"inlet":1, "walls":2, "outlet":3}, 
                    coefficients       = {"Re":Re, "Pe":0.7*Re, "gamma":1.4, "Ma":0.01}
)

if settings.restart_name == "baseflow":
    if not Re == 100:
        settings.restart_name = str(Re - 100) + "/baseflow"
    else:
        settings.restart_name = str(Re)

# Write set of coefficients to result folder.
np.savez(settings.result_folder + '/coefficients.npz',
         Re    = settings.coefficients["Re"],
         Pe    = settings.coefficients["Pe"],
         gamma = settings.coefficients["gamma"],
         Ma    = settings.coefficients["Ma"])

# Create discretization on mesh.
geometry = Geometry(settings=settings)

# Create equation object and set boundary conditions.
equation = BackwardFacingStepNSE(settings=settings, geometry=geometry)
equation.set_boundary_conditions_nonlinear_opreator(geometry=geometry, settings=settings)

# Set initial conditions if we do not restart from another steady state.
if(settings.restart_name == "no_restart"):
    assigner = FunctionAssigner(equation.q0.function_space(), equation.spaces)
    assigner.assign(equation.q0, [project(Expression(("1.0", "0.0"), degree = 2), geometry.VP2),
                                project(Expression("1.0",          degree = 1), geometry.UP1),
                                project(Expression("1.0",          degree = 1), geometry.UP2)])

# Load inital conditions if we restart.
else:
    try:
        if settings.restart_name == "baseflow":
            initial_condition_np = np.load(settings.result_folder + '/' + settings.restart_name + '.npz')
        else:
            initial_condition_np = np.load(settings.restart_name + '.npz')
    except:
        print("ERROR Could not load file: " + settings.result_folder + '/' + settings.restart_name + '.npz' + ". Does it exist? Try restart_name = no_restart.")
        exit()
    
    equation.q0_list[equation.dof["u"]].sub(0).vector().set_local(initial_condition_np["ux"])
    equation.q0_list[equation.dof["u"]].sub(1).vector().set_local(initial_condition_np["uy"])
    equation.q0_list[equation.dof["T"]].vector().set_local(initial_condition_np["T"])
    equation.q0_list[equation.dof["p"]].vector().set_local(initial_condition_np["p"])
    
    assigner = FunctionAssigner(equation.q0.function_space(), equation.spaces)
    assigner.assign(equation.q0, equation.q0_list)

# Solve for the steady state with a Newton solver.
fixpoint_solver = Fixpoint_solver()
fixpoint_solver.solve_fixpoint(equation=equation, geometry=geometry)

# Wirte solution to paraview.
assigner = FunctionAssigner(equation.spaces, equation.q0.function_space())
assigner.assign(equation.q0_list, equation.q0)

fields_to_write = {}
fields_to_write['u']  = equation.q0_list[equation.dof["u"]]
fields_to_write['p'] = equation.q0_list[equation.dof["p"]]
fields_to_write['T']  = equation.q0_list[equation.dof["T"]]

io_interface = Io()
io_interface.write_paraview(geometry = geometry,
                        settings     = settings, 
                        file_name    = settings.result_name,
                        fields       = fields_to_write)

# Write solution to numpy compressed file for loading it later into linear solver / restart.:
np.savez(settings.result_folder + '/' + str(settings.result_name),
    ux   = equation.q0_list[equation.dof["u"]].sub(0).vector().get_local(),
    uy   = equation.q0_list[equation.dof["u"]].sub(1).vector().get_local(),
    p   = equation.q0_list[equation.dof["p"]].vector().get_local(),
    T   = equation.q0_list[equation.dof["T"]].vector().get_local()) 

# Assemble and export matrices.
settings = Settings(mesh_folder     = "mesh_ref",
                    mesh_name       = "backward_facing_step",
                    result_folder   = str(Re),
                    baseflow_folder = str(Re),
                    baseflow_name   = "baseflow",
                    result_name     = "result",
                    boundary        = {"inlet":1, "walls":2, "outlet":3}
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
equation  = BackwardFacingStepNSE(settings=settings, geometry=geometry)

# Load base flow.
baseflow_np = np.load(settings.baseflow_folder + '/' + settings.baseflow_name + '.npz')

equation.q0_list[equation.dof["u"]].sub(0).vector().set_local(baseflow_np["ux"])
equation.q0_list[equation.dof["u"]].sub(1).vector().set_local(baseflow_np["uy"])
equation.q0_list[equation.dof["p"]].vector().set_local(baseflow_np["p"])
equation.q0_list[equation.dof["T"]].vector().set_local(baseflow_np["T"])

assigner = FunctionAssigner(equation.q0.function_space(), equation.spaces)
assigner.assign(equation.q0, equation.q0_list)

# Save used baseflow (for debuggin only).
io = Io()
io.write_paraview(geometry, settings, "baseflow",
                    {"u": equation.q0_list[equation.dof["u"]],
                     "p":equation.q0_list[equation.dof["p"]],
                     "T": equation.q0_list[equation.dof["T"]]})

# Important: This must be done after the baseflow is loaded, but before we build the weakforms.
equation.set_boundary_conditions_linear_operator_homogeneous(geometry, settings)

# Assemble mass and stiffness matrices.
equation.build_jacobian(geometry)
equation.build_A_matrix()

L = PETScMatrix()
A = PETScMatrix()
M = PETScMatrix()

L_mass = PETScMatrix()
L_mom  = PETScMatrix()
L_ene  = PETScMatrix()
L_spec = PETScMatrix()

print("Assemble L matrix (Jacobian).")
assemble(equation.L_weakform, tensor = L)

print("Assemble A matrix (Mass).")
assemble(equation.A_weakform, tensor = A)

[bc.apply(L) for bc in equation.bcs_linear_operator_homogeneous]
[bc.apply(A) for bc in equation.bcs_linear_operator_homogeneous]

L_pet = as_backend_type(L).mat()
A_pet = as_backend_type(A).mat()

# Assemble the forcing norm: L2 norm for axial- and normal body forcing.
print("Assemble forcing norm Qf.")
forcing_dofs = equation.VMixed.sub(equation.dof["u"]).dofmap().dofs()

Qforcform = (equation.v_mom[0]*equation.u[0] + equation.v_mom[1]*equation.u[1]) * dx
Bform     = (equation.v_mom[0]*equation.u[0] + equation.v_mom[1]*equation.u[1]) * dx

Qforc_full = assemble(Qforcform)
B_full     = assemble(Bform)

m = Qforc_full.size(0)
n = np.size(forcing_dofs)

# Build map from forcing subspace into full system state space.
P_csr = scipy.sparse.csr_matrix((np.ones(n), (forcing_dofs, np.arange(n))), (m,n))
P_pet = PETSc.Mat().createAIJ(size=P_csr.shape, csr=(P_csr.indptr, P_csr.indices, P_csr.data))

Qforc_full_csr = convert_2_csrmat(Qforc_full)
Qforc_csr      = P_csr.transpose()*Qforc_full_csr*P_csr
Qforc_pet      = PETSc.Mat().createAIJ(size=Qforc_csr.shape, csr=(Qforc_csr.indptr, Qforc_csr.indices, Qforc_csr.data))

# Assemble the extensor matrix.
print("Assemble extensor B.")
B_full_csr = convert_2_csrmat(B_full)
B_csr      = B_full_csr*P_csr
B_pet      = PETSc.Mat().createAIJ(size=B_csr.shape, csr=(B_csr.indptr, B_csr.indices, B_csr.data))

# Assemble the response norm: L2 norm for velocity (kinetic energy)
print("Assemble response norm Qr.")
Qrespform = (equation.v_mom[0]*equation.u[0] + equation.v_mom[1]*equation.u[1]) * dx
Qresp     = assemble(Qrespform)
Qresp_csr = convert_2_csrmat(Qresp)
Qresp_pet = PETSc.Mat().createAIJ(size=Qresp_csr.shape, csr=(Qresp_csr.indptr, Qresp_csr.indices, Qresp_csr.data))

# Save in scipy sparse format.
ai, aj, av = A_pet.getValuesCSR()
A = scipy.sparse.csr_matrix((av, aj, ai))
scipy.sparse.save_npz(settings.result_folder + '/A', A)

li, lj, lv = L_pet.getValuesCSR()
L = scipy.sparse.csr_matrix((lv, lj, li))
scipy.sparse.save_npz(settings.result_folder + '/L', L)

qri, qrj, qrv = Qresp_pet.getValuesCSR()
Qresp = scipy.sparse.csr_matrix((qrv, qrj, qri))
scipy.sparse.save_npz(settings.result_folder + '/Qresp', Qresp)
 
qfi, qfj, qfv = Qforc_pet.getValuesCSR()
Qforc = scipy.sparse.csr_matrix((qfv, qfj, qfi))
scipy.sparse.save_npz(settings.result_folder + '/Qforc', Qforc)
 
bi, bj, bv = B_pet.getValuesCSR()
B = scipy.sparse.csr_matrix((bv, bj, bi))
scipy.sparse.save_npz(settings.result_folder + '/B', B)
 
pi, pj, pv = P_pet.getValuesCSR()
P = scipy.sparse.csr_matrix((pv, pj, pi))
scipy.sparse.save_npz(settings.result_folder + '/P', P)

print("Done!")
