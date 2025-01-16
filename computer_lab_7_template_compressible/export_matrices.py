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
from equations.compressible.compressibleNSE import CompressibleNSE

from settings                       import Settings
from fixpoint                       import Fixpoint_solver
from geometry                       import Geometry
from utilities.io                   import Io

from dolfin import *
from settings import Settings
from petsc4py import PETSc    as PETSc
from slepc4py import SLEPc    as SLEPc
import numpy as np
import scipy

# Check anaconda environment. This script needs FEniCS loaded. Scalar type should be float64.
from petsc4py import PETSc
print("Scalar type: " + str(PETSc.ScalarType))

####################################################################
# Set boundary conditions.
####################################################################
class CylinderIncompressibleNSE(CompressibleNSE):

    def set_boundary_conditions_nonlinear_opreator(self, geometry, settings):
        BCmethod = "geometric"
        self.bcs_nonlinear_operator = [
            # inlet
            DirichletBC(self.VMixed.sub(self.dof["u"]).sub(0),  1.0,    geometry.boundaries, settings.boundary["inlet"],    BCmethod),
            DirichletBC(self.VMixed.sub(self.dof["u"]).sub(1),  0.0,    geometry.boundaries, settings.boundary["inlet"],    BCmethod),
            DirichletBC(self.VMixed.sub(self.dof["T"]),         1.0,    geometry.boundaries, settings.boundary["inlet"],    BCmethod),
            # cylinder
            DirichletBC(self.VMixed.sub(self.dof["u"]).sub(0),  0.0,    geometry.boundaries, settings.boundary["cylinder"], BCmethod),
            DirichletBC(self.VMixed.sub(self.dof["u"]).sub(1),  0.0,    geometry.boundaries, settings.boundary["cylinder"], BCmethod),
            DirichletBC(self.VMixed.sub(self.dof["T"]),         1.0,    geometry.boundaries, settings.boundary["cylinder"], BCmethod),
            # outlet
            #DirichletBC(self.VMixed.sub(self.dof["p"]),         0.0,    geometry.boundaries, settings.boundary["outlet"],   BCmethod),
            # symmetry
            DirichletBC(self.VMixed.sub(self.dof["u"]).sub(1),  0.0,    geometry.boundaries, settings.boundary["symmetry"], BCmethod),
            #DirichletBC(self.VMixed.sub(self.dof["p"]),         0.0,    geometry.boundaries, settings.boundary["symmetry"], BCmethod),
            DirichletBC(self.VMixed.sub(self.dof["T"]),         1.0,    geometry.boundaries, settings.boundary["symmetry"], BCmethod),
            
        ] 
        self.set_backflow_boundaries(geometry.ds(settings.boundary["outlet"]))

    def set_boundary_conditions_linear_operator_homogeneous(self, geometry, settings):
        BCmethod = "geometric"
        self.bcs_linear_operator_homogeneous = [
            # inlet
            DirichletBC(self.VMixed.sub(self.dof["u"]).sub(0),  0.0,    geometry.boundaries, settings.boundary["inlet"],    BCmethod),
            DirichletBC(self.VMixed.sub(self.dof["u"]).sub(1),  0.0,    geometry.boundaries, settings.boundary["inlet"],    BCmethod),
            DirichletBC(self.VMixed.sub(self.dof["T"]),         0.0,    geometry.boundaries, settings.boundary["inlet"],    BCmethod),
            # cylinder
            DirichletBC(self.VMixed.sub(self.dof["u"]).sub(0),  0.0,    geometry.boundaries, settings.boundary["cylinder"], BCmethod),
            DirichletBC(self.VMixed.sub(self.dof["u"]).sub(1),  0.0,    geometry.boundaries, settings.boundary["cylinder"], BCmethod),
            DirichletBC(self.VMixed.sub(self.dof["T"]),         0.0,    geometry.boundaries, settings.boundary["cylinder"], BCmethod),
            # outlet
            #DirichletBC(self.VMixed.sub(self.dof["p"]),         0.0,    geometry.boundaries, settings.boundary["outlet"],   BCmethod), 
            # symmetry
            DirichletBC(self.VMixed.sub(self.dof["u"]).sub(1),  0.0,    geometry.boundaries, settings.boundary["symmetry"], BCmethod),
            #DirichletBC(self.VMixed.sub(self.dof["p"]),         0.0,    geometry.boundaries, settings.boundary["symmetry"], BCmethod),
            DirichletBC(self.VMixed.sub(self.dof["T"]),         0.0,    geometry.boundaries, settings.boundary["symmetry"], BCmethod)
        ]
        self.set_backflow_boundaries(geometry.ds(settings.boundary["outlet"]))

######################################################################################
# Compute base flow
######################################################################################
parameters["refinement_algorithm"] = "plaza_with_parent_facets"
parameters["reorder_dofs_serial"]  = False
parameters["allow_extrapolation"]  = True

Re = 70

settings = Settings(mesh_folder        = "mesh_half",
                    mesh_base_name     = "cylinder_half_base", 
                    mesh_name          = "cylinder_half",
                    result_folder      = "baseflow", 
                    mesh_out_name      = "mesh_half",
                    result_name        = "baseflow",
                    restart_name       = "no_restart", #no_restart
                    boundary           = {"inlet":1, "outlet":2, "symmetry":3, "cylinder":4}, 
                    coefficients       = {"Re":Re, "Pe":0.7*Re, "gamma":1.4, "Ma":0.01}
)

# Write set of coefficients to result folder.
np.savez(settings.result_folder + '/coefficients.npz',
         Re    = settings.coefficients["Re"],
         Pe    = settings.coefficients["Pe"],
         gamma = settings.coefficients["gamma"],
         Ma    = settings.coefficients["Ma"])

# Create discretization on mesh.
geometry = Geometry(settings=settings)

# Create equation object and set boundary conditions.
equation = CylinderIncompressibleNSE(settings=settings, geometry=geometry)
equation.set_boundary_conditions_nonlinear_opreator(geometry=geometry, settings=settings)

# Set initial conditions if we do not restart from another steady state.
if(settings.restart_name == "no_restart"):
    assigner = FunctionAssigner(equation.q0.function_space(), equation.spaces)
    assigner.assign(equation.q0, [project(Expression(("1.0", "0.0"), degree = 2), geometry.VP2),
                                project(Expression("1.0",          degree = 1), geometry.UP1),
                                project(Expression("1.0",          degree = 1), geometry.UP2)])

# Load inital conditions if we restart.
else:
    inital_condition_np = np.load(settings.result_folder + '/' + settings.restart_name + '.npz')
    
    equation.q0_list[equation.dof["u"]].sub(0).vector().set_local(inital_condition_np["ux"])
    equation.q0_list[equation.dof["u"]].sub(0).vector().set_local(inital_condition_np["uy"])
    equation.q0_list[equation.dof["T"]].vector().set_local(inital_condition_np["T"])
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
                        file_name    = settings.result_name + "_half",
                        fields       = fields_to_write)

# Write solution to numpy compressed file for loading it later into linear solver / restart.:
np.savez(settings.result_folder + '/' + str(settings.result_name) + '_half.npz',
    ux  = equation.q0_list[equation.dof["u"]].sub(0).vector().get_local(),
    uy  = equation.q0_list[equation.dof["u"]].sub(1).vector().get_local(),
    p   = equation.q0_list[equation.dof["p"]].vector().get_local(),
    T   = equation.q0_list[equation.dof["T"]].vector().get_local()) 

# Mirror solution on full mesh.
print("Mirror solution on full mesh.")

mesh_full = Mesh("mesh_full/cylinder_full.xml")
domains   = MeshFunction('size_t', mesh_full, mesh_full.topology().dim())
P1_full   = FiniteElement('P', triangle, 1)
P2_full   = FiniteElement('P', triangle, 2)
UP1_full  = FunctionSpace(mesh_full, P1_full)
UP2_full  = FunctionSpace(mesh_full, P2_full)
VP2_full  = VectorFunctionSpace(mesh_full, P2_full, 2)

class MirrorScalar(UserExpression):
    def __init__(self, expression_to_mirror):
        super().__init__()
        self.expression_to_mirror = expression_to_mirror
    def eval(self, values, x):        
        if(x[1] >= 0):
            values[0] = self.expression_to_mirror(Point( x[0], x[1]))
        else:
            values[0] = self.expression_to_mirror(Point( x[0],-x[1]))
    def value_shape(self):
        return()

class MirrorVector(UserExpression):
    def __init__(self, expression_to_mirror):
        super().__init__()
        self.expression_to_mirror = expression_to_mirror
    def eval(self, values, x):        
        if(x[1] >= 0):
            values[0] = self.expression_to_mirror(Point( x[0], x[1]))[0]
            values[1] = self.expression_to_mirror(Point( x[0], x[1]))[1]
        else:
            values[0] =    self.expression_to_mirror(Point( x[0], -x[1]))[0]
            values[1] = -1*self.expression_to_mirror(Point( x[0], -x[1]))[1]
    def value_shape(self):
        return(2,)

p_full = interpolate(MirrorScalar(Expression("p", p = equation.q0_list[equation.dof["p"]],       degree=1)),  UP1_full)
T_full  = interpolate(MirrorScalar(Expression("T",  T  = equation.q0_list[equation.dof["T"]],        degree=2)), UP2_full)
u_full  = interpolate(MirrorVector(Expression(("ux","uy"),  ux  = equation.q0_list[equation.dof["u"]].sub(0), uy  = equation.q0_list[equation.dof["u"]].sub(1),        degree=2)), VP2_full)

xdmf_file_full = XDMFFile(mesh_full.mpi_comm(), settings.result_folder + "/" +str(settings.result_name) + "_full.xdmf")
xdmf_file_full.parameters["rewrite_function_mesh"] = False
xdmf_file_full.parameters["functions_share_mesh"]  = True

fields_to_write_full  = {}
fields_to_write_full["u"]   = u_full
fields_to_write_full["p"]  = p_full
fields_to_write_full["T"]   = T_full

for field_name in fields_to_write_full:
    fields_to_write_full[field_name].rename(str(field_name), str(field_name))
    xdmf_file_full.write(fields_to_write_full[field_name],0)
xdmf_file_full.close()

np.savez(settings.result_folder + '/' + str(settings.result_name) + '_full.npz',
    u  = u_full.vector().get_local(),
    p  = p_full.vector().get_local(),
    T  = T_full.vector().get_local()) 

###################################################################################
# Assemble and export matrices.
###################################################################################
settings = Settings(mesh_folder     = "mesh_full",
                    mesh_name       = "cylinder_full",
                    result_folder   = "result",
                    baseflow_folder = "baseflow",
                    baseflow_name   = "baseflow_full",
                    result_name     = "result",
                    boundary        = {"inlet":1, "outlet":2, "cylinder":3, "symmetry":5}
                   )

# Load coefficients of base flow.
coefficients = np.load(settings.baseflow_folder + '/coefficients.npz')
settings.coefficients = {"Re":    coefficients["Re"],
                         "Pe":    coefficients["Pe"],
                         "gamma": coefficients["gamma"],
                         "Ma":    coefficients["Ma"]}

from geometry                       import Geometry
from utilities.io                   import Io
from dolfin import *

# Make sure this flag is set if the base flow was computed with it!!!
parameters["reorder_dofs_serial"] = False

# Create discretization on mesh.
geometry = Geometry(settings)

# Setup equation.
equation  = CylinderIncompressibleNSE(settings=settings, geometry=geometry)

# Load base flow.
baseflow_np = np.load(settings.baseflow_folder + '/' + settings.baseflow_name + '.npz')

equation.q0_list[equation.dof["u"]].vector().set_local(baseflow_np["u"])
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

# Assemble.
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

# Save in scipy sparse format.
ai, aj, av = A_pet.getValuesCSR()
A = scipy.sparse.csr_matrix((av, aj, ai))
scipy.sparse.save_npz(settings.result_folder + '/A', A)

li, lj, lv = L_pet.getValuesCSR()
L = scipy.sparse.csr_matrix((lv, lj, li))
scipy.sparse.save_npz(settings.result_folder + '/L', L)