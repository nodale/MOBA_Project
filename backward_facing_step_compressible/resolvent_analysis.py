################################################################################################
### Computational Methods for Operator-Based Analysis
#Technical University of Munich, Professorship for Thermofluid Dynamics | Pr. Polifke, Ph.D
#
#Created: 06/2024 | P. Brokof, G. Varillon
################################################################################################

# Include packages
from petsc4py import PETSc
from slepc4py import SLEPc

import numpy as np
import scipy.sparse
from matplotlib import pyplot as plt
import sys

# Check anaconda environmnet. This script requires complex petsc/slepc. Scalar type should be complex128.
from petsc4py import PETSc
print("Scalar type: " + str(PETSc.ScalarType))

# Get case path from system arguments.
if len(sys.argv) > 1:
    case_path = str(sys.argv[1]) + "/"
    print("Using case path: " + case_path)
else:
    case_path = "100/"
    print("No case path provided. Using default case path: " + case_path)


######################################################################
# Configuration of Hermitian eigenvalue problem
######################################################################
def setupHermitianEigenvalueProblem(LEP):
    LEP.setProblemType(SLEPc.EPS.ProblemType.GHEP)
    LEP.setType(SLEPc.EPS.Type.KRYLOVSCHUR)
    LEP.setDimensions(10, PETSc.DEFAULT, PETSc.DEFAULT)
    LEP.setWhichEigenpairs(SLEPc.EPS.Which.LARGEST_REAL)
    LEP.setTolerances(1e-9, 50)
    
    LEP.setFromOptions
    LEP.view()

#####################################################################
# Matrix free representation of the left hand side.
#####################################################################
def setup_LHS(St, A_pet, L_pet, B_pet, Qresp_pet):
    
    # Build action of the resolvent.
    R_inv_pet = A_pet.copy()
    R_inv_pet.scale(2*np.pi*St*1j)
    R_inv_pet.axpy(1.0, L_pet)
    R = PETSc.KSP().create()
    R.setOperators(R_inv_pet)
    setupKSP(R)

    # Build action of the Hermitian of the resolvent.
    RH_inv_pet = R_inv_pet.copy()
    RH_inv_pet.conjugate().transpose()
    RH = PETSc.KSP().create()
    RH.setOperators(RH_inv_pet)
    setupKSP(RH)

    '''
    Computing the resolvent directly is computationally too expensive because of the involved matrix inversion. Instead, we define in 
    PETSc an operator with the same action as the LHS of the eigenvalue problem

    LHS * f = G^2 * Q_f *f
    with LHS = B^H R^H Q_q R B.

    Implement the matrix vector multiplication of the LHS in the function `mult(cls, J, x, y)` below. The function shoud return the 
    result of y = LHS * x. 

    Hint: To solve the linear system R_inv*z = w (which is equivalent to z = R*w) you can use the KSP routine `R.solve(w,z)`.
    '''
    # Build action of the LHS of the eigenvalue problem.
    w, z = R_inv_pet.getVecs()
    class LHS_class:

        # y = LHS x
        def mult(cls, J, x, y):
            # Perform multiplication step by step.
            ''' put your Python code here ... '''
            B_pet.mult(x,w)
            R.solve(w,z)
            Qresp_pet.mult(z,w)
            RH.solve(w,z)
            B_pet.multHermitian(z,y)

    
    LHS = PETSc.Mat().create()
    LHS.setSizes(Qforc_pet.size)
    LHS.setType(PETSc.Mat.Type.PYTHON)
    LHS.setPythonContext(LHS_class())
    LHS.setUp()

    return [R, RH, LHS]

def setupKSP(P):
    P.setType("preonly")
    PC = P.getPC(); PC.setType("lu")
    PC.setFactorSolverType("mumps")
    P.setFromOptions()

# Load matrices.
L     = scipy.sparse.load_npz(case_path + 'L.npz')
A     = scipy.sparse.load_npz(case_path + 'A.npz')
Qresp = scipy.sparse.load_npz(case_path + 'Qresp.npz')
Qforc = scipy.sparse.load_npz(case_path + 'Qforc.npz')
B     = scipy.sparse.load_npz(case_path + 'B.npz')
P     = scipy.sparse.load_npz(case_path + 'P.npz')
 
# Convert from scipy to petsc matrices.
L_pet     = PETSc.Mat().createAIJ(size=L.shape,     csr=(L.indptr,     L.indices,     L.data))
A_pet     = PETSc.Mat().createAIJ(size=A.shape,     csr=(A.indptr,     A.indices,     A.data))
Qresp_pet = PETSc.Mat().createAIJ(size=Qresp.shape, csr=(Qresp.indptr, Qresp.indices, Qresp.data))
Qforc_pet = PETSc.Mat().createAIJ(size=Qforc.shape, csr=(Qforc.indptr, Qforc.indices, Qforc.data))
B_pet     = PETSc.Mat().createAIJ(size=B.shape,     csr=(B.indptr,     B.indices,     B.data))
P_pet     = PETSc.Mat().createAIJ(size=P.shape,     csr=(P.indptr,     P.indices,     P.data))

n=20
St_list   = np.arange(0.02, 0.12, 0.02)
gain_list = []

for St in St_list:
    print("St: " + str(St))

    R, RH, LHS = setup_LHS(St, A_pet, L_pet, B_pet, Qresp_pet)

    # Set up the Hermitian eigenvalue problem.
    LEP = SLEPc.EPS().create()
    LEP.setOperators(LHS, Qforc_pet)
    setupHermitianEigenvalueProblem(LEP)

    # Solve.
    LEP.solve()

    # Access solution information.
    its           = LEP.getIterationNumber()
    eps_type      = LEP.getType()
    nev, ncv, mpd = LEP.getDimensions()
    tol, maxit    = LEP.getTolerances()
    nconv         = LEP.getConverged()

    print("Number of iterations of the method: %d" % its)
    print("Solution method: %s" % eps_type)
    print("Number of requested eigenvalues: %d" % nev)
    print("Dimension of subspace: %d" % ncv)
    print("Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit))
    print("Number of converged eigenpairs %d" % nconv)

    if nconv > 0:

        # Create the result vector.
        evR,_ = LHS.getVecs()
        w,_   = L_pet.getVecs()

        print()
        print("        k          ||Ax-kx||/||kx|| ")
        print("----------------- ------------------")

        for iEv in range(nconv):
    
            k     = LEP.getEigenvalue(iEv)
            error = LEP.computeError(iEv)
    
            print(" %12g      %12g" % (np.sqrt(k.real), error))
    
            # Acces dominant resolvent mode.
            if iEv == 0:
                gain_list.append(np.sqrt(k.real))
                LEP.getEigenvector(0,evR)

                # Compute optimal forcing in complete DOF space.
                P_pet.mult(evR,w)
                f = w.getArray()

                # Compute optimal response in complete DOF space.
                q,_ = Qresp_pet.getVecs()
                R.solve(w,q)
                q = q.getArray()

                # Write optimal forcing and response to hard disk.
                np.savez(case_path + 'fq_' + str("{:4f}".format(St)), f = f, q=q)

                # Save gain curve.
                np.savez(case_path + 'gain_curve', St = St_list, gain = gain_list)
                
        R.destroy()
        RH.destroy()
        LHS.destroy()


gain_list = np.array(gain_list)
St_list   = np.array(St_list)

###################################################
# Plot of gain curve
##################################################
plt.plot(St_list,gain_list, 'b-x')
plt.xlabel("St")
plt.ylabel("Optimal gain")
plt.savefig("gain_curve.svg")