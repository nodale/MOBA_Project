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

# Check anaconda environmnet. This script requires complex petsc/slepc. Scalar type should be complex128.
from petsc4py import PETSc
print("Scalar type: " + str(PETSc.ScalarType))

num_ev = 10
tolerance_ev = 1e-4
max_iters = 50


######################################################################
# Configuration of Hermitian eigenvalue problem
######################################################################
def setupHermitianEigenvalueProblem(LEP):
    LEP.setProblemType(SLEPc.EPS.ProblemType.GHEP)
    LEP.setType(SLEPc.EPS.Type.KRYLOVSCHUR)
    LEP.setDimensions(num_ev, PETSc.DEFAULT, PETSc.DEFAULT)
    LEP.setWhichEigenpairs(SLEPc.EPS.Which.LARGEST_REAL)
    LEP.setTolerances(tolerance_ev, max_iters)
    
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
            B_pet.mult(x, w)
            R.solve(w, z)
            Qresp_pet.mult(z, w)
            RH.solve(w, z)
            B_pet.multTranspose(z, y)
            
            
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
L     = scipy.sparse.load_npz('result/L.npz')
A     = scipy.sparse.load_npz('result/A.npz')
Qresp = scipy.sparse.load_npz('result/Qresp.npz')
Qforc = scipy.sparse.load_npz('result/Qforc.npz')
B     = scipy.sparse.load_npz('result/B.npz')
P     = scipy.sparse.load_npz('result/P.npz')
 
# Convert from scipy to petsc matrices.
L_pet     = PETSc.Mat().createAIJ(size=L.shape,     csr=(L.indptr,     L.indices,     L.data))
A_pet     = PETSc.Mat().createAIJ(size=A.shape,     csr=(A.indptr,     A.indices,     A.data))
Qresp_pet = PETSc.Mat().createAIJ(size=Qresp.shape, csr=(Qresp.indptr, Qresp.indices, Qresp.data))
Qforc_pet = PETSc.Mat().createAIJ(size=Qforc.shape, csr=(Qforc.indptr, Qforc.indices, Qforc.data))
B_pet     = PETSc.Mat().createAIJ(size=B.shape,     csr=(B.indptr,     B.indices,     B.data))
P_pet     = PETSc.Mat().createAIJ(size=P.shape,     csr=(P.indptr,     P.indices,     P.data))

St_list   = [0.1, 0.11, 0.12, 0.13, 0.14, 0.15]
gain_sum_list = []

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

    if nconv >= nev:

        # Create the result vector.
        evR,_ = LHS.getVecs()
        w,_   = L_pet.getVecs()

        print()
        print("        k          ||Ax-kx||/||kx|| ")
        print("----------------- ------------------")

        gain_list = []
        for i_Ev in range(nconv):
    
            k     = LEP.getEigenvalue(i_Ev)
            error = LEP.computeError(i_Ev)
    
            print(" %12g      %12g" % (np.sqrt(k.real), error))
    
            # Acces dominant resolvent mode.
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
            np.savez('result/fq_' + str("{:4f}".format(St)), f = f, q = q)

        gain_sum_list.append(np.sum(gain_list[:nev]))

        # Save gain curve.
        np.savez('result/gain_curve' + str(St), St = St, gain = gain_list)

    elif nconv > 0:
        print("ERROR: Not all requested eigenvalues converged! Request fewer eigenvalues or increase max_iters.")
        exit(1)
    elif nconv == 0:
        print("ERROR: No eigenvalues converged!")
        exit(1)

    np.savez('result/gain_sum_curve', St = St_list, gain = gain_sum_list)



###################################################
# Plot of gain curve
##################################################
plt.plot(St_list, gain_sum_list, 'b-x')
plt.xlabel("St")
plt.ylabel("Optimal gain")
plt.savefig("gain_curve.svg")
