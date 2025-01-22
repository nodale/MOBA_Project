################################################################################################
### Computational Methods for Operator-Based Analysis
#Technical University of Munich, Professorship for Thermofluid Dynamics | Pr. Polifke, Ph.D
#
#Created: 06/2024 | P. Brokof, G. Varillon
################################################################################################

# Import packages.
from petsc4py import PETSc
from slepc4py import SLEPc

import numpy as np
import scipy.sparse

from matplotlib import pyplot as plt
import csv
import sys

# Check anaconda environment. Script needs complex petsc/slepc. Scalar type should be complex128.
from petsc4py import PETSc
print("Scalar type: " + str(PETSc.ScalarType))

# Get re from system arguments.
if(len(sys.argv) > 1):
    case_path = str(sys.argv[1]) + "/"
    print("Using case path: " + case_path)
else:
    case_path = "40/"
    print("No case path provided. Using default case path: " + case_path)

#################################################################################
# Setup Generalized Non Hermitian Eigenvalue Problem
#################################################################################
def setupGeneralizedNonHermitianEPS(LEP, shift):
    LEP.setProblemType(SLEPc.EPS.ProblemType.GNHEP)
    LEP.setType(SLEPc.EPS.Type.KRYLOVSCHUR)
    LEP.setDimensions(10, PETSc.DEFAULT, PETSc.DEFAULT)
    LEP.setTarget(shift)
    LEP.setTwoSided(True)
    #LEP.setWhichEigenpairs(SLEPc.EPS.Which.LARGEST_REAL)
    LEP.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_MAGNITUDE)
    LEP.setTolerances(1e-9, 50)

    # Setting up linear problem solver for shift invert transformation.
    #We solve direct with MUMPS.
    st = LEP.getST()
    st.setType("sinvert")
    st.setShift(shift)
    ksp = st.getKSP()
    ksp.setType("preonly")
    pc = ksp.getPC()
    pc.setType("lu")
    pc.setFactorSolverType("mumps")

    LEP.setFromOptions
    LEP.view()

##################################################
#Solver
##################################################
# Load matrices.
L = scipy.sparse.load_npz(case_path + 'L.npz')
A = scipy.sparse.load_npz(case_path + 'A.npz')

# Important: Here equations are formulatet as (\sigma A + L)q = 0.
# However, slepc expects Lq = \sigma A q.
# => Multiply A matrix with -1, before stability analysis.
A = A * (-1)

# Convert from scipy to petsc matrices.
L_pet = PETSc.Mat().createAIJ(size=L.shape, csr=(L.indptr, L.indices, L.data))
A_pet = PETSc.Mat().createAIJ(size=A.shape, csr=(A.indptr, A.indices, A.data))

# Setup eigenvalue problem.
LEP = SLEPc.EPS().create()
LEP.setOperators(L_pet, A_pet)
setupGeneralizedNonHermitianEPS(LEP, shift = 0.45)#0.7*1j) #lambda=6+2pi*St

# Solve eigenvalue problem.
LEP.solve()

# Access solution.
tol, maxit = LEP.getTolerances()
nconv      = LEP.getConverged()
    
print("Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit))
print("Number of converged eigenpairs %d" % nconv)

spectrum = np.zeros(nconv, dtype=complex)

if nconv > 0:
    # Create vectors of dimensions of matrix.
    evR, evL = L_pet.getVecs()  
    
    print()
    print("        k          ||Ax-kx||/||kx|| ")
    print("----------------- ------------------")
    
    for iEv in range(nconv):
    
        k     = LEP.getEigenvalue(iEv)
        error = LEP.getErrorEstimate(iEv)
            
        LEP.getEigenvector(iEv, evR)
        LEP.getLeftEigenvector(iEv, evL)
    
        spectrum[iEv] = k
    
        if k.imag != 0.0:
            print(" %9f%+9f j %12g" % (k.real, k.imag, error))
        else:
            print(" %12f      %12g" % (k.real, error))
    
        evR_np = evR.getArray()
        evL_np = evL.getArray()
    
        # Save eigenvector as compressed numpy archive.
        np.savez(case_path + 'ev_' + str(iEv), evR = evR_np, evL = evL_np)
    
    # Save spectrum of case.
    np.savez(case_path + 'spectrum', spectrum = spectrum)

    # Plot spectrum.
    plt.scatter(np.imag(spectrum)/(2*3.415),np.real(spectrum))
    plt.ylabel("growth rate")
    plt.xlabel("St")
    plt.savefig(case_path + "spectrum.svg")
