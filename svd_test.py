#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Import packages.
from petsc4py import PETSc
from slepc4py import SLEPc

import numpy as np
import scipy

from matplotlib import pyplot as plt
import csv
import sys

# Check anaconda environment. Script needs complex petsc/slepc. Scalar type should be complex128.
from petsc4py import PETSc
print("Scalar type: " + str(PETSc.ScalarType))


# In[2]:


Re = 200
base_dir = './backward_facing_step_compressible/'
case_path = base_dir + str(Re) + '/'


# In[3]:


# Load matrices.
L = scipy.sparse.load_npz(case_path + 'L.npz')
A = scipy.sparse.load_npz(case_path + 'A.npz')

# Important: Here equations are formulatet as (\sigma A + L)q = 0.
# However, slepc expects Lq = \sigma A q.
# => Multiply A matrix with -1, before stability analysis.
A = A * (-1)


# In[4]:


# # Visualize the sparsity structure using spy plots.
# plt.figure(figsize=(10, 4))
# plt.subplot(1, 2, 1)
# plt.title("Sparsity pattern of A")
# plt.spy(A, markersize=2)

# plt.subplot(1, 2, 2)
# plt.title("Sparsity pattern of L")
# plt.spy(L, markersize=2)
# plt.tight_layout()
# plt.show()


# In[5]:


def normalize_matrix(M):
    """
    In-place normalization of a PETSc AIJ matrix M column-by-column using its CSR data.
    
    For each column j, compute:
      mean[j] = (sum_i M(i,j)) / m
      std[j]  = sqrt( (sum_i |M(i,j)|^2)/m - |mean[j]|^2 )
      
    Then update each stored nonzero entry:
      M(i,j) = (M(i,j) - mean[j]) / std[j]
      
    The matrix M is modified in-place.
    """
    # Get global sizes as integers.
    m_global, n_global = M.getSize()
    m_global = int(m_global)
    n_global = int(n_global)

    ia, ja, a = M.getValuesCSR()
    ja = ja.ravel()  # ensure column indices are 1D

    # Allocate per-column accumulators.
    col_sum   = np.zeros(n_global, dtype=a.dtype)
    col_sumsq = np.zeros(n_global, dtype=np.float64)

    # Accumulate sums from the stored nonzero entries.
    np.add.at(col_sum, ja, a)
    np.add.at(col_sumsq, ja, np.abs(a)**2)
    
    # Compute per-column mean and standard deviation.
    mean = col_sum / m_global
    std  = np.sqrt(col_sumsq / m_global - np.abs(mean)**2)
    
    # Update the nonzero entries in-place.
    a[:] = (a - mean[ja]) / std[ja]
    
    # Commit the new values back to M.
    M.setValuesCSR(ia, ja, a)
    M.assemble()
    return mean, std


# In[6]:


# Convert from scipy to petsc matrices.
L_pet = PETSc.Mat().createAIJ(size=L.shape, csr=(L.indptr, L.indices, L.data))
A_pet = PETSc.Mat().createAIJ(size=A.shape, csr=(A.indptr, A.indices, A.data))


# In[7]:


A_H_pet = A_pet.copy()
A_H_pet.hermitianTranspose()

L_H_pet = L_pet.copy()
L_H_pet.hermitianTranspose()


# In[8]:


# Normalize matrices A and L

meanA, stdA = normalize_matrix(A_pet)
meanL, stdL = normalize_matrix(L_pet)

meanA_H, stdA_H = normalize_matrix(A_H_pet)
meanL_H, stdL_H = normalize_matrix(L_H_pet)


# In[9]:


n_ev = 5

# Setup eigenvalue problem.
LEP = SLEPc.EPS().create()
LEP.setOperators(L_pet*L_H_pet, A_pet*A_H_pet)
LEP.setProblemType(SLEPc.EPS.ProblemType.GNHEP)
LEP.setType(SLEPc.EPS.Type.KRYLOVSCHUR)
LEP.setDimensions(n_ev, PETSc.DECIDE, PETSc.DECIDE)
# LEP.setTwoSided(True)
LEP.setWhichEigenpairs(SLEPc.EPS.Which.LARGEST_REAL)
# LEP.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_MAGNITUDE)
LEP.setTolerances(1e-9, 200)

st = LEP.getST()
ksp = st.getKSP()

ksp.setType("preonly")
pc = ksp.getPC()
pc.setType("jacobi")

LEP.setFromOptions()
LEP.view()

# Solve eigenvalue problem.
LEP.solve()

# Access solution.
tol, maxit = LEP.getTolerances()
nconv      = LEP.getConverged()
    
print("Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit))
print("Number of converged eigenpairs %d" % nconv)


# In[10]:


spectrum = np.zeros(nconv, dtype=complex)

# Create vectors of dimensions of matrix.
evR = L_pet.getVecs()[0].duplicate()
evR_array = []

print()
print("        k          ||Ax-kx||/||kx|| ")
print("----------------- ------------------")

for iEv in range(nconv):

    k     = LEP.getEigenvalue(iEv)
    error = LEP.getErrorEstimate(iEv)
        
    LEP.getEigenvector(iEv, evR)

    spectrum[iEv] = k

    if k.imag != 0.0:
        print(" %9g%+9g j %12g" % (k.real, k.imag, error))
    else:
        print(" %12f      %12g" % (k.real, error))

    evR_array.append(evR.copy())


# In[11]:


svdU = L_pet.getVecs()[0].duplicate()

svdU_array = []

for iEv in range(nconv):
    
    A_H_pet.mult(evR_array[iEv], svdU)
    # normalize eigenvector
    norm_val = svdU.norm(PETSc.NormType.NORM_2)
    svdU.scale(1.0 / norm_val)
    svdU_array.append(svdU.copy())


# In[12]:


def gram_schmidt(vec_list):
    """Perform modified Gram-Schmidt on a list of PETSc Vec objects."""
    orthonormal_vecs = []
    for v in vec_list:
        # Create a new copy to avoid modifying the original vector.
        v_copy = v.copy()
        for u in orthonormal_vecs:
            # Compute the projection coefficient alpha.
            alpha = v_copy.dot(u)
            # Remove the component in the direction u.
            v_copy.axpy(-alpha, u)
        # Normalize the vector.
        norm = v_copy.norm(PETSc.NormType.NORM_2)
        v_copy.scale(1.0 / norm)
        orthonormal_vecs.append(v_copy)
    return orthonormal_vecs

# svdU_array = gram_schmidt(svdU_array)


# In[13]:


#check for orthogonality
max_val = 0
for i in range(len(svdU_array)):
    for j in range(i+1, len(svdU_array)):
        val = svdU_array[i].dot(svdU_array[j])
        print("non-orthogonality: ", np.abs(val))
        max_val = max(max_val, np.abs(val.real))
print("Max non-orthogonality: ", max_val)


# In[14]:


#diagonal singular value matrix sigma
sigma = np.sqrt(spectrum)
np.savez(case_path + 'sigmas.npz', sigma=sigma)


# In[57]:


# Create a KSP solver for L_pet.
ksp = PETSc.KSP().create()
ksp.setOperators(L_pet)
ksp.setType("preonly")
pc = ksp.getPC()
pc.setType("lu")
ksp.setFromOptions()

# Create PETSc vectors for solution and right-hand-side.
svdV = L_pet.getVecs()[0].duplicate()
rhs = A_pet.getVecs()[0].duplicate()

svdV_array = []

# Solve for the i-th vector v[i].
for i in range(len(svdU_array)):

    # Compute rhs = sigma[i] * (A * svdU_array[i])
    A_pet.mult(svdU_array[i], rhs)
    rhs.scale(sigma[i])  # sigma stored as diagonal matrix

    # Solve L * v = rhs
    ksp.solve(rhs, svdV)
    #check for divergence
    if np.isnan(1 / svdV.getArray()[0]):
        break

    svdV_array.append(svdV.copy())


# In[58]:


#check for orthogonality
max_val = 0
for i in range(len(svdV_array)):
    for j in range(i+1, len(svdV_array)):
        val = svdV_array[i].dot(svdV_array[j])
        print("non-orthogonality: ", np.abs(val))
        max_val = max(max_val, np.abs(val.real))
print("Max non-orthogonality: ", max_val)


# # export SVD

# In[59]:


import os
os.makedirs(case_path + 'svd/', exist_ok=True)


# In[60]:


for iEv in range(nconv):
    np.savez(case_path + 'svd/' + str(iEv), svdU=svdU_array[iEv].getArray(), svdV=svdV_array[iEv].getArray())


# In[18]:


# for i in range(len(EWs)):
    # sv_data = np.load(case_path + '/svd/' + str(i) + '.npz')
    # svdU_np = sv_data['svdU']
    # svdV_np = sv_data['svdV']

    # equation.q_real_list[equation.dof["u"]].vector().set_local(np.real(svdV_np[equation.VMixed.sub(equation.dof["u"]).dofmap().dofs()]))
    # equation.q_imag_list[equation.dof["u"]].vector().set_local(np.imag(svdV_np[equation.VMixed.sub(equation.dof["u"]).dofmap().dofs()]))
    # equation.f_real_list[equation.dof["u"]].vector().set_local(np.real(svdU_np[equation.VMixed.sub(equation.dof["u"]).dofmap().dofs()]))
    # equation.f_imag_list[equation.dof["u"]].vector().set_local(np.imag(svdU_np[equation.VMixed.sub(equation.dof["u"]).dofmap().dofs()]))

    # fields_to_write = {}
    # fields_to_write["svdU_real"] = equation.f_real_list[equation.dof["u"]]
    # fields_to_write["svdU_imag"] = equation.f_imag_list[equation.dof["u"]]
    # fields_to_write["svdV_real"] = equation.q_real_list[equation.dof["u"]]
    # fields_to_write["svdV_imag"] = equation.q_imag_list[equation.dof["u"]]

    # # Write eigenmodes and singular vectors.
    # io = Io()
    # io.write_paraview(geometry, settings, "SVD_" + str(i), fields_to_write)


# ## reconstruct M from SVD

# In[65]:


# M = U * Sigma * V^T
# export Mk for each k
# Mk = U_k * Sigma_k * V_k^T


# Create a PETSc shell matrix representing M = Σₖ σₖ (svdUₖ ⊗ svdVₖ)

# Define a Python context class for the shell matrix.
class SVDMat(object):
    def __init__(self, nconv, svdU_array, svdV_array, sigma):
        self.nconv = nconv
        self.svdU_array = svdU_array
        self.svdV_array = svdV_array
        self.sigma = sigma

    def mult(self, mat, x, y):
        # y = Σₖ sigmaₖ * (svdVₖ^T * x) * svdUₖ
        y.zeroEntries()
        for i in range(self.nconv):
            # Compute dot product v_k^T * x.
            dot = self.svdV_array[i].dot(x)
            coeff = self.sigma[i] * dot
            y.axpy(coeff, self.svdU_array[i])
        y.assemble()

# Determine global sizes from one of the PETSc Vec objects.
m = svdU_array[0].getSize()
n = svdV_array[0].getSize()

# Create the shell matrix using PETSc's Python matrix type.
M_shell = PETSc.Mat().create(PETSc.COMM_WORLD)
M_shell.setSizes([[m, None], [n, None]])
M_shell.setType('python')
# Set the Python context that provides the matrix-vector multiply.
M_shell.setPythonContext(SVDMat(nconv, svdU_array, svdV_array, sigma))
M_shell.setUp()

# Test the shell matrix multiplication.
x = PETSc.Vec().createSeq(n, comm=PETSc.COMM_SELF)
x.setRandom()
y = PETSc.Vec().createSeq(m, comm=PETSc.COMM_SELF)
M_shell.mult(x, y)
print("Shell matrix multiplication completed.")


# In[66]:


y.getArray()


# In[89]:


baseflow_ux = np.load(case_path + 'baseflow.npz')['ux']
baseflow_uy = np.load(case_path + 'baseflow.npz')['uy']
P_mat = np.load(case_path + 'P.npz')
P_csr = scipy.sparse.csr_matrix((P_mat['data'], P_mat['indices'], P_mat['indptr']), shape=P_mat['shape'])


# In[94]:


baseflow_ux


# In[90]:


baseflow = P_csr*baseflow_ux + P_csr*baseflow_uy


# In[91]:


baseflow


# In[92]:


# draw sparsity pattern of baseflow #should be only in velocity dofs


# In[2]:


get_ipython().system('jupyter nbconvert --to script svd_test.ipynb')


# In[ ]:




