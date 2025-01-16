import numpy as np
import matplotlib.pyplot as plt

results_dir = './backward_facing_step_compressible/result/'

# Qforc_full_csr = convert_2_csrmat(Qforc_full)
# Qforc_csr      = P_csr.transpose()*Qforc_full_csr*P_csr
# Qforc_pet      = PETSc.Mat().createAIJ(size=Qforc_csr.shape, csr=(Qforc_csr.indptr, Qforc_csr.indices, Qforc_csr.data))

# # Assemble the response norm: L2 norm for velocity (kinetic energy)
# print("Assemble response norm Qr.")
# Qrespform = (equation.v_mom[0]*equation.u[0] + equation.v_mom[1]*equation.u[1]) * dx
# Qresp     = assemble(Qrespform)
# Qresp_csr = convert_2_csrmat(Qresp)
# Qresp_pet = PETSc.Mat().createAIJ(size=Qresp_csr.shape, csr=(Qresp_csr.indptr, Qresp_csr.indices, Qresp_csr.data))

# qfi, qfj, qfv = Qforc_pet.getValuesCSR()
# Qforc = scipy.sparse.csr_matrix((qfv, qfj, qfi))
# scipy.sparse.save_npz(settings.result_folder + '/Qforc', Qforc)


# sparse matrices
Qforc = np.load(results_dir + 'Qforc.npz')
Qresp = np.load(results_dir + 'Qresp.npz')

print(Qforc.files)
#['indices', 'indptr', 'format', 'shape', 'data']
print(Qforc['indices'][:20])
print(Qforc['indptr'][:20])