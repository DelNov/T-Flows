import numpy as np
from scipy.io import mmread
from scipy.sparse.linalg import eigsh

# 1. Read the matrix from MatrixMarket file
A = mmread("t-flows.mtx").tocsr()

print("Matrix shape:", A.shape)
print("NNZ:", A.nnz)

# 2. Compute largest eigenvalue
# 'LM' = Largest Magnitude
lam_max, _ = eigsh(A, k=1, which='LM')
lam_max = lam_max[0]

# 3. Compute smallest eigenvalue
# 'SM' = Smallest Magnitude (for SPD this is the smallest positive eigenvalue)
lam_min, _ = eigsh(A, k=1, which='SM')
lam_min = lam_min[0]

# 4. Print results
print("lambda_max =", lam_max)
print("lambda_min =", lam_min)
print("kappa      =", lam_max / lam_min)

