"""
Conjugate Gradient for 2D Laplace system:
Solves Ax = b where A is the 2D Laplace matrix (SPD), b = Ax* for random x*.
Compares CG with and without incomplete Cholesky preconditioning.
Plots convergence (residual norm per iteration).
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

def laplace_2d(n):
    K1d = sp.diags([-1, 2, -1], [-1, 0, 1], shape=(n, n), format='csr')
    I   = sp.eye(n, format='csr')
    A = sp.kron(I, K1d) + sp.kron(K1d, I)
    return A

n = 100  # system size n^2
A = laplace_2d(n)
N = n * n

x_true = np.random.rand(N)
b = A @ x_true
x0 = np.zeros(N)

# CG without preconditioning
residuals_noprec = []
def cb_noprec(rk): residuals_noprec.append(np.linalg.norm(rk))

x_cg, info = spla.cg(A, b, x0=x0, tol=1e-8, callback=cb_noprec, maxiter=1000)

# CG with incomplete Cholesky (spilu)
try:
    ilu = spla.spilu(A.tocsc(), drop_tol=1e-4, fill_factor=10)
    Mx = lambda x: ilu.solve(x)
    M = spla.LinearOperator((N, N), Mx)
    residuals_prec = []
    def cb_prec(rk): residuals_prec.append(np.linalg.norm(rk))
    x_cg_prec, info_prec = spla.cg(A, b, x0=x0, tol=1e-8, M=M, callback=cb_prec, maxiter=1000)
except Exception as e:
    print("Incomplete Cholesky preconditioning failed:", e)
    residuals_prec = None

# Plot
plt.figure(figsize=(7,5))
plt.semilogy(residuals_noprec, label="CG (no preconditioner)")
if residuals_prec is not None:
    plt.semilogy(residuals_prec, label="CG (incomplete Cholesky)")
plt.xlabel('Iteration')
plt.ylabel('Residual Norm')
plt.title(f'CG Convergence (2D Laplacian, N={N})')
plt.legend()
plt.grid(True, which='both', ls=':')
plt.tight_layout()
plt.show()