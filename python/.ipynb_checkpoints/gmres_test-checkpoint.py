"""
GMRES (Restarted) for unsymmetric systems, with and without ILU preconditioning.
Runs on:
    - 2D convection-diffusion system (unsymmetric)
    - random large unsymmetric diagonally dominant matrix
Reports number of iterations for each restart and preconditioning setting.
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

def convection_diffusion_2d(n, diffusion=1.0, convection=10.0):
    N = n * n
    main = 4 * diffusion * np.ones(N)
    ew_diff = -diffusion * np.ones(N - 1)
    ns = -diffusion * np.ones(N - n)
    bias = convection / (2 * (n - 1))
    ew_conv_plus  = -bias * np.ones(N - 1)
    ew_conv_minus =  bias * np.ones(N - 1)
    ew_plus  = ew_diff + ew_conv_plus    # offset +1
    ew_minus = ew_diff + ew_conv_minus   # offset -1
    diags = [main, ew_minus, ew_plus, ns, ns]
    offs  = [0,     -1,       +1,      -n,  +n]
    A = sp.diags(diags, offs, shape=(N, N), format='csr')
    for i in range(1, n):
        A[i*n,     i*n - 1] = 0
        A[i*n - 1, i*n    ] = 0
    return A

def build_random_dd(size, density=1e-3, seed1=1, seed2=2):
    R1 = sp.rand(size, size, density=density, format='csr', random_state=seed1)
    R2 = sp.rand(size, size, density=density, format='csr', random_state=seed2)
    return R1 - R2 + sp.diags(10 * np.ones(size), format='csr')

def test_gmres(A, b, restarts, tol=1e-6):
    results = []
    try:
        ilu = spla.spilu(A.tocsc(), drop_tol=1e-4)
        M = spla.LinearOperator(A.shape, matvec=lambda x: ilu.solve(x))
        preconds = [("None", None), ("ILU", M)]
    except Exception:
        preconds = [("None", None)]

    for name, M in preconds:
        for m in restarts:
            count = {'iters': 0}
            def cb(rk):
                count['iters'] += 1
            _, info = spla.gmres(
                A, b,
                restart=m,
                tol=tol,
                M=M,
                callback=cb,
                callback_type='pr_norm'
            )
            results.append((name, m, count['iters']))
    return results

def main():
    problems = {
        "Convection-Diffusion (50×50)": convection_diffusion_2d(50, diffusion=1.0, convection=50.0),
        "Random DD Sparse (5000×5000)": build_random_dd(5000, density=1e-3)
    }
    restarts = [10, 20, 40, 80]
    tol = 1e-6

    for prob_name, A in problems.items():
        b = np.ones(A.shape[0])
        print(f"\n=== {prob_name} ===")
        results = test_gmres(A, b, restarts, tol)
        print(f"{'Precond':>8s}  {'m':>4s}  {'Iters':>6s}")
        print("-" * 24)
        for pre, m, iters in results:
            print(f"{pre:>8s}  {m:4d}  {iters:6d}")

if __name__ == "__main__":
    main()