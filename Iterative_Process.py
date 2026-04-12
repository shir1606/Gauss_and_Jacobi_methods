import numpy as np

def jacobi(A, b, X):
    n = len(X)
    X_new = np.zeros(n)
    for i in range(n):
        sum_ax = sum(A[i][j] * X[j] for j in range(n) if j != i)
        X_new[i] = (b[i] - sum_ax) / A[i][i]
    return X_new

def gauss_seidel(A, b, X):
    n = len(X)
    X_new = X.copy()
    for i in range(n):
        sum_ax = sum(A[i][j] * X_new[j] for j in range(i)) + sum(A[i][j] * X[j] for j in range(i+1, n))
        X_new[i] = (b[i] - sum_ax) / A[i][i]
    return X_new

def iterative_solver(method, A, b, X0, eps=0.001, max_iter=1000):
    XR = X0.copy()
    for r in range(max_iter):
        print(f"Iteration {r}: X = {XR}")
        if method == 'jacobi':
            XRnew = jacobi(A, b, XR)
        elif method == 'gauss_seidel':
            XRnew = gauss_seidel(A, b, XR)
        else:
            raise ValueError("Method must be 'jacobi' or 'gauss_seidel'")
        # Check convergence
        if np.linalg.norm(XRnew - XR) < eps:
            print(f"Converged at iteration {r+1}: X = {XRnew}")
            return XRnew
        XR = XRnew
    print(f"Did not converge after {max_iter} iterations. Final X = {XR}")
    return XR
