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

def iterative_solve(method, A, b, X0, eps=1e-6, max_iter=1000):
    X = X0.copy()
    for r in range(max_iter):
        print(f"Iteration {r}: X = {X}")
        if method == 'jacobi':
            X_next = jacobi(A, b, X)
        elif method == 'gauss_seidel':
            X_next = gauss_seidel(A, b, X)
        else:
            raise ValueError("Method must be 'jacobi' or 'gauss_seidel'")
        # Check convergence
        if np.linalg.norm(X_next - X) < eps:
            print(f"Converged at iteration {r+1}: X = {X_next}")
            return X_next
        X = X_next
    print(f"Did not converge after {max_iter} iterations. Final X = {X}")
    return X
