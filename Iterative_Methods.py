TOL = 1e-5
MAX_ITER = 100


# =====================================


# Jacobi Method
def jacobi(A, b):
    n = len(A)

    x_old = [0.0] * n
    x_new = [0.0] * n

    for iteration in range(MAX_ITER):

        for i in range(n):
            sum_not_i = 0.0

            for j in range(n):
                if j != i:
                    sum_not_i += A[i][j] * x_old[j]

            x_new[i] = (b[i][0] - sum_not_i) / A[i][i]

        print("Iteration", iteration + 1, ":", x_new)

        # תנאי עצירה
        max_diff = 0.0
        for i in range(n):
            diff = abs(x_new[i] - x_old[i])
            if diff > max_diff:
                max_diff = diff

        if max_diff < TOL:
            print("\nConverged after", iteration + 1, "iterations")
            return x_new

        x_old = x_new.copy()

    print("\nDid not converge")
    return x_old


# =====================================


# Gauss-Seidel Method
def gauss_seidel(A, b):
    n = len(A)

    x = [0.0] * n

    for iteration in range(MAX_ITER):

        x_old = x.copy()

        for i in range(n):

            for j in range(i):
                sum_before += A[i][j] * x[j]

            for j in range(i + 1, n):
                sum_after += A[i][j] * x_old[j]

            x[i] = (b[i][0] - sum_before - sum_after) / A[i][i]

        print("Iteration", iteration + 1, ":", x)

        max_diff = 0.0
        for i in range(n):
            diff = abs(x[i] - x_old[i])
            if diff > max_diff:
                max_diff = diff

        if max_diff < TOL:
            print("\nConverged after", iteration + 1, "iterations")
            return x

    print("\nDid not converge")
    return x