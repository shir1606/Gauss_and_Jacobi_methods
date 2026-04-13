from xml.etree.ElementPath import xpath_tokenizer_re

import numpy as np
TOL = 1e-5
MAX_ITER = 1000


MaT =[[60,62,77,76,26],[11,51,1,59,95],[20,44,12,34,97],[3,65,4,13,48],[3,31,83,5,36]]
Size=5
b=[1,2,3,4,5]
#Main
def task():
    pivoting(MaT,Size)#calling for pivoting the matrix

    #calling for iteration
    method='jacobi'
    X0=[1,2,3,4,5]
    iterative_solver(MaT, b, X0,method )
    method='gauss_seidel'
    iterative_solver( MaT, b, X0, method)


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

# =====================================

def pivoting(mat,size):
    #the function ordering by max item under the diagonal the matrix
    for col in range(size):#for each level in diagonal checks the max row and ordering it
        max=findmaxrow(mat,size,col)
        if (col!= max):
            changeroworder(mat, size, col, max)

def findmaxrow(mat,size,col):
    max=mat[col][col]#starts with the item in diagonal
    maxrow=col#remember the first row of iteration
    for row in range(col+1,size):#for each item under diagonal
        if abs(mat[row][col])>abs(max):#if absolute item under diagonal is bigger then max then update
            max = mat[row][col]
            maxrow = row
    return maxrow

def changeroworder(mat, size, row, maxrow):
    #swaps the items in maxrow to the row currently working with

    for i in range(row , size):
        temp = mat[i][row]
        mat[i][row] = mat[maxrow][row]
        mat[maxrow][row] = temp

# =====================================

def diagonally_dominant(matrix):
    #the function checks if the diagonal value is bigger then the other objects on its row
    n = len(matrix)
    for i in range(n):
        diagonal_val = abs(matrix[i][i])
        off_diagonal_sum = sum(abs(matrix[i][j]) for j in range(n) if i != j)

        if diagonal_val <= off_diagonal_sum:
            return False
    return True

# =====================================

def iterative_solver(A, b, X0,  method='jacobi'):
    """
    Main algorithm loop to find the solution to Ax = b.
    """
    max_iter = 1000
    epsilon = 1e-6
    # Ensure inputs are float arrays to prevent integer division issues
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    Xr = np.array(X0, dtype=float)

    print(f"--- Starting {method.upper()} Method ---")

    for iteration in range(max_iter):
        # Print the current vector Xr
        print(f"Iteration {iteration:3}: Xr = {Xr}")

        # Compute the next approximation vector (Xr+1)
        if method == 'jacobi':
            Xr_next = jacobi(A, b)
        elif method == 'gauss_seidel':
            Xr_next = gauss_seidel(A, b)
        else:
            raise ValueError("Method must be 'jacobi' or 'gauss_seidel'")

        # Check stopping condition: |Xr+1 - Xr| < ε
        # Using the infinity norm (maximum absolute difference)
        diff = np.max(np.abs(Xr_next - Xr))

        if diff < epsilon:
            print(f"\nConvergence achieved after {iteration + 1} iterations.")
            print(f"Final Solution: X = {Xr_next}")
            break

        # Update: Xr to Xr+1
        Xr = Xr_next

    has_dominant = diagonally_dominant(A)

    print("\n--- Final Result ---")

    if diff < epsilon:
        if has_dominant:
            print("Matrix is diagonally dominant and the system converged.")
        else:
            print("Matrix is not diagonally dominant, but the system converged.")

        print("Solution:", Xr_next)
        print("Iterations:", iteration + 1)
    else:
        print("The system did not converge after", max_iter, "iterations.")
        print("Last approximation:", Xr)

    return Xr


# =====================================
task()
