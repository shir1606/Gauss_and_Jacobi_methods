# Iterative Methods for Solving Linear Systems (Jacobi & Gauss-Seidel)

## Project Description

This project implements numerical methods for solving systems of linear equations of the form:

**A·x = b**

where:
- `A` is a square matrix (n x n)
- `b` is a column vector (right-hand side)
- `x` is the solution vector

The project includes implementation of two iterative methods:
1. Jacobi Method
2. Gauss-Seidel Method

In addition, the project includes preprocessing steps such as pivoting and checking diagonal dominance to ensure convergence.



## Implemented Features

### 1. System Solver Function
A main function that receives:
a. A square matrix `A`
b. A solution vector `b`

and solves the system using iterative methods.


### 2. Pivoting (Row Rearrangement)
Before applying the iterative methods, the matrix is processed using partial pivoting.

For each column `i`:
- The algorithm searches for the maximum absolute value in column `i`, starting from row `i` downward.
- Rows are swapped to place the maximum element on the diagonal position.

Important rule:
Only rows below (or equal to) the current index are considered for swapping.


### 3. Diagonal Dominance Check
After pivoting, the program verifies whether the matrix is diagonally dominant.

For each row `i`:

\[
|A[i][i]| \ge \sum_{j \ne i} |A[i][j]|
\]

If this condition holds for all rows, the matrix is considered diagonally dominant, which increases the likelihood of convergence.

---

### 4. Initialization of Solution Vectors

- `Xr`: initial approximation vector, initialized to zeros
- `Xr+1`: updated solution vector computed in each iteration

Stopping criterion:
- ε (epsilon) = 0.001


### 5. Iterative Process

The algorithm runs up to **1000 iterations** or until convergence.

In each iteration:
1. Compute the next approximation vector (`Xr+1`) using either:
   - Jacobi method OR
   - Gauss-Seidel method (implemented in separate functions)
2. Print the current vector `Xr`
3. Update:
   - `Xr to Xr+1`
4. Check stopping condition:
   - If `|Xr+1 - Xr| < ε`, stop
   - Otherwise continue


### 6. Stopping Condition & Output

The program stops when:
- The solution converges (difference < epsilon), OR
- The number of iterations reaches 1000

Final output includes:
- Final solution vector
- Whether the matrix is diagonally dominant
- Whether the method converged successfully

  
## Project Structure
├── main.cpp                # Main program flow
├── jacobi.cpp / .h         # Jacobi implementation
├── gauss_seidel.cpp / .h   # Gauss-Seidel implementation
├── pivoting.cpp / .h       # Pivoting logic
├── utils.cpp / .h          # Helper functions (norms, checks)
└── README.md               # Project documentation
