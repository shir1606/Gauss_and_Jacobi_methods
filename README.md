# Solving Linear Systems Using Iterative Methods

### Jacobi & Gauss-Seidel | Python Implementation

---

## 👩‍💻 Project Overview

This project presents a full implementation of two classical **iterative numerical methods** for solving systems of linear equations:

* **Jacobi Method**
* **Gauss-Seidel Method**

The implementation follows the assignment requirements and emphasizes:

* Correct mathematical logic
* Code modularity and reusability
* Handling of edge cases (non-diagonally dominant matrices)
* Clear and informative output

---

## 🎯 Assignment Requirements – Implementation Summary

Implement two functions:

* `jacobi(...)`
* `gauss_seidel(...)`

Each method:

* Receives matrix **A** and vector **b**
* Checks for **diagonal dominance**
* Performs **row pivoting** if needed
* Attempts convergence even if dominance is not achieved

Output includes:

* Iteration-by-iteration results
* Convergence status
* Number of iterations

Stopping condition:

```python
TOL = 1e-5
```

No user input required – system defined in code

---

## Mathematical Background

We solve a linear system:

[
Ax = b
]

Where:

* ( A ) – coefficient matrix
* ( x ) – unknown vector
* ( b ) – constants vector

---

### Jacobi Method

Each iteration computes:

[
x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j \ne i} a_{ij} x_j^{(k)} \right)
]

* Uses **only previous iteration values**
* Always computes a completely new vector

---

### Gauss-Seidel Method

Each iteration computes:

[
x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j < i} a_{ij} x_j^{(k+1)} - \sum_{j > i} a_{ij} x_j^{(k)} \right)
]

* Uses **updated values immediately**
* Typically converges faster than Jacobi

---

##  System Definition

The system is hardcoded as required:

```python
matrixA = [[4,2,0],
           [2,10,4],
           [0,4,5]]

vectorB = [[2],
           [6],
           [5]]
```

*(The implementation also supports larger matrices as shown in the code.)*

---

## 🔄 Algorithm Flow

1. **Matrix Initialization**
2. **Pivoting Attempt**

   * Reorders rows to improve diagonal dominance
3. **Diagonal Dominance Check**
4. **Iterative Solution**

   * Jacobi / Gauss-Seidel
5. **Convergence Check**
6. **Result Output**

---

## 🔍 Diagonal Dominance

A matrix is diagonally dominant if:

[
|a_{ii}| > \sum_{j \ne i} |a_{ij}|
]

### Behavior:

*  If satisfied → algorithm expected to converge
*  If not → algorithm still runs

  * If converges → prints:

    > "Despite no diagonal dominance, the system converged"
  * If not → prints:

    > "The system did not converge"

---

## 🔁 Iterative Solver Design

The function `iterative_solver(...)`:

* Handles both methods dynamically
* Avoids code duplication
* Uses **NumPy norm** for convergence check:

```python
np.linalg.norm(X_new - X_old) < eps
```

---

## 📊 Output Format

For each iteration:

```
Iteration k: X = [...]
```

Final output includes:

* Solution vector
* Number of iterations
* Convergence status

---

##  Code Structure

```
├── task()                      # Main program
├── iterative_solver()          # Controls iteration process
├── jacobi()                    # Jacobi step
├── gauss_seidel()              # Gauss-Seidel step
├── pivoting()                  # Row reordering
├── findmaxrow()                # Helper for pivoting
├── changeroworder()            # Swap rows
├── diagonally_dominant()       # Check dominance
```

---

## ⚙️Installation & Setup

Get the project up and running in under a minute:

1. Clone the repository
git clone <your-repo-link>
cd <project-folder>
2. Install dependencies
pip install numpy
3. Run the project
python main.py
```

---

## 💡 Design Decisions

* Separation between iteration logic and method implementation
* Reuse of helper functions to prevent redundancy
* Use of NumPy for numerical stability
* Support for non-ideal matrices (robustness)

---

## ⚠️ Limitations

* Convergence is not guaranteed for all matrices
* Pivoting improves but does not ensure diagonal dominance
* No dynamic user input (as required)

---

## Conclusion

This project demonstrates:

* Strong understanding of **numerical linear algebra**
* Correct implementation of **iterative solvers**
* Ability to handle **real-world edge cases**
* Writing **clean, modular, and maintainable Python code**

---

