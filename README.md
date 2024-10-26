# Numerical-Methods

**HEAD Console application for solving different types of mathematical equations using numerical methods.**

=======

Console application for solving different types of mathematical equations using numerical methods.  
At first, we are offered multiple options from which we the user have to choose which type of system or problem we want to deal with, i.e., which type of equation we are trying to solve.

---

## Linear Equations:

### **Jacobi Iterative Method**:

An approximate solution vector is obtained. Many times the zero vector is used as an approximate solution vector.

- **Iteration**: Each variable of the system at each iteration uses values from the previous iterations.
- **Checking for Convergence**: Continuing iterations past convergence where the difference in values for successive iterations is less than a pre-specified tolerance.
- **Input**: Number of variables `n` in the system of linear equations.
- **Coefficients**: A 2-D vector of floats to store coefficients of the linear equations as a matrix.
- **Constants**: A vector for storing the constants on the right side of the equations.
- **Current and Previous**: Vectors corresponding to the current and previous value of the variables in every iteration.
- **Tolerance**: Tolerance value of the error which defines whether the solution has converged.
- **Error**: This is a vector of the error in the value of each variable at each stage; this is needed to test for convergence.
- **Diagonal Dominance**: For best convergence results, matrix `A` should be diagonally dominant, meaning every element on the diagonal has to be greater than the summation of all other elements' absolute value on the same row. If it's not diagonally dominant, Jacobi's method does not necessarily converge.

---

### **Gauss-Seidal Method**:

This is an iterative method to find successive approximations to the solution of a system of linear equations. This differs from Jacobi, where each variable is independently calculated with respect to others in an iteration; in the Gauss-Seidal method, it uses the most updated value in the iteration while updating each variable in place.

- **Initialization**: Guess the initial solution vector solution
- **Iterative Update**: Update each variable, one after another at each step, using the remaining variables.
- **Convergence Check**: Stop iterating when the maximum change in the variable values between two consecutive iterations is smaller than some tolerance level.

---

### **Gaussian Elimination**

The Gaussian elimination algorithm converts a system of linear equations into upper triangular matrix form using row operations.

- **gauss()**: The gauss() method contains the implementation of the Gauss Elimination Algorithm.
- **pivotSwaper()**: The pivotSwaper was used so that there occur no situation such 0 in the diagonal element.
- **printV()**: This method was used to print the 2D Vector.

---

### **Gauss-Jordan Elimination**

Gauss-Jordan elimination is the extension of Gaussian elimination that fully reduces the matrix to reduced row echelon form, RREF.

- **gauss()**: The gauss() method contains the implementation of the Gauss Elimination Algorithm.
- **jordan()**: The Jordan() method contains the implementation of the Gauss-Jordan Elimination Algorithm.
- **pivotSwaper()**: The pivotSwaper was used so that there occur no situation such 0 in the diagonal element.
- **printV()**: This method was used to print the 2D Vector.

---

### **LU Factorization**

The LU factorization of the given square matrix `A` expresses `A` as the product of two matrices: a lower triangular matrix `L` and an upper triangular matrix `U`.

- **LU_Factorization()**: The LU_Factorization() method contains the implementation of the Gauss Elimination Algorithm.
- **printV()**: This method was used to print the 2D Vector.

---

## Non-Linear Equations:

### **Bisection Method**

- Enter the Polynomial Equation: In this section, the user should define a polynomial function of some degree. The user will input coefficients, which afterwards will be stored in a vector coefficient in which the highest degree goes first.
- Specify Tolerance and Initial Interval.
- Checking for Convergence: The iterations continue until |f(x)|≤tolerance.

---

### **False Position Method**

- Input the Polynomial Equation.
- Define Tolerance and Initial Interval.
- Approximate Root x: function FV(pos, neg, coefficients) calculates new point \\x using formula for False Position.
- Check for Convergence: Loop is continued until ∣f(x)∣≤tolerance.

---

### **Newton-Raphson Method**

The Newton-Raphson Strategy is a more effective and speedier strategy that begins with an introductory figure and iteratively moves forward.

---

### **Secant Method**

The Secant method is comparable to the Newton-Raphson method. Instep, it approximates the value by utilizing two past values.

---

## Differential Equations:

### **Runge-Kutta Method**

1. Characterize Beginning Conditions or intervals
2. Begin with the starting condition y(x0) = y0 and step estimate h.
3. Calculate the Four formulas:
    - For each step `xn`, calculate the `kn`:
        - `k1 = h*f(x(n),y(n))`
        - `k2 = h*f(x(n)+(h/2),y(n)+(k1/2))`
        - `k3 = h*f(x(n)+(h/2),y(n)+(k2/2))`
        - `k4 = h*f(x(n)+h,y(n)+k3)`
4. Update `Yn`, `Xn`: `y(n+1) = y(n) + (1/6)(k1 + 2k2 + 2*k3 + k4)`
5. Rehash for each `Xn`: Set `xn+1 = xn + h` and repeat the method.

---

## Matrix Inversion

Matrix inversion by Gaussian-Jordan elimination.

- **gaussI()**: The gaussI() method contains the implementation of the Gauss Elimination Algorithm.
- **jordanI()**: The jordanI() method contains the implementation of the Gauss-Jordan Elimination Algorithm.
- **pivotSwaper()**: The pivotSwaper was used so that there occur no situation such 0 in the diagonal element.
- **printVI()**: This method was used to print the 2D Vector.
