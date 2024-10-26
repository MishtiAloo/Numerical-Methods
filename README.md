# Numerical-Methods

**HEAD Console application for solving different types of mathematical equations using numerical methods.**

=======

Console application for solving different types of mathematical equations using numerical methods.  
At first, we are offered multiple options from which the user has to choose which type of system or problem they want to deal with,  
i.e., which type of equation they are trying to solve.

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

- **Initialization**: Guess the initial solution vector solution.
- **Iterative Update**: Update each variable, one after another at each step, using the remaining variables.
- **Convergence Check**: Stop iterating when the maximum change in the variable values between two consecutive iterations is smaller than some tolerance level.
- **Coefficients**: 2D vector containing coefficients of every equation.
- **Constants**: Constants appearing on the right-hand side of the equations.
- **Current and Previous**: Vectors holding the current and previous values of the variables.  
  That is, the threshold of each variable value at which the solution is considered converged.

---

### **Gaussian Elimination**:

The Gaussian elimination algorithm converts a system of linear equations into upper triangular matrix form using row operations. It replaces all entries below each pivot—the first nonzero element in each row—by the zero gotten by subtracting off a multiple of the pivot row from each subsequent row, working from top to bottom. This continues until all rows beneath the pivot have zeros in that column to yield an upper triangular matrix. This is then solved for each variable by back-substitution in a sequence from the bottom row upwards, substituting each solved variable into the rows above to work out the remaining unknowns.

**Method Explanation**:

- **gauss()**: The gauss() method contains the implementation of the Gauss Elimination Algorithm.
- **pivotSwaper()**: The pivotSwaper was used so that there occur no situation such 0 in the diagonal element.
- **printV()**: This method was used to print the 2D Vector.

---

### **Gauss-Jordan Elimination**:

Gauss-Jordan elimination is the extension of Gaussian elimination that fully reduces the matrix to reduced row echelon form, RREF, where each pivot element is 1 and all elements above and below the pivot are zeros. Once an upper triangular form is reached, it works its way from the bottom up, making zeros above each pivot by subtracting off appropriately scaled versions of each row. This results in a diagonal matrix whereby each row carries the value of only one variable, hence giving the direct solution of all the unknowns without the use of back-substitution.

**Method Explanation**:

- **gauss()**: The gauss() method contains the implementation of the Gauss Elimination Algorithm.
- **jordan()**: The Jordan() method contains the implementation of the Gauss-Jordan Elimination Algorithm.
- **pivotSwaper()**: The pivotSwaper was used so that there occur no situation such 0 in the diagonal element.
- **printV()**: This method was used to print the 2D Vector.

---

### **LU Factorization**:

The LU factorization of the given square matrix `A` expresses `A` as the product of two matrices: a lower triangular matrix `L` and an upper triangular matrix `U`, such that `A = LU`. In the process, `L` has ones on its diagonal while its elements are below the diagonal, and `U` has elements on and above the diagonal. It does this by iteratively eliminating entries below the main diagonal, scaling and subtracting rows while storing the scaling factors in `L` and the results in `U`. This allows one to solve the system `Ax = b` by first solving the system `Ly = b` using forward substitution and then solving `Ux = y` with backward substitution to obtain the solution vector `x`.

**Method Explanation**:

- **LU_Factorization()**: The LU_Factorization() method contains the implementation of the Gauss Elimination Algorithm.
- **printV()**: This method was used to print the 2D Vector.

# Non-Linear Equations:

### **Bisection Method**:

- **Enter the Polynomial Equation**: In this section, the user should define a polynomial function of some degree. The user will input coefficients, which will be stored in a vector `coefficient` where the highest degree goes first.
- **Specify Tolerance and Initial Interval**:
  - **Tolerance**: Acceptable error for approximating the root.
  - **First, pos and neg**: Define an interval [a, b] within which f(x) takes on values with opposite signs, meaning there exists a root between them.
- **Testing the Sign Condition**: This tests if the signs of f(pos) and f(neg) differ. If they do not, there is no guarantee that the algorithm will have a root within that interval.
- **Find Midpoint**: In each iteration, find the midpoint. This midpoint becomes a candidate to be the root by default.
  - **Evaluation**: The function f(x) is evaluated.
  - **Interval Bisection**: Based on the sign of f(x), either `pos` or `neg` is updated to the midpoint, halving the size of the interval.
- **Checking for Convergence**: The iterations continue until |f(x)| ≤ tolerance. When the absolute value of f(x) is less than the tolerance value, the program stops, and x is printed out as the approximate root.

---

### **False Position Method**:

- **Input the Polynomial Equation**: This section allows the user to specify a polynomial function of some order by providing the constants. The constants are stored in a vector called `coefficients`, with the highest order first.
- **Define Tolerance and Initial Interval**:
  - **Tolerance**: Defines the acceptable error in approximating the root.
  - **pos and neg**: An interval [a, b] for which the function f(x) has values with opposite signs, guaranteeing a root is between them.
- **Check the Sign of Condition**: This checks that f(pos) and f(neg) have different signs; otherwise, the algorithm can't be guaranteed to converge to a root within that interval.
- **Approximate Root x**: Function `FV(pos, neg, coefficients)` calculates the new point x using the formula for False Position:
  
  **x = b - (f(b) * (a - b)) / (f(a) - f(b))**
  
  Where `a` and `b` are the current values of `neg` and `pos` respectively.
- **Review f(x)**: Evaluate function value f(x) and check which bound gets updated based on the sign of f(x).
  - **Update Interval**: Depending on the positivity or negativity of f(x), either `pos` or `neg` gets updated to the midpoint, reducing the interval size.
- **Check for Convergence**: The loop continues until |f(x)| ≤ tolerance. When the absolute value of f(x) falls within the tolerance, the algorithm stops, and x is represented as the approximate root.

---

# Newton-Raphson Method for Finding Roots

The Newton-Raphson method is an iterative method used in finding the roots of a function f(x). The code prompts the user to enter the degree and coefficients of the polynomial under consideration, followed by successive approximations in determining the equation's roots using the formula guiding the Newton-Raphson method. Convergence is checked when the difference between successive approximations falls below a certain tolerance. This process is repeated until all the roots are computed.

---

### **Formula:**
The method uses successive approximations starting from an initial guess `X_0` by the formula:

**X_{n+1} = X_n - f(X_n) / f'(X_n)**

where **f'(X_n)** is the derivative of **f(x)**. Iterations are performed until the difference between successive approximations is less than a given tolerance. The method converges fast, provided the initial guess is close enough to the actual root. However, in case **f'(x) = 0** at any stage of the iterations, this method will fail because of division by zero.

---

### **1. Input:**
The program gets the **degree** and **coefficients** of the polynomial from the user.

---

### **2. Functions:**
- **fx**: Calculation of polynomial value for given `x`.
- **f_x**: The derivative of the polynomial for given `x`.
- **interval_Calculator**: In accordance with polynomial coefficients, an interval to search for roots will be determined.
- **find_ab**: It searches for an interval `[a, b]` where a root may exist by investigating the changes in the sign of the function.

---

### **3. Newton-Raphson Iteration:**
- Taking an initial guess (**Xo**) from the interval `[a, b]`, updates **Xo** with the formula of the Newton-Raphson method:
  
  **Xn = Xo - fx(Xo) / f_x(Xo)**

- Continues the iteration until the difference between the new guess (**Xn**) and the old guess (**Xo**) is less than a predefined tolerance (**Et**), meaning it has converged to a root.

---

### **4. Multiple Roots:**
- The loop executes until all roots—equal to the degree of the polynomial—are calculated.

---

### **Error Handling:**
- If a division by zero takes place—in other words, the derivative **f_x(Xo)** equals zero—the program is interrupted and a **"Divide by zero"** message is flashed.


---

# Secant Method for Finding Roots

The Secant method is a numerical method to find the roots of a non-linear equation. It iteratively approximates the root by using two initial guesses and a formula based on the function’s values at those points. Unlike Newton’s method, it doesn’t need the derivative calculation. It generates a new point using the secant line (the line between the two points) and repeats the process until the root is found with a specified tolerance. It’s faster than the Bisection method but can fail if the initial guesses are not close to the actual root.

---

### **Function: `Secant()`**

The `Secant()` function finds the roots of a polynomial equation using the Secant method.

### **1. Polynomial Setup**

- The user inputs the **degree** and **coefficients** of the polynomial.

### **2. Helper Functions**

- **fx**: Evaluates the polynomial at a given value `n`.
- **interval_Calculator**: Calculates the search range for the root based on the polynomial’s coefficients.
- **find_ab**: Finds initial values `a` and `b` where the function changes sign (i.e., between which a root lies).

### **3. Root Calculation**

- The method iterates to find each root, using the Secant formula to approximate the next `Xn`.
- The method stops when the difference between consecutive `Xn` is smaller than a tolerance value `Et`.

### **4. Multiple Roots**

- The function will keep searching for roots until the number of roots is equal to the degree of the polynomial.


# Differential Equation:

### **Runge-Kutta Method**:

**Steps**:
1. **Characterize Beginning Conditions or Intervals**.
2. **Begin with the starting condition** y(x0) = y0 and step size h.
3. **Calculate the Four formulas**:
   - For each step xn, calculate the kn:
     - k1 = h * f(x(n), y(n))
     - k2 = h * f(x(n) + (h/2), y(n) + (k1/2))
     - k3 = h * f(x(n) + (h/2), y(n) + (k2/2))
     - k4 = h * f(x(n) + h, y(n) + k3)

4. **Update Yn, Xn**:
   - y(n+1) = y(n) + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)

5. **Rehash for Each Xn**:
   - Set xn+1 = xn + h and rehash the method.

**Step size- Accuracy tradeoff**:
- If a higher value of h is chosen, the method fails to measure accurately in steep curve changes.
- If the step size is chosen to be very small, the method takes much longer to compute the points.

---

### **Matrix Inversion**:

Matrix inversion by **Gaussian-Jordan elimination** is a step-by-step process that transforms any given square matrix A into an identity matrix while simultaneously transforming an identity matrix of the same size into the inverse of A. 

- This is done by appending an identity matrix onto matrix A, and then using row operations on the resultant augmented matrix to turn all entries below and above each pivot element into zero. 
- The process starts with the reduced row echelon form of the upper triangular form of A.
- When A has been reduced to the identity matrix, the appended identity matrix will be A<sup>-1</sup>.
- The approach may involve interchanging rows if necessary to avoid division by zero, ensuring that A is non-singular (i.e., it has an inverse).

**Method Explanation**:
- **gaussI()**: The `gaussI()` method contains the implementation of the Gauss Elimination Algorithm.
- **jordanI()**: The `jordanI()` method contains the implementation of the Gauss-Jordan Elimination Algorithm.
- **pivotSwaper()**: The `pivotSwaper()` method ensures that there are no situations where 0 occurs in the diagonal element.
- **printVI()**: This method is used to print the 2D Vector.

