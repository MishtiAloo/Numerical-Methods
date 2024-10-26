# Numerical-Methods
<<<<<<< HEAD
Console application for solving different types of mathematical equations using numerical methods.
=======
Console application for solving different types of mathematical equations using numerical methods.
At first we are offered multiple options from which we the user have to choose which type of system or problem we want to deal with 
i.e. which type of equation we are trying to solve


Gaussian Elimination

The Gaussian elimination algorithm converts a system of linear equations into upper triangular matrix form using row operations. It replaces all entries below each pivot-the first nonzero element in each row-by the zero gotten by subtracting off a multiple of the pivot row from each subsequent row, working from top to bottom. This continues until all rows beneath the pivot have zeros in that column to yield an upper triangular matrix. This is then solved for each variable by back-substitution in a sequence from the bottom row upwards, substituting each solved variable into the rows above to work out the remaining unknowns.

Method Explanation 

gauss() : The gauss() method contains the implementation of the Gauss Elimination Algorithm. 
pivotSwaper(): The pivotSwaper was used so that there occur no situation such 0 in the diagonal element.
printV(): This method was used to print the 2D Vector


GaussJordan Elimination

Gauss-Jordan elimination is the extension of Gaussian elimination that fully reduces the matrix to reduced row echelon form, RREF, where each pivot element is 1 and all elements above and below the pivot are zeros. Once an upper triangular form is reached, it works its way from the bottom up, making zeros above each pivot by subtracting off appropriately scaled versions of each row. This results in a diagonal matrix whereby each row carries the value of only one variable, hence giving the direct solution of all the unknowns without the use of back-substitution.

Method Explanation 

gauss() : The gauss() method contains the implementation of the Gauss Elimination Algorithm. 
jordan() : The Jordan() method contains the implementation of the Gauss-Jordan Elimination Algorithm.
pivotSwaper(): The pivotSwaper was used so that there occur no situation such 0 in the diagonal element.
printV(): This method was used to print the 2D Vector

LU Factorization

The LU factorization of the given square matrix A expresses A as the product of two matrices: a lower triangular matrix L and an upper triangular matrix U, such that A = LU. In the process, L has ones on its diagonal while its elements are below the diagonal, and U has elements on and above the diagonal. It does this by iteratively eliminating entries below the main diagonal, scaling and subtracting rows while storing the scaling factors in L and the results in U. This allows one to solve the system Ax = b by first solving the system Ly = b using forward substitution and then solving Ux = y with backward substitution to obtain the solution vector x.

Method Explanation 

LU_Factorization() : The LU_Factorization() method contains the implementation of the Gauss Elimination Algorithm. 

printV(): This method was used to print the 2D Vector


Matrix Inversion

Matrix inversion by Gaussian-Jordan elimination is a step-by-step process that will eventually transform any given square matrix A into an identity matrix and, simultaneously, do the reverse with another identity matrix of the same size, yielding the inverse of A. It does this by appending an identity matrix on matrix A, and then uses row operations on the resultant augmented matrix to turn all entries below and above each pivot element into zero, starting with the reduced row echelon form of the upper triangular form of A. When A has been reduced to the identity matrix, the appended identity matrix will be A-1. This approach employs the method of interchanging rows if need be to avoid division by zero and essentially takes advantage of the fact that A is non-singular, that is, it has an inverse.

gaussI() : The gaussI() method contains the implementation of the Gauss Elimination Algorithm. 

jordanI() : The jordanI() method contains the implementation of the Gauss-Jordan Elimination Algorithm.

pivotSwaper(): The pivotSwaper was used so that there occur no situation such 0 in the diagonal element.

printVI(): This method was used to print the 2D Vector


Jacobi Iterative Method:

An approximate solution vector  is obtained. Many times the zero vector is used as an approximate solution vector.

-Iteration : Each variable of the system at each iteration uses values from the previous iterations.
-Checking for Convergence: Continuing iterations past convergence where the difference in values for successive iterations is less than a pre-specified tolerance.
-input: Number of variables n in the system of linear equations.
-coefficients: A 2-D vector of floats to store coefficients of the linear equations as a matrix. 
-constants: A vector for storing the constants on the right side of the equations.
-current and previous: Vectors corresponding to the current and previous value of the variables in every iteration.
-Tolerance: Tolerance value of the error which defines whether the solution has converged.
-Error: This is a vector of the error in the value of each variable at each stage; this is needed to test for convergence.
-The diagonal dominance, i.e., for best convergence results, matrix A should be diagonally dominant.That, by itself, means every element on the diagonal has 
 to be greater than the summation of   all other elements' absolute value on the same row.
 In case it is not diagonally dominant, the Jacobi's method is not necessary to converge.

Gauss Seidal Method :

This is an iterative method to find successive approximations to the solution of a system of linear equations. This differs from Jacobi, where each variable is independently calculated with respect to others in an iteration; in the Gauss-Seidal method, it uses the most updated value in the iteration while updating each variable in place.
Initialization: Guess the initial solution vector solution
• Iterative Update: Update each variable, one after another at each step, using the remaining variables.
• Convergence Check: Stop iterating when the maximum change in the variable values between two consecutive iterations is smaller than some tolerance level.
coefficients: 2D vector containing coefficients of every equation
constants: constants appearing on the right-hand side of the equations
current and previous : vectors holding the current and previous values of the variables.
That is, the threshold of each variable value at which the solution is considered converged.

Bisection Method:

Enter the Polynomial Equation: In this section, the user should define a polynomial function of some degree. The user will input coefficients, which afterwards will be stored in a vector coefficient in which the highest degree goes first.
Specify Tolerance and Initial Interval:
tolerance acceptable error for approximating the root.
First, pos and neg define an interval [a, b] within which f(x) takes on values with opposite signs. That is to say, there exists a root between them. Testing the Sign Condition: This tries to test if the signs of f( pos) and f(neg) differ. If they do not, there is no guarantee that the algorithm will have a root within that interval. Find Midpoint: In each iteration, find the midpoint. That by default becomes a candidate to be the root.
-Evaluation: the function f(x) is evaluated.
- Interval bisection: Based on the sign of f(x), either pos or neg is updated to the midpoint and, in effect halves the size of the interval.
Checking for Convergence: The iterations continue until |f(x)|≤tolerance.
When the absolute value of f(x) is less than the value of the tolerance, the program is stopped and
x is printed out as the approximate root.
False Position Method :

• Input the Polynomial Equation: This section allows the user to specify a polynomial function of some order by providing the constants. The constants will be stored in a vector called coefficients, with the highest order first.
• Define Tolerance and Initial Interval:
Tolerance: defines the acceptable error in approximating the root.
pos and neg: an interval [a,b] for which the function f(x) has values with opposite signs, guaranteeing a root is between them.
Check the sign of condition:
This checks that f (pos) and f (neg) have different signs else algorithm can't be guaranteed to converge to a root within that interval.
Approximate Root x: function FV(pos, neg, coefficients) calculates new point \\x using formula for False Position
x= b -(f(b)*(ab)) / f(a)-f(b)
Where a and b are the current values of neg and pos respectively.
Review: f(x): Evaluate function value f(x) and check sign of f(x) tells which bound is updated.
-Update Interval: Depending on the positivity or negativity of f(x), either pos or neg gets updated to midpoint, therefore interval size gets reduced by half.
Check for Convergence: Loop is continued until ∣f(x)∣≤tolerance.
When the absolute value of f(x) falls into the tolerance, the algorithm stops, and
x is represented as the approximate root.

3. Newton-Raphson Strategy
The Newton-Raphson Strategy may be a more effective and speedier strategy that begins with an introductory figure and iteratively moves forward it. This strategy employments the subsidiary of the work to discover progressively way better approximations to the root.

The equation is:
X_i = X_i-1 - ( f(X_i-1) / f'(X_i-1) )

Steps:
a.Begin with an beginning figure (X_0) where (X_0 = (a + b) / 2).
b. Compute the another guess utilizing the equation.
c. Rehash until the work esteem (f(X_i)) is adequately near to zero.

Secant Strategy:
The Secant Strategy is comparable to the Newton-Raphson strategy but does not require the computation of the subsidiary. Instep, it approximates the subsidiary by utilizing two past focuses. 

The equation is:
X_n+1 = X_n - ( f(X_n) * (X_n - X_n-1) ) / f(X_n) - f(X_n-1)

Steps:
1. Begin with an introductory figure X_n and X_n-1 where X_n = a and X_n-1 = b.
2. Compute the another estimation utilizing the equation.
3. Rehash until the work esteem f(X_n+1) is adequately near to zero.

Runge-Kutta method:
Steps:
1. Characterize Beginning Conditions or intervals
2. Begin with the starting condition y(x0) = y0 and step estimate h.
3. Calculate the Four formulas:
    For each step xn, calculate the kn:
    k1=h*f(x(n),y(n))
    k2=h*f(x(n)+(h/2),y(n)+(k1/2))
    k3=h*f(x(n)+(h/2),y(n)+(k2/2))
    k4=h*f(x(n)+h,y(n)+k3)

4. update Yn, Xn:
    y(n+1) = y(n) + (1/6)(k1 + 2k2 + 2*k3 + k4)

5. Rehash for Each Xn:
    Set xn+1 = xn + h and rehash the method.

Step size- Accuracy tradeoff:
If higher value of h is chosen then the methods fails to measure accurately in the steep curve changes, if step size is chosen to be very small, the method takes much time to give the points.
 



