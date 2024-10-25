
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
