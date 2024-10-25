#include <bits/stdc++.h>
// #include <iostream>
// #include <vector>
// #include <cmath>

using namespace std;


vector<float> backSubstitution(vector<vector<float>> A, int n);
vector<float> extractSolution(const vector<vector<float>>& A, int n);
vector<vector <float> > gauss(vector<vector <float> > A, int n);
vector<vector <float> > jordan(vector<vector <float> > A, int n);
void pivotSwaper(vector<vector <float> > &A, int n, int i);
int sign(float a);
void printV(vector<vector <float> > A);


 bool isDiagonallyDominant(const vector<vector<float>>& coefficients, int n) ;
 void Gauss_Seidel();
 void Jacobi_Iteration();
int main () {
    cout<<"HELLO world ";
    int choice;
    while (1) {
        cout << "What do you want to solve?\n\n1.Linear system of equations\n2.Non-Linear equation\n3.Differential equation\n4.Matrix inversion\n5. Exit\n\nYour choice: ";
        cin >> choice;        

        if (choice < 1 || choice > 5) {
            cout << "\n\nInvalid input!!! Try again.\n\n";
        }
        else break;
    }

    switch (choice) {
        int method;
    case 1: {
        
        while (1) {
            cout << "Chose your method\n\n1. Jacobi iterative method\n2. Gauss-Seidel iterative method\n3. Gauss elimination\n4. Gausss-Jordan elimination\n5. LU factorization\n6. Default\n7. Exit\n\nYour choice: ";
            cin >> method;

            if (choice < 1 || choice > 7) {
                cout << "\n\nInvalid input!!! Try again.\n\n";
            }
            else break;
        }

        switch (method) {
        case 1: {
            // Jacobi
            Jacobi_Iteration();
            break;
        }
        case 2: {
            // Gauss Seidel
            Gauss_Seidel();
            break;
        }
        case 3: {
            // GE
            vector<vector<float>> A; 
            int n; 
            cout << "Enter the number of variables: ";
            cin >> n;
            
            cout << "\nEnter the coefficients and constant of each equation in matrix format:\n\n";
            
            // Input matrix in augmented form, including the constants in the last column
            for(int i = 0; i < n; i++) {
                vector<float> tmp; 
                for(int j = 0; j <= n; j++) { 
                    float t;
                    cin >> t;
                    tmp.push_back(t); 
                }
                A.push_back(tmp); 
            }
            
            // Ensure the pivot elements are non-zero by swapping if necessary
            for(int i = 0; i < n; i++) {
                pivotSwaper(A, n, i); 
            }
            
            A = gauss(A, n);
            // Use back-substitution to solve for each variable in the upper triangular form
            vector<float> solution = backSubstitution(A, n);
            
            
            cout << "\nSolution after Gaussian elimination:\n";
            
            for(int i = 0; i < n; i++) {
                cout << "x" << i + 1 << "\t=\t" << solution[i] << endl;
            }

            break;
        }
        case 4: {
            // GJ
            vector<vector<float>> A; 
            int n; 
            cout << "Enter the number of variables: ";
            cin >> n;
            
            cout << "\nEnter the coefficients and constant of each equation in matrix format:\n\n";
            
            // Input matrix in augmented form, including the constants in the last column
            for(int i = 0; i < n; i++) {
                vector<float> tmp; 
                for(int j = 0; j <= n; j++) { 
                    float t;
                    cin >> t;
                    tmp.push_back(t); 
                }
                A.push_back(tmp); 
            }
            
            
            for(int i = 0; i < n; i++) {
                pivotSwaper(A, n, i); 
            }
            
            
            A = gauss(A, n);
            
            vector<float> solution ;

            A = jordan(A, n);
            
            // Extract the final solution directly from the reduced row echelon form
            solution = extractSolution(A, n);
            
            cout << "\nSolution after Gauss-Jordan elimination:\n";
            for(int i = 0; i < n; i++) {
                cout << "x" << i + 1 << "\t=\t" << solution[i] << endl;
            }
    
            break;
        }
        case 5: {
            // LU
            break;
        }
        case 6: {
            // Whatever
            break;
        }
        case 7: {
            return -1;
        }
        }

        break;
    }

    case 2: {
        
        while (1) {
            cout << "Chose your method\n\n1. Bisection method\n2. Fakse Position method\n3. Secant\n4. Newton Raphson method\n5. Default\n6. Exit\n\nYour choice: ";
            cin >> method;

            if (choice < 1 || choice > 6) {
                cout << "\n\nInvalid input!!! Try again.\n\n";
            }
            else break;
        }

        switch (method) {
        case 1: {
            // BS
            break;
        }
        case 2: {
            // FP
            break;
        }
        case 3: {
            // Secant
            break;
        }
        case 4: {
            // NR
            break;
        }
        case 5: {
            // Whatever
            break;
        }
        case 6: {
            return -1;
        }
        }

        break;
    }

    case 3: {
        // RK
        break;
    }

    case 4: {
        // Inveerse
        break;
    }

    case 5: {
        return -1;
    } 
    }
}


// Gaussian Elimination to transform matrix into upper triangular form
vector<vector<float>> gauss(vector<vector<float>> A, int n) {
    for(int i = 0; i < n - 1; i++) { 
        for(int j = n - 1; j > i; j--) { 
            float a = A[j][i]; 
            int l = j - 1;    
            if (a == 0) continue; 

            // Find a non-zero element above as the pivot if current row is zero
            while (A[l][i] == 0) {
                l--;
                if (l < 0) { // If no non-zero pivot found, matrix is singular
                    cout << "Invalid Matrix.\n";
                    return A;
                }
            }

            // Scale rows to create zeros below the pivot
            float b = A[l][i];
            for(int k = 0; k <= n; k++) {
                A[j][k] = A[j][k] * b - A[l][k] * a;
                if (fabs(A[j][k]) < 1e-6) A[j][k] = 0; // Prevent floating-point errors
            }
        }
    }
    return A;
}

// Gauss-Jordan Elimination to transform matrix into reduced row echelon form
vector<vector<float>> jordan(vector<vector<float>> A, int n) {
    for(int i = n - 1; i > 0; i--) { 
        for(int j = 0; j < i; j++) { 
            float a = A[j][i]; 
            int l = j + 1;     
            if (a == 0) continue; 

            // Find a non-zero element below as the pivot if needed
            while (A[l][i] == 0) {
                l++;
                if (l > n - 1) { // If no non-zero pivot found, matrix might be singular
                    cout << "Invalid Matrix.\n";
                    return A;
                }
            }

            // Scale rows to create zeros above the pivot
            float b = A[l][i];
            for(int k = 0; k <= n; k++) {
                A[j][k] = A[j][k] * b - A[l][k] * a;
                if (fabs(A[j][k]) < 1e-6) A[j][k] = 0; // Prevent floating-point errors
            }
        }
    }
    return A;
}


void printV(vector<vector<float>> A) {
    for(int i = 0; i < A.size(); i++) {
        for(int j = 0; j <= A.size(); j++) {
            cout << A[i][j] << "\t";
        }
        cout << endl;
    }
}

// Finds and swaps the pivot to avoid division by zero during elimination
void pivotSwaper(vector<vector<float>>& A, int n, int i) {
    if (A[i][i] != 0) return; 

    // Search for a row below with a non-zero element in the pivot column
    for (int j = i + 1; j < n; j++) {
        if (A[j][i] != 0) {
            swap(A[j], A[i]); 
            return;
        }
    }

    // If no valid pivot is found, the matrix might be singular
    cout << "No valid pivot found in column " << i << ". Matrix might be singular.\n";
}

// Returns the sign of a float: -1 for negative, 1 for non-negative
int sign(float a) {
    if (a < 0) return -1;
    else return 1;
}

// Performs back-substitution to solve for each variable in upper-triangular matrix form
vector<float> backSubstitution(vector<vector<float>> A, int n) {
    vector<float> x(n, 0); 
    for (int i = n - 1; i >= 0; i--) {
        x[i] = A[i][n]; 
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j]; 
        }
        x[i] /= A[i][i]; // Divide by the pivot element to get the solution for x[i]
    }
    return x;
}

// Extracts the final solution from the matrix in reduced row echelon form
vector<float> extractSolution(const vector<vector<float>>& A, int n) {
    vector<float> solution(n);
    for (int i = 0; i < n; i++) {
        solution[i] = A[i][n] / A[i][i]; 
    }
    return solution;
}

bool isDiagonallyDominant(const vector<vector<float>>& coefficients, int n) {
    for (int i = 0; i < n; i++) {
        float sum = 0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                sum += fabs(coefficients[i][j]);
            }
        }
        if (fabs(coefficients[i][i]) <= sum) {
            return false;  // Not diagonally dominant if diagonal element is not greater than sum
        }
    }
    return true;
}


 void Gauss_Seidel()
 {
     cout << "Enter the number of variables : ";
    int n;
    cin >> n;

  

    vector<vector<float>> coefficients(n, vector<float>(n));
    vector<float> constants(n);
    vector<float> current(n, 0), previous(n, 0);
    float tolerance;

    cout << "Enter the tolerance level: ";
    cin >> tolerance;

    // Input coefficients and constants for each equation
    for (int i = 0; i < n; i++) {
        cout << "Enter coefficients and constant term for equation " << i + 1 << ": ";
        for (int j = 0; j < n; j++) {
            cin >> coefficients[i][j];
        }
        cin >> constants[i];
    }


    if (!isDiagonallyDominant(coefficients, n)) {
        cout << "Warning: The matrix is not diagonally dominant. The Jacobi method may not converge.\n";

    }
    else{
    int iteration = 1;
    vector<float> errors(n, 0.0);

    // (Gauss-Seidel)
    do {
        cout << "Iteration: " << iteration << endl;
        for (int i = 0; i < n; i++) {
            float sum = constants[i];
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    sum -= coefficients[i][j] * current[j];
                }
            }
            previous[i] = current[i];
            current[i] = sum / coefficients[i][i];
            errors[i] = fabs(current[i] - previous[i]);
            cout << "Variable " << i + 1 << " = " << current[i] << endl;
        }
        iteration++;
    } while (*max_element(errors.begin(), errors.end()) > tolerance);

    cout << "\nRoots:- \n";
    for (int i = 0; i < n; i++) {
        cout << "Variable " << i + 1 << " = " << current[i] << endl;
    }
 }
 }
 void Jacobi_Iteration()
 {
      cout << "Enter the number of variables : ";
    int n;
    cin >> n;

  

    vector<vector<float>> coefficients(n, vector<float>(n));
    vector<float> constants(n);
    vector<float> current(n, 0), previous(n, 0);
    float tolerance;
     vector<float> errors(n, 0.0);

    cout << "Enter the tolerance level: ";
    cin >> tolerance;

    
    for (int i = 0; i < n; i++) {
        cout << "Enter coefficients and constant term for equation " << i + 1 << ": ";
        for (int j = 0; j < n; j++) {
            cin >> coefficients[i][j];
        }
        cin >> constants[i];
    }


    if (!isDiagonallyDominant(coefficients, n)) {
        cout << "Warning: The matrix is not diagonally dominant. The Jacobi method may not converge.\n";

    }
else {
    int iteration = 1;
    

    do {
    
        cout << "Iteration " << iteration << ":\n";

    
        for (int i = 0; i < n; i++) {
          

            float sum = constants[i];
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    sum -= coefficients[i][j] * previous[j];
                }
            }
            current[i] = sum / coefficients[i][i];

            
            errors[i] = fabs(current[i] - previous[i]);
           

            cout << "x" << i + 1 << " = " << fixed << setprecision(6) << current[i] << " ";
        }
        cout << endl;

        
        previous = current;
        iteration++;
    } while (*max_element(errors.begin(), errors.end()) > tolerance);

    cout << "\nSolution:\n";
    for (int i = 0; i < n; i++) {
        cout << "x" << i + 1 << " = " << fixed << setprecision(6) << current[i] << endl;
    }
}
 }