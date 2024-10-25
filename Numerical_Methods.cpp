#include <bits/stdc++.h>

using namespace std;
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
            break;
        }
        case 2: {
            // GS
            Gauss_Seidel();
            break;
        }
        case 3: {
            // GE
            break;
        }
        case 4: {
            // GJ
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