#include <iostream>

using namespace std;
 
int main () {
    cout<<"HELLO world ok";
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