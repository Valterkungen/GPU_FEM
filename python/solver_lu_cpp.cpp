#include <iostream>
#include <cmath>

using namespace std;
/*
void LUfactorisation(int n, float *L, float *U, float *A){
    for(int i = 0; i < n; i++){
        L[i*n + i] = 1.0; // set diagonal of L to 1
        for(int j = i; j<n; j++){
            float sum = 0.0;
            for(int k = 0; k < i; k++){
                sum += L[i * n + k] * U[k * n + j]; 
            }
            U[i * n + j] = A[i * n + j] - sum; //Computes U[i,j]
        }
        for(int j = i + 1; j < n; j++ ){
            float sum = 0.0;
            for(int k = 0; k < i; k++){
                sum += L[j*n + k] * U[k*n + i];
            }
            L[j*n + i] = (A[j*n + i] - sum) / U[i*n + i]; //Computes L[i,j]
        }
    } 
}*/

void LUfactorisation(int n, float L[3][3], float U[3][3], float A[3][3]) {
    for (int i = 0; i < n; i++) {
        L[i][i] = 1.0; // set diagonal of L to 1
        for (int j = i; j < n; j++) {
            float sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += L[i][k] * U[k][j];
            }
            U[i][j] = A[i][j] - sum; // Computes U[i,j]
        }
        for (int j = i + 1; j < n; j++) {
            float sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += L[j][k] * U[k][i];
            }
            L[j][i] = (A[j][i] - sum) / U[i][i]; // Computes L[j,i]
        }
    }
}

/*
float* forward_substitution(int n, float *L, float *b){
    float *y = new float[n]; 
    for(int i = 0; i < n; i++){
        float sum = 0.0;
        for(int j = 0; j < i; j++)
            sum += L[i*n + j] * y[j];
        y[i] = (b[i] - sum) / L[i*n + i];
    }
    return y;
}*/

float* forward_substitution(int n, float L[3][3], float b[3]) {
    float *y = new float[n];  // Dynamically allocate memory for the result vector y
    for (int i = 0; i < n; i++) {
        float sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];  // Access elements using two-dimensional indexing
        }
        y[i] = (b[i] - sum) / L[i][i];  // Access diagonal element of L directly
    }
    return y;  // Return the dynamically allocated array
}

/*float* backward_substitution(int n, float *U, float *y){
        float *x = new float[n];
        for(int i = n-1; i >= 0; i--){
            float sum = 0.0; 
            for(int j = i+1; j < i; j++) {
                sum += U[n*(n-1) - i*n + j] * y[n-j];
            x[i] = (y[i] - sum) / U[i*(n-1) + i];
            }
        }
    return x;
}*/

float* backward_substitution(int n, float U[3][3], float *y) {
    float *x = new float[n];  // Dynamically allocate memory for the result vector x
    for (int i = n - 1; i >= 0; i--) {
        float sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += U[i][j] * x[j];  // Correct indexing to access elements in U and y
        }
        x[i] = (y[i] - sum) / U[i][i];  // Directly access the diagonal element U[i][i]
    }
    return x;  // Return the dynamically allocated array
}


void printMatrix(const string& name, float matrix[3][3], int rows, int cols) {
    cout << name << ":\n";
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}


void printVector(const string& name, float vector[3], int rows){
    cout << name << ":\n";
    for(int i = 0; i<rows; i++){
        cout << vector[i] <<" ";
    }
    cout << endl;
}


int main() {

    int n = 3;
    float A_cpu[3][3] = {
        {4, -1, 1},
        {-1, 4, -2},
        {1, -2, 4}
    }; // Example matrix

    // Initializing L_cpu with zeros
    float L_cpu[3][3] = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}
    };

    // Initializing U_cpu with zeros
    float U_cpu[3][3] = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}
    };

    // b_cpu, x0_cpu, and x remain one-dimensional arrays
    float b_cpu[3] = {12, -1, 5}; // Example vectors
    float x0_cpu[3] = {1, 2, 1};  // Initial guess for some iterative methods
    float x[3] = {0, 0, 0};       // Resultant vector, initialized to zero


    LUfactorisation(n, L_cpu, U_cpu, A_cpu);
    printMatrix("L_cpu", L_cpu, n, n);
    printMatrix("U_cpu", U_cpu, n, n);

    float* y_vector = forward_substitution(n, L_cpu, b_cpu);
    printVector("Solution y", y_vector, n);

    float* x_vector = backward_substitution(n, U_cpu, y_vector);
    printVector("Solution x", x_vector, n);



}
