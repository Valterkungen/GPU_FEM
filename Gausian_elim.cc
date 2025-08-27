#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

void partial_pivot(vector<vector<double> >& A, int n) {
    for (int i = 0; i < n; i++) {
        int pivot_row = i;
        for (int j = i+1; j < n; j++) {
            if (abs(A[j][i]) > abs(A[pivot_row][i])) {
                pivot_row = j;
            }
        }
        if (pivot_row != i) {
            for (int j = i; j <= n; j++) {
                swap(A[i][j], A[pivot_row][j]);
            }
        }
        for (int j = i+1; j < n; j++) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k <= n; k++) {
                A[j][k] -= factor * A[i][k];
            }
        }
    }
}

void back_substitute(const vector<vector<double> >& A, vector<double>& x) {
    int n = A.size();
    for (int i = n-1; i >= 0; i--) {
        double sum = 0;
        for (int j = i+1; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        x[i] = (A[i][n] - sum) / A[i][i];
    }
}

void gaussian_elimination(vector<vector<double> >& A, vector<double>& x) {
    int n = A.size();
    partial_pivot(A, n);
    back_substitute(A, x);
}

int main() {
    vector<vector<double> > A = {
        {3.0, 2.0, -4.0, 3.0},
        {2.0, 3.0, 3.0, 15.0},
        {5.0, -3, 1.0, 14.0}
    };
    int n = A.size();
    vector<double> x(n);

    gaussian_elimination(A, x);

    cout << "Solution for the system:\n";
    for (int i = 0; i < n; i++) {
        cout << "x[" << i << "] = " << x[i] << endl;
    }

    return 0;
}
