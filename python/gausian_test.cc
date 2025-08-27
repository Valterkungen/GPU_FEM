#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

void partial_pivot(vector<vector<double>>& A, int n) {
    for (int i = 0; i < n; i++) {
        int pivot_row = i;
        for (int j = i + 1; j < n; j++) {
            if (abs(A[j][i]) > abs(A[pivot_row][i])) {
                pivot_row = j;
            }
        }
        if (pivot_row != i) {
            for (int j = i; j <= n; j++) {
                swap(A[i][j], A[pivot_row][j]);
            }
        }
        for (int j = i + 1; j < n; j++) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k <= n; k++) {
                A[j][k] -= factor * A[i][k];
            }
        }
    }
}

void back_substitute(const vector<vector<double>>& A, vector<double>& x) {
    int n = A.size();
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        x[i] = (A[i][n] - sum) / A[i][i];
    }
}

void gaussian_elimination(vector<vector<double>>& A, vector<double>& x) {
    int n = A.size();
    partial_pivot(A, n);
    back_substitute(A, x);
}

vector<double> F(double t, const vector<double>& Z, const vector<vector<double>>& M, const vector<vector<double>>& A) {
    int n = Z.size() / 2;
    vector<double> X(Z.begin(), Z.begin() + n);
    vector<double> Y(Z.begin() + n, Z.end());
    vector<double> dXdt = Y;

    // Prepare matrix for Gaussian elimination
    vector<vector<double>> MA(M.size(), vector<double>(M[0].size() + 1));
    for (size_t i = 0; i < M.size(); ++i) {
        for (size_t j = 0; j < M[i].size(); ++j) {
            MA[i][j] = M[i][j];
        }
    }

    // Multiply A*X
    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < n; ++j) {
            sum += A[i][j] * X[j];
        }
        MA[i][M[i].size()] = -sum;
    }

    // Solve for dY/dt using Gaussian elimination
    vector<double> dYdt(n);
    gaussian_elimination(MA, dYdt);

    vector<double> dZdt;
    dZdt.insert(dZdt.end(), dXdt.begin(), dXdt.end());
    dZdt.insert(dZdt.end(), dYdt.begin(), dYdt.end());
    return dZdt;
}

vector<double> rk4_step(double t, const vector<double>& Z, double h, const vector<vector<double>>& M, const vector<vector<double>>& A) {
    vector<double> k1 = F(t, Z, M, A);
    vector<double> Z_temp = Z;
    for (size_t i = 0; i < Z.size(); ++i) {
        Z_temp[i] = Z[i] + 0.5 * h * k1[i];
    }
    vector<double> k2 = F(t + 0.5 * h, Z_temp, M, A);

    for (size_t i = 0; i < Z.size(); ++i) {
        Z_temp[i] = Z[i] + 0.5 * h * k2[i];
    }
    vector<double> k3 = F(t + 0.5 * h, Z_temp, M, A);

    for (size_t i = 0; i < Z.size(); ++i) {
        Z_temp[i] = Z[i] + h * k3[i];
    }
    vector<double> k4 = F(t + h, Z_temp, M, A);

    vector<double> next_Z = Z;
    for (size_t i = 0; i < Z.size(); ++i) {
        next_Z[i] = Z[i] + (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
    return next_Z;
}

int main() {
    int n = 4; // Size of the system
    vector<vector<double>> M = {{2, 1, 1, 0}, {4, 3, 3, 1}, {8, 7, 9, 5}, {6, 7, 9, 8}};
    vector<vector<double>> A = {{1, 0, 2, -1}, {1, 3, 1, 0}, {1, 2, 2, 1}, {3, 1, 1, 2}};
    vector<double> x = {1, 2, 3, 4};
    vector<double> Z(2 * n);
    copy(x.begin(), x.end(), Z.begin());
    fill(Z.begin() + n, Z.end(), 0); // Fill the second half with zeros

    double h = 0.01; // Time step
    double t = 0; // Initial time
    int n_steps = 100; // Number of RK4 steps

    for (int step = 0; step < n_steps; ++step) {
        Z = rk4_step(t, Z, h, M, A);
        t += h;
    }

    cout << "Results after RK4 integration:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "X[" << i << "] = " << Z[i] << ", Y[" << i << "] = " << Z[i + n] << endl;
    }

    return 0;
}
