#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

vector<vector<double>> gaussianElimination(vector<vector<double>>& A) {
    int n = A.size();

    // Augment the matrix with an identity matrix
    vector<vector<double>> B(n, vector<double>(2 * n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            B[i][j] = A[i][j];
        }
        B[i][i + n] = 1.0;
    }

    // Applying Gauss Jordan Elimination
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                double ratio = B[j][i] / B[i][i];
                for (int k = 0; k < 2 * n; ++k) {
                    B[j][k] -= ratio * B[i][k];
                }
            }
        }
    }

    // Scaling all values to make diagonal elements 1
    for (int i = 0; i < n; ++i) {
        double divisor = B[i][i];
        for (int j = 0; j < 2 * n; ++j) {
            B[i][j] /= divisor;
        }
    }

    // Extracting the right half of the augmented matrix which is the inverse
    vector<vector<double>> result(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i][j] = B[i][j + n];
        }
    }

    return result;
}






int prod(int x, int depth) {
    int res = 1;
    for (int i = 0; i < depth; i++) {
        res *= x;
        x = x - 1;
    }
    return res;
}

int poly(int x, int deg, int derr) {
    if (derr == 0) {
        return pow(x, deg);
    } else {
        int coeff = prod(deg, derr);
        return coeff * pow(x, abs(deg - derr));
    }
}

vector<vector<double>> genH(int n) {
    vector<vector<double>> H(n, vector<double>(n, 0.0)); // Define an n x n matrix filled with zeros
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            H[i][j] = poly(i / 2, j, i % 2);
        }
    }
    cout << "Matrix H:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << H[i][j] << " ";
        }
        cout << endl;
    }
    vector<vector<double>> Hinv = gaussianElimination(H);
    return Hinv;
}

vector<double> getPoly(int index, const vector<vector<double>>& M) {
    vector<double> poly;
    int n = M.size();
    for (int i = 0; i < n; i++) {
        poly.push_back(M[i][index]);
    }
    return poly;
}
