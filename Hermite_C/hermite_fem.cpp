#include <iostream>
#include <vector>
#include "hermite_bases.h"
#include "hermite_polynomials.h"

using namespace std;

pair<vector<vector<double>>, vector<vector<double>>> GenMatrices(int dim, const vector<double>& xvec, const vector<double>& f, const string& BC) {
    int N = xvec.size();
    int nel = N - 1;
    int ndof = 2 * N; // might be dim dependent 2 -> (dim/2)
    double h = 1.0 / nel;

    vector<vector<double>> Msub(dim, vector<double>(dim, 0));
    vector<vector<double>> Asub(dim, vector<double>(dim, 0));

    //Gen H

    vector<vector<double>> Hinv = genH(dim);
    
    for (int i = 0; i < dim; i++) {
        cout << "Before getPoly" << endl;
        vector<double> coeffs1 = getPoly(i, Hinv); // Corrected argument
        cout << "After getPoly" << endl;
        Polynomial p1(coeffs1);
        Polynomial p1prim = p1.derive();
        for (int j = 0; j < dim; j++) {
            vector<double> coeffs2 = getPoly(j, Hinv); // Corrected argument
            Polynomial p2(coeffs2);
            Polynomial p2prim = p2.derive();
            Msub[i][j] = (p1 * p2).integrate(0.0, 1.0);
            Asub[i][j] = (p1prim * p2prim).integrate(0.0, 1.0);
        }
    }

    cout << "Msub" << endl;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            cout << Msub[i][j] << " ";
        }
        cout << endl;
    }
    cout << "Asub" << endl;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            cout << Asub[i][j] << " ";
        }
        cout << endl;
    }

    vector<vector<double>> M(ndof, vector<double>(ndof, 0));
    vector<vector<double>> A(ndof, vector<double>(ndof, 0));
    vector<double> F(ndof, 0);

    // Assembly loop
    for (int k = 0; k < nel; k++) {
        int kl = 2 * k;
        int ku = 2 * k + dim;

        // Loop over the rows and columns of the local matrices
        for (int i = kl; i < ku; i++) {
            for (int j = kl; j < ku; j++) {
                // Update the corresponding elements of the global matrices
                M[i][j] += Msub[i - kl][j - kl] * (1.0 / h);
                A[i][j] += Asub[i - kl][j - kl] * (1.0 / h);
            }
        }
    }
    if (BC == "dirichlet") {
        // Implement the Dirichlet boundary condition handling here
    }

    else if (BC == "periodic") {
        int st = ndof - dim / 2;
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                M[(st + i) % (ndof)][(st + j) % (ndof)] += Msub[i][j] * (1.0 / h);
                A[(st + i) % (ndof)][(st + j) % (ndof)] += Asub[i][j] * (1.0 / h);
            }
        }
    }

    return make_pair(M, A);
}

int main() {
    // Sample input vectors
    vector<double> xvec = {0, 0.25, 0.5, 0.75, 1};
    // Initialize f with 100 elements of value 0
    vector<double> f(4, 0);
    string BC = "periodic";
    int dim = 4;
    // Call the GenMatrices function
    pair<vector<vector<double>>, vector<vector<double>>> matrices = GenMatrices(dim, xvec, f, BC);
    vector<vector<double>> M = matrices.first;
    vector<vector<double>> A = matrices.second;

    cout << "Mass Matrix" << endl;
    for (const auto& row : M) {
        for (double val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
    cout << "Stiffness Matrix" << endl;
    for (const auto& row : A) {
        for (double val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
    return 0;
}

