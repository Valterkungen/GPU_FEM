#include <iostream>
#include <vector>
#include <armadillo>
#include "hermite_bases.cuh"
#include "hermite_polynomials.cuh"
#include <cuda_runtime.h>

using namespace std;
using namespace arma;

__global__ void GenMatrices_kernel(int dim, double h, double* Msub_dev, double* Asub_dev) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid < dim * dim) {
        int i = tid / dim;
        int j = tid % dim;
        mat H = genA(dim);
        vec coeffs1 = getPoly(i, H);
        Polynomial p1(coeffs1);
        Polynomial p1prim = p1.derive();
        vec coeffs2 = getPoly(j, H);
        Polynomial p2(coeffs2);
        Polynomial p2prim = p2.derive();
        Msub_dev[tid] = (p1 * p2).integrate(0.0, 1.0);
        Asub_dev[tid] = (p1prim * p2prim).integrate(0.0, 1.0);
    }
}

pair<mat,mat> GenMatrices(int dim, const vec& xvec, const vec& f, const string& BC) {
    int N = xvec.size();
    int nel = N - 1;
    int ndof = 2 * N; //might be dim dependent 2 -> (dim/2)
    double h = 1.0 / nel;

    mat Msub(dim, dim, fill::zeros);
    mat Asub(dim, dim, fill::zeros);

    double* Msub_dev;
    double* Asub_dev;
    cudaMalloc((void**)&Msub_dev, dim * dim * sizeof(double));
    cudaMalloc((void**)&Asub_dev, dim * dim * sizeof(double));
    dim3 block_shape(256);
    dim3 grid_shape((dim * dim + block_shape.x - 1) / block_shape.x);
    GenMatrices_kernel<<<grid_shape, block_shape>>>(dim, h, Msub_dev, Asub_dev);
    cudaMemcpy(Msub.memptr(), Msub_dev, dim * dim * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(Asub.memptr(), Asub_dev, dim * dim * sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(Msub_dev);
    cudaFree(Asub_dev);

    mat M(ndof, ndof, fill::zeros);
    mat A(ndof, ndof, fill::zeros);
    vector<double> F(ndof, 0);

    // Assembly loop
    for (int k = 0; k < nel; k++) {
        int kl = 2 * k;
        int ku = 2 * k + dim;

        // Loop over the rows and columns of the local matrices
        for (int i = kl; i < ku; i++) {
            for (int j = kl; j < ku; j++) {
                // Update the corresponding elements of the global matrices              
                M(i, j) += Msub(i - kl, j - kl)*(1.0 / h);
                A(i, j) += Asub(i - kl, j - kl)*(1.0 / h);
            }
        }
    }
    if (BC == "dirichlet") {
        // Implement the Dirichlet boundary condition handling here
    }
    
    else if (BC == "periodic"){
        int st = ndof - dim/2;        
        for (int i = 0; i < dim; i++){
            for (int j = 0; j < dim; j++){
                M((st+i)%(ndof),(st+j)%(ndof)) += Msub(i,j)*(1.0/h);                            
                A((st+i)%(ndof),(st+j)%(ndof)) += Asub(i,j)*(1.0/h);  
            }
        }
    }
    
    return make_pair(M,A);
}

int main() {
    // Sample input vectors
    vec xvec = linspace<vec>(0, 1, 5);
    // Initialize f with 100 elements of value 0
    vec f = zeros<vec>(4);
    string BC = "periodic"; 
    int dim = 4;
    // Call the GenMatrices function
    pair<mat, mat> matrices = GenMatrices(dim, xvec, f, BC);
    mat M = matrices.first;
    mat A = matrices.second;
    cout << "Mass Matrix" << endl;
    cout << M << endl;
    cout << "Stiffness Matrix" << endl;
    cout << A << endl;
    
    return 0;
}
