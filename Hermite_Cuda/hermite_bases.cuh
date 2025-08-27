#include <iostream>
#include <armadillo>
#include <cuda_runtime.h>

using namespace std;
using namespace arma;

__device__ int dev_prod(int x, int depth) {
    int res = 1;
    for (int i = 0; i < depth; i++) {
        res *= x;
        x = x - 1;
    }
    return res;
}

__device__ int dev_poly(int x, int deg, int derr) {
    if (derr == 0) {
        return pow(x, deg);
    } else {
        int coeff = dev_prod(deg, derr);
        return coeff * pow(x, abs(deg - derr));
    }
}

__global__ void genA_kernel(mat::fixed<100,100>* A, int n) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid < n * n) {
        int i = tid / n;
        int j = tid % n;
        (*A)(i, j) = dev_poly(i / 2, j, i % 2);
    }
}

mat genA(int n) {
    mat A(n, n, fill::zeros);
    mat::fixed<100,100>* dev_A;
    cudaMalloc((void**)&dev_A, sizeof(mat::fixed<100,100>));
    dim3 block_shape(256);
    dim3 grid_shape((n * n + block_shape.x - 1) / block_shape.x);
    genA_kernel<<<grid_shape, block_shape>>>(dev_A, n);
    cudaMemcpy(&A(0,0), dev_A, sizeof(mat::fixed<100,100>), cudaMemcpyDeviceToHost);
    cudaFree(dev_A);
    return inv(A);
}

vec getPoly(int index, const mat& M) {
    return M.col(index);
}


