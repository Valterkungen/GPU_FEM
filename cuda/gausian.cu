#include <iostream>
#include <vector>
#include <cmath>
#include <cuda_runtime.h>

using namespace std;

__global__ void partial_pivot_kernel(double* A, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    // Pivot Selection
    int pivot = i;
    double maxElement = fabs(A[i * n + i]);
    for (int row = i + 1; row < n; row++) {
        double element = fabs(A[row * n + i]);
        if (element > maxElement) {
            maxElement = element;
            pivot = row;
        }
    }

    if (pivot != i) {
        // Swap rows
        for (int j = 0; j <= n; j++) {
            double tmp = A[i * (n+1) + j];
            A[i * (n+1) + j] = A[pivot * (n+1) + j];
            A[pivot * (n+1) + j] = tmp;
        }
    }

    __syncthreads();

    // Row Reduction
    for (int row = i + 1; row < n; row++) {
        double factor = A[row * (n+1) + i] / A[i * (n+1) + i];
        for (int col = i; col <= n; col++) {
            A[row * (n+1) + col] -= A[i * (n+1) + col] * factor;
        }
    }
}

__global__ void back_substitution_kernel(double* A, double* x, int n) {
    int i = n - 1 - (blockIdx.x * blockDim.x + threadIdx.x);
    if (i < 0) return;

    double sum = 0.0;
    for (int j = i + 1; j < n; j++) {
        sum += A[i * (n+1) + j] * x[j];
    }
    x[i] = (A[i * (n+1) + n] - sum) / A[i * (n+1) + i];
}

int main() {
    const int n = 3;
    double h_A[n * (n+1)] = {
        3, 2, -4, 3,
        2, 3, 3, 15,
        5, -3, 1, 14
    };

    double* d_A;
    cudaMalloc(&d_A, n * (n+1) * sizeof(double));
    cudaMemcpy(d_A, h_A, n * (n+1) * sizeof(double), cudaMemcpyHostToDevice);

    double* d_x;
    cudaMalloc(&d_x, n * sizeof(double));

    int blockSize = 256;
    int numBlocks = (n + blockSize - 1) / blockSize;

    // Launch kernel to perform partial pivoting and row reduction
    partial_pivot_kernel<<<numBlocks, blockSize>>>(d_A, n);
    cudaDeviceSynchronize();

    // Launch kernel to perform back substitution
    back_substitution_kernel<<<numBlocks, blockSize>>>(d_A, d_x, n);
    cudaDeviceSynchronize();

    // Copy result back to host
    double h_x[n];
    cudaMemcpy(h_x, d_x, n * sizeof(double), cudaMemcpyDeviceToHost);

    cout << "Solution:" << endl;
    for (int i = 0; i < n; i++) {
        cout << "x[" << i << "] = " << h_x[i] << endl;
    }

    cudaFree(d_A);
    cudaFree(d_x);
    return 0;
}
