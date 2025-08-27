#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib> 
#include <cmath>
#include <cuda_runtime.h>


__global__ void matrix_vector_product(const float *M, const float *v, int vec_size, float *y) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < vec_size) {
        float temp = 0.0;
        for (int i = 0; i < vec_size; i++) {
            temp += M[idx * vec_size + i] * v[i];
        }
        y[idx] = temp;
    }
}

__device__ float vector_dot_product(const float *v, const float *u, int n, float y) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x; 
    if (idx < n ) {
        y += v[idx] * u[idx]
    }
}


__global__ void LUfactorisation(int n, float *L, float *U, float *A){
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if(i < n){
        L[i*n + i] = 1.0; // set diagonal of L to 1
        for(int j = 0; j<i; j++ ){
            float sum = 0.0;
            for(int k = 0; k < i, k++){
                sum += L[i * n + k] * U[k * n + j]; 
            }
            U[i * n + j] = A[i * n + j] - sum; //Computes U[i,j]
        }
        for(int j = i + 1; j < n; j++ ){
            float sum = 0.0;
            for(int k = 0; k < i; k++){
                sum += L[j*n + k] * U[k*n +1];
            }
            L[j*n + i] = (A[j*n + i] - sum) / U[i*n + i]; //Computes L[i,j]
        }
    }
}

float* forward_substitution(int n, float *L, float *b){
    float *y = new float[n]; 
    for(int i = 0; i < n; i++){
        float sum = 0.0;
        for(int j = 0; j < i; j++)
            sum += L[i*n + j] * y[j];
        y[i] = (b[i] - sum) / L[i*n + i]
    return *y
    }
}

float* backward_substitution(int n, float *U, float *y){
        float *x = new float[n];
        for(int i = n-1; i < -1; i--){
            float sum = 0.0; 
            for(int j = i+1; j < i; k++) {
                sum += U[n*(n-1) - i*n + k] * y[n-k];
            x[i] 
            }
        }
}

int main(){

    int n = 3;
    // Allocate and initialize memory on host
    float A_cpu[9] = {4, -1, 1, -1, 4, -2, 1, -2, 4}; // Example matrix
    float L_cpu[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    float U_cpu[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    float b_cpu[3] = {12, -1, 5}; // Example vectors
    float x0_cpu[3] = {1, 2, 1};
    float x[3] = {0, 0, 0};
    dim3 block_shape(256);
    dim3 grid_shape(())

    // Allocate memory for L, U and A matrices 
    float *L_gpu, *U_gpu, *A_gpu;
    cudaMalloc(&L_gpu, n * n * sizeof(float));
    cudaMalloc(&U_gpu, n * n * sizeof(float));
    cudaMalloc(&A_gpu, n * n * sizeof(float));

    
    // Copy matrice A to GPU
    cudaMemcpy(L_gpu, L_cpu, n * n * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(U_gpu, U_cpu, n * n * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(A_gpu, A_cpu, n * n * sizeof(float), cudaMemcpyHostToDevice);


    //Run LU factorisation program 

    LUfactorisation<<<grid_shape, block_shape>>>(n, L_gpu, U_gpu, A_gpu);

    // Copy matrice L and U to CPU 
    cudaMemcpy(L_cpu, L_gpu, n * n * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(U_cpu, U_gpu, n * n * sizeof(float), cudaMemcpyDeviceToHost);

    //     



    //allocate memory for initial guess and result
    float *x0_gpu, *xnew_gpu;
    cudaMalloc(&x0_gpu, n * sizeof(float));
    cudaMalloc(&xnew_gpu, n * sizeof(float));

    // Copy initial guess and 

}