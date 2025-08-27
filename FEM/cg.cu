#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <cuda_runtime.h>


__global__ void vector_vector_product(const float *v1, const float *v2, int vec_size, float *result){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < vec_size) {
        atomicAdd(result, v1[idx] * v2[idx]);
    }

}
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

__global__ void vector_addition(const float *x, const float *y, int vec_size, float alpha, float *result) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < vec_size) {
        result[idx] = x[idx] + alpha * y[idx];
    }
}

__global__ void vector_subtraction(const float *x, const float *y, int vec_size, float alpha, float *result) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < vec_size) {
        result[idx] = x[idx] - alpha * y[idx];
    }
}

__global__ void scalar_vector_product(float scalar, const float *v, int vec_size, float *result) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < vec_size) {
        result[idx] = scalar * v[idx];
    }
}

__global__ void vector_difference(const float *v1, const float *v2, int vec_size, float *result) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < vec_size) {
        result[idx] = v1[idx] - v2[idx];
    }
}

__global__ void vector_norm_squared(const float *v, int vec_size, float *result) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < vec_size) {
        atomicAdd(result, v[idx] * v[idx]);
    }
}

void conjugate_gradient(const float *A_cpu, const float *b_cpu, float *x_cpu, int vec_size, float *x0 = NULL, float tol = 1e-6) {
    float *A_gpu, *b_gpu, *x_gpu;
    float *r, *p, *Ax, *w, *diff, *x_new;
    float *norm_x_new_squared, *norm_diff_squared;
    float h_norm_x_new_squared, h_norm_diff_squared, err, rho, rho_prev;
    float *rho_gpu; 
    dim3 block_shape(256);
    dim3 grid_shape((vec_size + block_shape.x - 1) / block_shape.x);

    // Allocating Memory on GPU
    cudaMalloc(&A_gpu, vec_size * vec_size * sizeof(float));
    cudaMalloc(&b_gpu, vec_size * sizeof(float));
    cudaMalloc(&x_gpu, vec_size * sizeof(float));
    cudaMalloc(&r, vec_size * sizeof(float));
    cudaMalloc(&p, vec_size * sizeof(float));
    cudaMalloc(&Ax, vec_size * sizeof(float));
    cudaMalloc(&w, vec_size * sizeof(float));
    cudaMalloc(&diff, vec_size * sizeof(float));
    cudaMalloc(&x_new, vec_size * sizeof(float));
    cudaMalloc(&norm_diff_squared, sizeof(float));
    cudaMalloc(&norm_x_new_squared, sizeof(float));
    cudaMalloc(&norm_diff_squared, sizeof(float));
    cudaMalloc(&rho_gpu, sizeof(float));
    cudaMemset(rho_gpu, 0, sizeof(float));


    cudaMemcpy(A_gpu, A_cpu, vec_size * vec_size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(b_gpu, b_cpu, vec_size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(x_gpu, x_cpu, vec_size * vec_size * sizeof(float), cudaMemcpyHostToDevice);

    // Initialize vectors
    cudaMemcpy(x_new, x0, vec_size * sizeof(float), cudaMemcpyHostToDevice); // Assume initial guess x=x0

    //Initalize Ax, residual, copy, rho
    matrix_vector_product<<<grid_shape, block_shape>>>(A_gpu, x_new, vec_size, Ax);
    vector_subtraction<<<grid_shape, block_shape>>>(b_gpu, Ax, vec_size, 1.0, r);
    vector_vector_product<<<grid_shape, block_shape>>>(r, r, vec_size, rho_gpu);
    cudaMemcpy(&rho, rho_gpu, sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(p, r, vec_size * sizeof(float), cudaMemcpyDeviceToDevice);
    // Set initial error
    err = 2*tol;


    //For testing
    int idx = 1;
        while (err > tol){

        cudaMemcpy(x_gpu, x_new, vec_size * sizeof(float), cudaMemcpyDeviceToDevice);
        
        // Compute Ap
        matrix_vector_product<<<grid_shape, block_shape>>>(A_gpu, p, vec_size, w);
        cudaDeviceSynchronize();
        
        // Compute dot products needed for alpha: rho_new = r^T r and dot_pAp = p^T Ap
        float dot_pAp;
        float *d_dot_pAp;
        cudaMalloc(&d_dot_pAp, sizeof(float));
        cudaMemset(d_dot_pAp, 0, sizeof(float));
        vector_vector_product<<<grid_shape, block_shape>>>(p, w, vec_size, d_dot_pAp);
        cudaMemcpy(&dot_pAp, d_dot_pAp, sizeof(float), cudaMemcpyDeviceToHost);

        if (dot_pAp == 0){
            std::cout << "Anorm_p_squared is zero" << std::endl;
            break;

        }
        
        float alpha = rho / dot_pAp;

        
        // x_new = x + alpha * p
        vector_addition<<<grid_shape, block_shape>>>(x_gpu, p, vec_size, alpha, x_new);
        cudaDeviceSynchronize();
        
    

        // r = r - alpha * Ap
        vector_subtraction<<<grid_shape, block_shape>>>(r, w, vec_size, alpha, r);
        cudaDeviceSynchronize();

        // Update Rho
        rho_prev = rho;
        cudaMemset(rho_gpu, 0, sizeof(float));
        vector_norm_squared<<<grid_shape, block_shape>>>(r, vec_size, rho_gpu);
        cudaMemcpy(&rho, rho_gpu, sizeof(float), cudaMemcpyDeviceToHost);


        // Update error
        cudaMemset(norm_diff_squared, 0, sizeof(float));
        cudaMemset(norm_x_new_squared, 0, sizeof(float));
        
        vector_difference<<<grid_shape, block_shape>>>(x_new, x_gpu, vec_size, diff);
        cudaDeviceSynchronize();

        vector_norm_squared<<<grid_shape, block_shape>>>(diff, vec_size, norm_diff_squared);
        cudaDeviceSynchronize();
        vector_norm_squared<<<grid_shape, block_shape>>>(x_new, vec_size, norm_x_new_squared);
        cudaDeviceSynchronize();
        cudaMemcpy(&h_norm_diff_squared, norm_diff_squared, sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpy(&h_norm_x_new_squared, norm_x_new_squared, sizeof(float), cudaMemcpyDeviceToHost);
        
        // p = r + rho_new/rho * p
        float Omega = rho/rho_prev;
        
        vector_addition<<<grid_shape, block_shape>>>(r, p, vec_size, Omega, p);
        cudaDeviceSynchronize();    
        
        err = sqrt(h_norm_diff_squared) / sqrt(h_norm_x_new_squared);

        // Cleanup
        cudaFree(d_dot_pAp);
        cudaFree(rho_gpu);

        idx = idx + 1;
}
    // Free all allocated memory
    cudaMemcpy(x_cpu, x_new, vec_size * sizeof(float), cudaMemcpyDeviceToDevice);
    cudaFree(r);
    cudaFree(p);
    cudaFree(Ax);
    cudaFree(w);
    cudaFree(diff);
    cudaFree(norm_x_new_squared);
    cudaFree(norm_diff_squared);
    cudaFree(x_new);
    std::cout << "Converge after " << idx << " iterations" << std::endl;

}

int main(int argc, char **argv) {

    int vec_size = 3;
    // Allocate and initialize memory on host
    float A_cpu[9] = {4, -1, 1, -1, 4, -2, 1, -2, 4}; // Example matrix
    float b_cpu[3] = {12, -1, 5}; // Example vectors
    float x0_cpu[3] = {1, 2, 1};
    float x[3] = {0, 0, 0};

    // Allocate memory for the initial guess and the result
    float *x0_gpu, *x_new_gpu;
    cudaMalloc(&x0_gpu, vec_size * sizeof(float));
    cudaMalloc(&x_new_gpu, vec_size * sizeof(float));
    
    // Copy the initial guess to the GPU
    cudaMemcpy(x0_gpu, x0_cpu, vec_size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(x_new_gpu, x, vec_size * sizeof(float), cudaMemcpyHostToDevice);
    std::cout << "Run Conjugate Gradient" << std::endl; 

    // Run Conjugate Gradient Solver
    conjugate_gradient(A_cpu, b_cpu, x_new_gpu, vec_size);

    // Copy results back to host
    float *x_new_cpu;
    x_new_cpu = (float *)malloc(vec_size * sizeof(float));
    cudaMemcpy(x_new_cpu, x_new_gpu, vec_size * sizeof(float), cudaMemcpyDeviceToHost);

    // Print results
    std::cout << "Solution vector x:" << std::endl;
    for (int i = 0; i < vec_size; i++) {
        std::cout << x_new_cpu[i] << " ";
    }

    // Free allocated memory
    cudaFree(x0_gpu);
    cudaFree(x_new_gpu);

    return 0;
}
