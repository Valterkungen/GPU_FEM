#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib> 
#include <cmath>
#include <cuda_runtime.h>


__global__ void assembleStiffnessMatrix(int numElements, int numNodes, float h, const float* localMass, const float *localStiffness, const int *connectivity, float *globalMass, float *globalStiffness) {
    
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (elemIdx < numElements) {

        int node1 = connectivity[2 * elemIdx];     // Node i
        int node2 = connectivity[2 * elemIdx + 1]; // Node j

        // Atomic add local stiffness to global stiffness
        atomicAdd(&globalStiffness[node1 * numNodes + node1], (1.f/h) * localStiffness[0]); 
        atomicAdd(&globalStiffness[node1 * numNodes + node2], (1.f/h) * localStiffness[1]); 
        atomicAdd(&globalStiffness[node2 * numNodes + node1], (1.f/h) * localStiffness[2]); 
        atomicAdd(&globalStiffness[node2 * numNodes + node2], (1.f/h) * localStiffness[3]); 

        // Atomic add local stiffness to global mass
        atomicAdd(&globalMass[node1 * numNodes + node1], (1.f/h) * localMass[0]); 
        atomicAdd(&globalMass[node1 * numNodes + node2], (1.f/h) * localMass[1]); 
        atomicAdd(&globalMass[node2 * numNodes + node1], (1.f/h) * localMass[2]); 
        atomicAdd(&globalMass[node2 * numNodes + node2], (1.f/h) * localMass[3]); 
    }
}

__global__ void setupConnectivity(int numElements, int *connectivity) {
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < numElements) {
        int node1 = idx;               // current node
        int node2 = (idx + 1) % numElements; // next node, wraps around

        // Each element connects two nodes
        connectivity[2 * idx] = node1;
        connectivity[2 * idx + 1] = node2;
    }
}

__global__ void setMatrices(const int matrix_size, float *matrix){
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < matrix_size){
        matrix[idx] = 0.f;
    }
}


int main(int argc, char **argv) {
    
    int numNodes = 5;
    int numElements = numNodes;	
    float h = (1.f / numNodes);

    float localMass[2][2] = {{1.0f/3, 1.0f/6}, {1.0f/6, 1.0f/3}};
    float localStiffness[2][2] = {{1, -1}, {-1, 1}};
    float *localMass_gpu, *localStiffness_gpu;
    float *mass_matrix, *stiff_matrix;
    float *mass_matrix_gpu, *stiff_matrix_gpu;
    int *connectivity_gpu;

    // Allocation size for grids and blocks
    dim3 block_shape(256);
    dim3 grid_shape((numElements + block_shape.x - 1) / block_shape.x);

    // Allocate memory on GPU
    cudaMalloc(&mass_matrix_gpu, numNodes * numNodes * sizeof(float));
    cudaMalloc(&stiff_matrix_gpu, numNodes * numNodes * sizeof(float));
    cudaMalloc(&connectivity_gpu, 2 * numElements * sizeof(int));
    cudaMalloc(&localMass_gpu, 4 * sizeof(float));
    cudaMalloc(&localStiffness_gpu, 4 * sizeof(float));

    // Allocate memory on host
    mass_matrix = (float *)malloc(numNodes * numNodes * sizeof(float));
    stiff_matrix = (float *)malloc(numNodes * numNodes * sizeof(float));

    // Initialize matrices on GPU
    setMatrices<<<grid_shape, block_shape>>>(numNodes * numNodes, mass_matrix_gpu);
    setMatrices<<<grid_shape, block_shape>>>(numNodes * numNodes, stiff_matrix_gpu);
    setupConnectivity<<<grid_shape, block_shape>>>(numElements, connectivity_gpu);

    // Copy local matrices to GPU
    cudaMemcpy(localMass_gpu, localMass, 4 * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(localStiffness_gpu, localStiffness, 4 * sizeof(float), cudaMemcpyHostToDevice);

    // Run kernel to assemble matrices
    assembleStiffnessMatrix<<<grid_shape, block_shape>>>(numElements, numNodes, h,  localMass_gpu, localStiffness_gpu, connectivity_gpu, mass_matrix_gpu, stiff_matrix_gpu);

    // Copy results back to host
    cudaMemcpy(mass_matrix, mass_matrix_gpu, numNodes * numNodes * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(stiff_matrix, stiff_matrix_gpu, numNodes * numNodes * sizeof(float), cudaMemcpyDeviceToHost);

    // Print stiffness matrix
    std::cout << "Stiffness Matrix:" << std::endl;
    for (int i = 0; i < numNodes; i++) {
        for (int j = 0; j < numNodes; j++) {
            printf("%.2f ", stiff_matrix[i * numNodes + j]);
        }
        printf("\n");
    }

    // Print mass matrix
    std::cout << "Mass Matrix:" << std::endl;
    for (int i = 0; i < numNodes; i++) {
	for (int j = 0; j < numNodes; j++) {
            printf("%.2f ", mass_matrix[i * numNodes + j]);
	}
	    printf("\n");
    }

    // Free all allocated memory
    cudaFree(mass_matrix_gpu);
    cudaFree(stiff_matrix_gpu);
    cudaFree(connectivity_gpu);
    cudaFree(localMass_gpu);
    cudaFree(localStiffness_gpu);
    free(mass_matrix);
    free(stiff_matrix);
}
