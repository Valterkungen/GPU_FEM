#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib> 


__global__ void CSR_Assembly(const int size, const float h , float *row_offsets, float *col_indices, float *Stiffness, float *Mass) {
    
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (idx < size) {
        int row_idx = 3 * idx;  // Varje rad har 3 element
        int nextIdx = (idx + 1) % size;  // Hantera Periodiska randvillkor
        int prevIdx = (idx - 1 + size) % size;  

        row_offsets[idx] = row_idx;

        // Kolumn index för 3 närliggande element
        col_indices[row_idx] = prevIdx;   
        col_indices[row_idx + 1] = idx;   
        col_indices[row_idx + 2] = nextIdx; 

        //Sätt Stiffnes matrix
        Stiffness[row_idx] = -1/h;
        Stiffness[row_idx + 1] = 2/h;
        Stiffness[row_idx + 2] = -1/h;

        //Sätt Mass Matrix
        Mass[row_idx] = h/6;
        Mass[row_idx + 1] = 4*h/6;
        Mass[row_idx + 2] = h/6;
    }

    // Ta hand om sista i rows
    if (idx == size) {
        row_offsets[size] = 3 * size;  
    }
}


int main(int argc, char **argv){

    float *row_offset, *row_offset_gpu;
    float *col_indicies, *col_indicies_gpu;
    float *Stiffness, *Stiffness_gpu;
    float *Mass, *Mass_gpu;

    int row_size = 5;
    int val_size = 3 * row_size;
    float h = 1;
    dim3 block_shape = dim3(32, 32);

    dim3 grid_shape = dim3(max(1.0, ceil((float) row_size / (float) block_shape.x)), max(1.0, ceil((float) row_size / (float) block_shape.y)));

    cudaMalloc(&row_offset_gpu, (row_size + 1) * sizeof(float));
    cudaMalloc(&col_indicies_gpu, val_size * sizeof(float));
    cudaMalloc(&Stiffness_gpu, val_size * sizeof(float));
    cudaMalloc(&Mass_gpu, val_size * sizeof(float));

    row_offset= (float *) malloc((row_size + 1) * sizeof(float));
    col_indicies = (float *) malloc(val_size * sizeof(float));
    Stiffness = (float *) malloc(val_size * sizeof(float));
    Mass = (float *) malloc(val_size * sizeof(float));
        
    CSR_Assembly<<<grid_shape, block_shape>>>(row_size, h, row_offset_gpu, col_indicies_gpu, Stiffness_gpu, Mass_gpu);

    cudaMemcpy(row_offset, row_offset_gpu, (row_size + 1) * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(col_indicies, col_indicies_gpu, val_size * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(Stiffness, Stiffness_gpu, val_size * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(Mass, Mass_gpu, val_size * sizeof(float), cudaMemcpyDeviceToHost);
    

    std::cout << "Row offsets" << std::endl;
    for (int i = 0; i < row_size + 1; i++){
        printf("%.2f\n", row_offset[i]);
        }

    std::cout << "Column indicies" << std::endl;

    for (int i = 0; i < val_size; i++){
        printf("%.2f\n", col_indicies[i]);
        }

    std::cout << "Stiffness" << std::endl;

    for (int i = 0; i < val_size; i++){
        printf("%.2f\n", Stiffness[i]);
        }

    std::cout << "Mass" << std::endl;

    for (int i = 0; i < val_size; i++){
        printf("%.2f\n", Mass[i]);
        }

    free(row_offset);
    free(col_indicies);
    free(Stiffness);
    free(Mass);

    cudaFree(row_offset_gpu);
    cudaFree(col_indicies_gpu);
    cudaFree(Stiffness_gpu);
    cudaFree(Mass_gpu);

    return 0;
}