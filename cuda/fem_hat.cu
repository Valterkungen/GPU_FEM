#include "fem_hat.cuh"
#include "main.cuh"

static void CSR_Assembly_CPU(const int size, const float h, float *row_offsets, float *col_indices, float *stiffness, float *mass) {
    
    float h_inv = 1.f / h;

    for (int idx = 0; idx <= size; idx++) { 
        if (idx < size) {
            int row_idx = 3 * idx;  
            int nextIdx = (idx + 1) % size; 
            int prevIdx = (idx - 1 + size) % size;

            row_offsets[idx] = row_idx;
            
            col_indices[row_idx] = prevIdx;
            col_indices[row_idx + 1] = idx;
            col_indices[row_idx + 2] = nextIdx;

            stiffness[row_idx] = -1 * h_inv;
            stiffness[row_idx + 1] = 2 * h_inv;
            stiffness[row_idx + 2] = -1 * h_inv;

            mass[row_idx] = h / 6;
            mass[row_idx + 1] = 4 * h / 6;
            mass[row_idx + 2] = h / 6;
        }

        if (idx == size) {
            row_offsets[size] = 3 * size;
        }
    }
}


__global__ static void CSR_Assembly_GPU(const int size, const float h , float *row_offsets, float *col_indices, float *stiffness, float *mass) {
    
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    float h_inv = 1.f / h;
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
        stiffness[row_idx] = -1 * h_inv;
        stiffness[row_idx + 1] = 2 * h_inv;
        stiffness[row_idx + 2] = -1 * h_inv;

        //Sätt mass Matrix
        mass[row_idx] = h/6;
        mass[row_idx + 1] = 4*h/6;
        mass[row_idx + 2] = h/6;
    }

    // Ta hand om sista i rows
    if (idx == size) {
        row_offsets[size] = 3 * size;  
    }
}

void FEMHatDiscretization::CPUdiscretization() {
    // Add CPU implementation here
    float *row_offset_cpu;
    float *col_indicies_cpu, *col_indicies_cpu_copy;
    float *stiffness_cpu;
    float *mass_cpu;

    int num_elements = 3;
    int row_size = x.getLenght();
    int col_size = row_size;
    int val_size = num_elements * row_size;
    float h = 1.f/row_size;

    row_offset_cpu = (float *) malloc((row_size + 1) * sizeof(float));
    col_indicies_cpu = (float *) malloc(val_size * sizeof(float));
    col_indicies_cpu_copy = (float *) malloc(val_size * sizeof(float));
    stiffness_cpu = (float *) malloc(val_size * sizeof(float));
    mass_cpu = (float *) malloc(val_size * sizeof(float));

        
    CSR_Assembly_CPU(row_size, h, row_offset_cpu, col_indicies_cpu, stiffness_cpu, mass_cpu);
    memcpy(col_indicies_cpu_copy, col_indicies_cpu, val_size * sizeof(float));

    this->M = SparseMatrix<T>(Device::CPU, row_size, col_size, num_elements, stiffness_cpu, col_indicies_cpu);
    this->A = SparseMatrix<T>(Device::CPU, row_size, col_size, num_elements, mass_cpu, col_indicies_cpu_copy);
    //this->F = ...;
}

void FEMHatDiscretization::GPUdiscretization() {
    // Add GPU implementation here
    float *row_offset_gpu;
    float *col_indicies_gpu, *col_indicies_gpu_copy;
    float *stiffness_gpu;
    float *mass_gpu;

    int num_elements = 3;
    int row_size = x.getLength()
    int val_size = 3 * row_size;
    float h = 1.f / row_size;

    cudaMalloc(&row_offset_gpu, (row_size + 1) * sizeof(float));
    cudaMalloc(&col_indicies_gpu, val_size * sizeof(float));
    cudaMalloc(&col_indicies_gpu_copy, val_size * sizeof(float));
    cudaMalloc(&stiffness_gpu, val_size * sizeof(float));
    cudaMalloc(&mass_gpu, val_size * sizeof(float));

        
    CSR_Assembly_GPU<<<numberOfBlocks, threadsPerBlock>>>(row_size, h, row_offset_gpu, col_indicies_gpu, stiffness_gpu, mass_gpu);
    cudaMemcpy(col_indicies_gpu_copy, col_indicies_gpu, val_size * sizeof(float), cudaMemcpyDeviceToDevice);

    this->M = SparseMatrix<T>(Device::GPU, row_size, col_size, num_elements, stiffness_gpu, col_indicies_gpu);
    this->A = SparseMatrix<T>(Device::GPU, row_size, col_size, num_elements, mass_gpu, col_indicies_gpu_copy);
    //this->F = ...;
    
}

FEMHatDiscretization::FEMHatDiscretization(Device device, Vector<T> x, function<T(T)> f): device(device) {
    /**
    * Constructor for FEM Hat Discretization. 
    *
    * @param device Device, either ``cpu`` or ``gpu``. 
    * @param x Spatial vector. 
    * @param f Forcing function. 
    */

    // Struntar i Load tillsvidare, se den som noll 
    this->F = Vector<T>(x.getLength());
    F.setZeros();

    switch(this->device) {
        case CPU:
            this->CPUdiscretization();
            break;
        case GPU:
            this->GPUdiscretization();
            break;
    }
}