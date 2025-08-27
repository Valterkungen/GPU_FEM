#include "fem_hermite.cuh"
#include "main.cuh"

void FEMHermiteDiscretization::CPUdiscretization() {
    // Add CPU implementation here
    //this->M = ...;
    //this->A = ...;
    //this->F = ...;
}

__global__ void FEMHermiteDiscretization::GPUdiscretization() {
    // Add GPU implementation here
    //this->M = ...;
    //this->A = ...;
    //this->F = ...;
    
}

FEMHermiteDiscretization::FEMHermiteDiscretization(Device device, Vector<T> x, function<T(T)> f): device(device) {
    /**
    * Constructor for FEM Hermite Discretization. 
    *
    * @param device Device, either ``cpu`` or ``gpu``. 
    * @param x Spatial vector. 
    * @param f Forcing function. 
    */

    this->M = SparseMatrix<T>(x.getLength(), x.getLength());
    this->A = SparseMatrix<T>(x.getLength(), x.getLength());
    this->F = Vector<T>(x.getLength());

    switch(this->device) {
        case CPU:
            this->CPUdiscretization();
            break;
        case GPU:
            this->GPUdiscretization<<<numberOfBlocks, threadsPerBlock>>>();
            break;
    }
}