#include "solver_cg.cuh"
#include "main.cuh"

Vector<T> &SolverDirect::CPUsolve(Vector<T> &rhs) {
    // TODO: Add CPU implementation

    return Vector<T>();
}

__global__ Vector<T> &SolverDirect::GPUsolve(Vector<T> &rhs) {
    // TODO: Add GPU implementation
    
    return Vector<T>();
}

SolverDirect::SolverDirect(FEMDiscretization &discretization): discretization(discretization) {
    this->device = discretization.device;

    switch(this->device) {
        case CPU:
            this->CPUdecomp<T>(this->L, this->U);
            break;
        case GPU:
            this->GPUdecomp<T><<<numberOfBlocks, threadsPerBlock>>>(this->L, this->U);
            break;
    }
}

Vector<T> &SolverDirect::solve(Vector<T> rhs) {
    switch(this->device) {
        case CPU:
            this->CPUsolve<T>(rhs);
            break;
        case GPU:
            this->GPUsolve<T><<<numberOfBlocks, threadsPerBlock>>>(rhs);
            break;
    }
}