#ifndef SOLVER_LU_H
#define SOLVER_LU_H

#include "linalg.cuh"
#include "interfaces.cuh"

template <typename T> class SolverDirect: public Solver {
    private:
        Vector<T> &CPUsolve(Vector<T> &rhs);
        __global__ Vector<T> &GPUsolve(Vector<T> &rhs);
    public:
        SolverDirect(FEMDiscretization &discretization);
        Vector<T> &solve(Vector<T> rhs);
}; 

#endif
