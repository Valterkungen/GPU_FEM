#ifndef SOLVER_LU_H
#define SOLVER_LU_H

#include "linalg.cuh"
#include "interfaces.cuh"

template <typename T> class SolverLU: public Solver {
    private:
        DenseMatrix<T> L, U;
        void CPUdecomp();
        __global__ void LUdecomp(int n, float *L, float *U, float *A);
        void GPUdecomp();
        Vector<T> &CPUsolve(Vector<T> &rhs);
        Vector<T> &GPUsolve(Vector<T> &rhs);
    public:
        SolverLU(FEMDiscretization &discretization);
        Vector<T> &solve(Vector<T> rhs);
}; 

#endif
