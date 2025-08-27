#ifndef SOLVER_LU_H
#define SOLVER_LU_H

#include "linalg.cuh"
#include "interfaces.cuh"

template <typename T> class SolverCG: public Solver {
    private:
        Vector<T> &CPUsolve(Vector<T> &rhs);
        Vector<T> &GPUsolve(Vector<T> &rhs);
    public:
        SolverCG(FEMDiscretization &discretization);
        Vector<T> &solve(Vector<T> rhs);
}; 

#endif
