#ifndef INTERFACES_H
#define INTERFACES_H

#include "device.cuh"
#include "linalg.cuh"

template <typename T> class FEMDiscretization {
    protected:
        Device device;
    public:
        SparseMatrix<T> M, A;
        Vector<T> F;
        FEMDiscretization();
        typedef T type;
};

template <typename T> class Solver {
    protected:
        Device device;
        FEMDiscretization<T> &discretization;
    public:
        Solver();
        virtual Vector<T> &solve(SparseMatrix<T> rhs);
};

#endif