#ifndef RK4_H
#define RK4_H

#include "linalg.cuh"
#include "interfaces.cuh"

template <typename T> Vector<T> rk4(
    int steps,
    Vector<T> &x, 
    Vector<T> &v, 
    const FEMDiscretization &discretization, 
    const Solver solver,
    T dt
);

#endif