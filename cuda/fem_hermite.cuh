#ifndef FEM_HAT_H
#define FEM_HAT_H

#include "interfaces.cuh"

template <typename T> class FEMHermiteDiscretization: public FEMDiscretization<T> {
    private:
        void CPUdiscretization();
        __global__ void GPUdiscretization();
    public:
        FEMHermiteDiscretization(
            Device device, 
            Vector<T> x, 
            function<T(T)> f
        );
}; 

#endif
