#ifndef FEM_HAT_H
#define FEM_HAT_H

#include <functional>
#include "interfaces.cuh"

template <typename T> class FEMHatDiscretization: public FEMDiscretization<T> {
    private:
        void CPUdiscretization();
        void GPUdiscretization();
    public:
        FEMHatDiscretization(
            Device device, 
            Vector<T> x, 
            std::function<T(T)> f
        );
};

#endif
