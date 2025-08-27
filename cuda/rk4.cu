#include "rk4.cuh"

template <typename T> void rk4(int steps, Vector<T> &x, Vector<T> &v, const FEMDiscretization &discretization, const Solver &solver, T dt) {
    /**
    * RK4 time integrator. 
    *
    * @param x Refrence to spatial vector. 
    * @param v Reference to derivative of ``x``. 
    * @param discretization Refrence to discretization. 
    * @param solver Refrence to solver. 
    */
    
    for(int i = 0; i < steps; i++) {
        auto k1_x = v;
        auto k1_v = solver.solve(discretization.F - discretization.A*x);

        auto k2_x = v + .5*dt*dt*k1_v;
        auto k2_v = solver.solve(discretization.F - discretization.A*(x + .5*dt*k1_x));

        auto k3_x = v + 0.5*dt*k2_v;
        auto k3_v = solver.solve(discretization.F - discretization.A*(x + .5*dt*k2_x));

        auto k4_x = v + dt*k3_v;
        auto k4_v = solver.solve(discretization.F - discretization.A*(x + dt*k3_x));

        x += dt*(k1_x + 2*k2_x + 2*k3_x + k4_x) / 6;
        v += dt*(k1_v + 2*k2_v + 2*k3_v + k4_v) / 6;
    }
}
