#define _USE_MATH_DEFINES

#include <iostream>
#include <stdexcept>
#include <sstream>
#include <cmath>

#include "main.cuh"
#include "rk4.cuh"
#include "interfaces.cuh"
#include "linalg.cuh"
#include "fem_hat.cuh"
#include "fem_hermite.cuh"
#include "solver_lu.cuh"
#include "solver_cg.cuh"
#include "solver_direct.cuh"

using namespace std;

int main(int argc, char *argv[]) {
    // Handle input arguments
    if(argc != 10) {
        cerr << "Incorrect number of input arguemnts. " << endl;
        cout << "Usage: ./app <plot> <n steps> <grid size> <dt> <blocks> <threads> <basis function> <solver> <device>" << endl;
        cout << "            plot: Export plotable data? (yes/no)" << endl;
        cout << "         n steps: Number of time steps" << endl;
        cout << "       grid size: Number of grid points. " << endl;
        cout << "              dt: Length of time steps." << endl;
        cout << "          blocks: Number of CUDA blocks. " << endl;
        cout << "         threads: Number of threads per CUDA block. " << endl;
        cout << "  basis function: Type of basis function. (hat/hermite) " << endl;
        cout << "          solver: Type of solver. (direct/cg/lu)" << endl;
        cout << "          device: Type of device. (cpu/gpu)" << endl;
        return -1; 
    }

    bool plot;
    int nsteps, gridsize;
    double dt;
    FEMDiscretization discretization;
    Solver solver;
    Device device;

    Vector<double> x;
    Vector<double> v;

    try {
        // Arg 1: Plot
        plot = string(argv[1]).compare("yes") == 0;

        // Arg 2: N steps
        istringstream s2(argv[2]);
        s2 >> nsteps;

        // Arg 3: Gridsize
        istringstream s3(argv[3]);
        s3 >> gridsize;
        x = Vector<double>(gridsize);
        v = Vector<double>(gridsize);
        for(int i = 0; i < gridsize; i++) x[i] = i*2/gridsize - 1; // x = [-1, 1)

        // Arg 4: Time step (dt)
        istringstream s4(argv[4]);
        s4 >> dt;

        // Arg 5: Number of blocks
        istringstream s5(argv[5]);
        s5 >> numberOfBlocks;

        // Arg 6: Threads per block
        istringstream s6(argv[6]);
        s6 >> threadsPerBlock;

        // Arg 9: Device
        string s9(argv[9]);

        if(s9.compare("cpu") == 0) {
            device = Device.CPU;
        } else if(s9.compare("gpu") == 0) {
            device = Device.GPU;
        } else {
            cerr << "Invalid device \"" << s9 << "\". " << endl;
            return -1; 
        }

        // Arg 7: Basis function
        string s7(argv[7]);
        auto forcing_function = [] (double x) -> double {
            return 0;
        }

        if(s7.compare("hat") == 0) {
            discretization = FEMHatDiscretization(device, x, forcing_function);
        } else if(s7.compare("hermite") == 0) {
            discretization = FEMHermiteDiscretization(device, x, forcing_function);
        } else {
            cerr << "Invalid hat function \"" << s7 << "\". " << endl;
            return -1; 
        }

        // Arg 8: Solver
        string s8(argv[8]);

        if(s8.compare("direct") == 0) {
            solver = SolverDirect(discretization);
        } else if(s8.compare("lu") == 0) {
            solver = SolverLU(discretization);
        } else if(s8.compare("cg") == 0) {
            solver = SolverCG(discretization);
        } else {
            cerr << "Invalid solver \"" << s8 << "\". " << endl;
            return -1; 
        }
    } catch(invalid_argument const &ex) {
        cerr << "Invalid input argument. " << endl;
        return -1;
    } catch(out_of_range const &ex) {
        cerr << "Input argument out of range" << endl;
        return -1;
    }

    // Initialize vectors
    for(int i = 0; i < gridsize; i++) {
        x[i] = 1f/sqrt(2*M_PI)*exp(-x[i]*x[i]/2);
        v[i] = 0;
    }

    // Time step
    if(plot) {
        const int stepsPerFrame = 30/dt;
        const int frames = steps/stepsPerFrame;

        for(int i = 0; i < frames; i++) {
            rk4<double>(stepsPerFrame, x, v, discretization, solver, dt);

            // TODO: Export x (and v) to file for each frame. 
        }
    } else {
        rk4<double>(steps, x, v, discretization, solver, dt);
    }
}
