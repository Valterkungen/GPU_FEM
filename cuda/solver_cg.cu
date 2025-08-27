#include "solver_cg.cuh"
#include "main.cuh"
#include <cmath> // 


Vector<T> &SolverCG::CPUsolve(Vector<T> &rhs) {
    // TODO: Add CPU implementation
    return Vector<T>();
}

Vector<T>& SolverCG::GPUsolve(Matrix<T>& A, Vector<T>& rhs, Vector<T>* x0 = nullptr, float tol = 1e-6) {
    // Conjugate Gradient ska lösa ekv Ax = b
    // I FEM context är A = M (Mass matrix) & b är RHS 
    int vec_size = rhs.getLength();
    Vector<T> x(vec_size), x_new(vec_size), r(vec_size), p(vec_size), Ax(vec_size), w(vec_size);
    float rho, rho_prev, alpha, beta;

    // Initial Gissning eller sätt noll vektor
    if (x0 == nullptr) {
        x_new.setZero();  // Assuming setZero initializes all elements to zero
    } else {
        x_new = *x0;
    }

    // Beräkna residual r = b - Ax and initial riktning p = r
    Ax = A * x;  // Ax = A*x
    r = rhs - Ax;        // r = b - Ax
    p = r;               // p = r
    rho = r*r;      // rho = r^T * r

    float err = 2 * tol;

    while (err > tol) {
        x = x_new;
        w = A*p;  // w = A*p
        float dot_pAp = p*A;  // dot_pAp = p^T * Ap

        if (dot_pAp == 0) {
            std::cout << "Error: p^T * Ap = 0, stopping iteration." << std::endl;
            break;
        }

        alpha = rho / dot_pAp; // alpha = rho / (p^T * Ap)
        x_new = x + p * alpha;     // x = x + alpha * p
        r = r - w * alpha;    // r = r - alpha * Ap
        
        rho_prev = rho;
        rho = r*r;   // rho_new = r^T * r
        beta = rho / rho_prev;
        p = r + p * beta;     // p = r + beta * p

        // Update error
        err = sqrt((x_new - x) * (x_new - x)) / sqrt(x_new * x_new); //norm(x_new - x) / norm(x_new)

    }

    return x_new;
}

SolverCG::SolverCG(FEMDiscretization &discretization): discretization(discretization) {
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

Vector<T> &SolverCG::solve(Vector<T> rhs) {
    switch(this->device) {
        case CPU:
            this->CPUsolve<T>(rhs);
            break;
        case GPU:
            this->GPUsolve<T>(rhs);
            break;
    }
}