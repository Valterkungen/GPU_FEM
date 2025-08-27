#include "solver_lu.cuh"
#include "main.cuh"


void SolverLU::CPUdecomp() {
    int n = this->discretization.M.getRows();
    this->L = DenseMatrix(Device.CPU, this->discretization.M.getRows(), this->discretization.M.getCols());
    this->U = DenseMatrix(Device.CPU, this->discretization.M.getRows(), this->discretization.M.getCols());
    this->M = this->discretization.M;
    auto M = this->M.toDense();
    for (int i = 0; i < n; i++){
        L[i][i] = 1.0; // set diagonal of L to 1
        for (int j = i; j < n; j++) {
            float sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += L[i][k] * U[k][j];
            }
            U[i][j] = M[i][j] - sum; // Computes U[i,j]
        }
        for (int j = i + 1; j < n; j++) {
            float sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += L[j][k] * U[k][i];
            }
            L[j][i] = (M[j][i] - sum) / U[i][i]; // Computes L[j,i]
        }
    }
}

__global__ void LUdecomp(int n, float *L, float *U, float *A){
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if(i < n){
        L[i*n + i] = 1.0; // set diagonal of L to 1
        for(int j = 0; j<i; j++ ){
            float sum = 0.0;
            for(int k = 0; k < i; k++){
                sum += L[i * n + k] * U[k * n + j]; 
            }
            U[i * n + j] = A[i * n + j] - sum; //Computes U[i,j]
        }
        for(int j = i + 1; j < n; j++ ){
            float sum = 0.0;
            for(int k = 0; k < i; k++){
                sum += L[j*n + k] * U[k*n +1];
            }
            L[j*n + i] = (A[j*n + i] - sum) / U[i*n + i]; //Computes L[i,j]
        }
    }
}

void SolverLU::GPUdecomp() { 
    // Update L and U. No need to return a value since we pass a refrence as a parameter. 
    int n = this->discretization.M.getRows();
    this->L = DenseMatrix(Device.GPU, this->discretization.M.getRows(), this->discretization.M.getCols());
    this->U = DenseMatrix(Device.GPU, this->discretization.M.getRows(), this->discretization.M.getCols());
    this->M = this->discretization.M; 
    auto M = this->M.toDense();
    dim3 block_shape(256);
    dim3 grid_shape((n + block_shape.x - 1) / block_shape.x);

    LUdecomp<<<grid_shape,block_shape>>>(n,&L.data,&U.data,&M.data);
}

Vector<T> &SolverLU::CPUsolve(Vector<T> &rhs) {
    // TODO: Add CPU implementation
    int n = this->discretization.M.getRows(); 
    Vector<T> y = Vector(Device.CPU, n); 
    Vector<T> x = Vector(Device.CPU, n);  

    //Forward substitution to solve Ly = rhs
    for (int i = 0; i < n; i++) {
        float sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];  // Access elements using two-dimensional indexing
        }
        y[i] = (rhs[i] - sum) / L[i][i];  // Access diagonal element of L directly
    }

    //Backward substitution to solve Ux = y
    for (int i = n - 1; i >= 0; i--) {
        float sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += U[i][j] * x[j];  // Correct indexing to access elements in U and y
        }
        x[i] = (y[i] - sum) / U[i][i];  // Directly access the diagonal element U[i][i]
    }
    return x;
}



Vector<T> &SolverLU::GPUsolve(Vector<T> &rhs) {
    // TODO: Add GPU implementation
    int n = rhs.getRows();
    Vector<T> y = Vector(Device.GPU, n);
    Vector<T> x = Vector(Device.GPU, n); 

    //Forward substitution to solve Ly = rhs
    for (int i = 0; i < n; i++) {
        //float sum = 0.0;
        //for (int j = 0; j < i; j++) {
        //    sum += L[i][j] * y[j];  // Access elements using two-dimensional indexing
        
        y[i] = (rhs[i] - RowVector(L,i)(0,i)*y(0:i)) / L[i][i];  // Access diagonal element of L directly
    }

    //Backward substitution 
    for (int i = n - 1; i >= 0; i--) {
        //float sum = 0.0;
        //for (int j = i + 1; j < n; j++) {
            //sum += U[i][j] * x[j];  // Correct indexing to access elements in U and y
        // }

        x[i] = (y[i] - RowVector(U,i)(i,n)*x(i,n)) / U[i][i];  // Directly access the diagonal element U[i][i]
    }
    
    return x;
}

SolverLU::SolverLU(FEMDiscretization &discretization): discretization(discretization) {
    this->device = discretization.device;

    switch(this->device) {
        case CPU:
            this->CPUdecomp<T>(this->L, this->U);
            break;
        case GPU:
            this->GPUdecomp<T>(this->L, this->U);
            break;
    }
}

Vector<T> &SolverLU::solve(Vector<T> rhs) {
    switch(this->device) {
        case CPU:
            this->CPUsolve<T>(rhs);
            break;
        case GPU:
            this->GPUsolve<T>(rhs);
            break;
    }
}