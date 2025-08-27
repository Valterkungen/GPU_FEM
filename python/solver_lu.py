

import scipy
import scipy.linalg
from interfaces import FEMDiscretization, Solver
import numpy as np
from scipy.sparse import csr_matrix


class LUSolver(Solver):
    def __init__(self, discretization: FEMDiscretization):
        """
        LU solver. 

        ## Parameters
         - ``discretixation`` A FEM Discretization. 
        """
        self.discretization = discretization

        # Do other setup if requried here

    def step(self, n_steps: int, x: np.ndarray | csr_matrix, dt: float) -> np.ndarray | csr_matrix:
        """
        Step solution ``n_steps`` steps with timestep ``dt``. 
        
        ## Parameters
         - ``n_steps`` Number of steps. 
         - ``dt`` Timestep. 
        """
        M, A, F = self.discretization.matrices


        # Function for computing LU decomposition of matrix A
        def lu_factorization(A):
            n = len(A)
            L = np.zeros((n, n))
            U = np.zeros((n, n))

            for i in range(n):
                L[i, i] = 1
                for j in range(i, n):
                    U[i, j] = A[i, j] - np.dot(L[i, :i], U[:i, j])
                for j in range(i + 1, n):
                    L[j, i] = (A[j, i] - np.dot(L[j, :i], U[:i, i])) / U[i, i]

            return L, U
        
        # Function to solve the lower triangular system Ly = b
        def forward_substitution(L, b):
            n = L.shape[0]
            y = np.zeros(n)
            for i in range(n):
                y[i] = (b[i] - np.dot(L[i, :i], y[:i])) / L[i, i]
            return y

        # Function to solve the upper triangular system Ux = y
        def backward_substitution(U, y):
            n = U.shape[0]
            x = np.zeros(n)
            for i in range(n - 1, -1, -1):
                x[i] = (y[i] - np.dot(U[i, i+1:], x[i+1:])) / U[i, i]
            return x

        
        #Find LU factorisation
        L, U = lu_factorization(M)

        #Initial values
        x_0 = x
        v_0 = 0


        for _ in range(n_steps):
            # Do loop iteration here
            #the system to solve is LUvprim = -A_v
            # compute rhs b = -A@v_0
            b = -A@v_0
            #solving the system Ly = b 
            y = forward_substitution(L,b)
            #solving the system Uv_prim = y
            v_prim = backward_substitution(U, y)
            v = v_0 +dt*v_prim  #do rk4 instead
            x = x_0 + dt*v #do rk4 instead

            v_0 = v
            x_0 = x
            

        return x




def lu_factorization(A):
    n = len(A)
    L = np.zeros((n, n))
    U = np.zeros((n, n))

    for i in range(n):
        L[i, i] = 1
        for j in range(i, n):
            U[i, j] = A[i, j] - np.dot(L[i, :i], U[:i, j])
            #print('U=', U)
        for j in range(i + 1, n):
            L[j, i] = (A[j, i] - np.dot(L[j, :i], U[:i, i])) / U[i, i]
            #print('L=', L)
    return L, U

# Function to solve the lower triangular system Ly = b
def forward_substitution(L, b):
    n = L.shape[0]
    y = np.zeros(n)
    for i in range(n):
        y[i] = (b[i] - np.dot(L[i, :i], y[:i])) / L[i, i]
    return y

# Function to solve the upper triangular system Ux = y
def backward_substitution(U, y):
    n = U.shape[0]
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - np.dot(U[i, i+1:], x[i+1:])) / U[i, i]
    return x

# Example usage:
#A = np.random.rand(5, 5)
A = np.array([[4,-1,1],[-1,4,-2],[1,-2,4]])
#b = np.random.rand(5)
b = np.array([12,-1,5])
x_test = scipy.linalg.solve(A,b)
L, U = lu_factorization(A)
y = forward_substitution(L,b)
x = backward_substitution(U,y)
print(x,x_test)



