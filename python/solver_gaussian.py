from interfaces import *
import numpy
class GaussianSolver(Solver):
    def __init__(self, discretization: FEMDiscretization):
        """
        Gausian solver. 

        ## Parameters
         - ``discretization`` A FEM Discretization. 
        """
        self.discretization = discretization

        # Do other setup if requried here
    # Function to perform partial pivot for Gaussian elimination

    def step(self, n_steps: int, x: np.ndarray | csr_matrix, dt: float) -> np.ndarray | csr_matrix:
        # Function to perform partial pivot for Gaussian elimination
        def gaussian_elimination(A,b):
            def partial_pivot(A, n):
                # Iterate through each row in the matrix
                for i in range(n):
                    pivot_row = i
                    # Find the row with the maximum absolute value in the current column
                    for j in range(i + 1, n):
                        if abs(A[j][i]) > abs(A[pivot_row][i]):
                            pivot_row = j
                    # Swap the current row with the row having the maximum absolute value
                    if pivot_row != i:
                        A[[i, pivot_row]] = A[[pivot_row, i]]
                    # Perform Gaussian elimination on the matrix
                    for j in range(i + 1, n):
                        factor = A[j][i] / A[i][i]
                        A[j] -= factor * A[i]

            # Function to perform back substitution to solve the system of equations
            def back_substitute(A, n):
                x = np.zeros(n)
                # Iterate through each row in reverse order
                for i in range(n - 1, -1, -1):
                    sum_val = sum(A[i][j] * x[j] for j in range(i + 1, n))
                    # Solve for x[i] using the previously calculated values of x
                    x[i] = (A[i][n] - sum_val) / A[i][i]
                return x
            n = len(A)
            # Coefficient matrix augmented with the constant terms
            A = np.concatenate((A,b.reshape(-1, 1)), axis=1)
            # Perform Gaussian elimination with partial pivot
            partial_pivot(A, n)
            x = back_substitute(A, n)
            return x
        """
        Step solution ``n_steps`` steps with timestep ``dt``. 

        ## Parameters
         - ``n_steps`` Number of steps. 
         - ``dt`` Timestep. 
        """
        h = 0.01
        t = 0
        M, A, F = self.discretization.matrices


        def rk4_step(t, Z, h, F, M, A):
            k1 = F(t, Z, M, A)
            k2 = F(t + 0.5 * h, Z + 0.5 * h * k1, M, A)
            k3 = F(t + 0.5 * h, Z + 0.5 * h * k2, M, A)
            k4 = F(t + h, Z + h * k3, M, A)
            return Z + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

        def F(t, Z, M, A):
            n = len(Z) // 2
            X, Y = Z[:n], Z[n:]
            dXdt = Y
            dYdt = - gaussian_elimination(M,A @ X)
            return np.concatenate([dXdt, dYdt])
        
        Z = np.concatenate([x,np.zeros_like(x)])

        with Timer('Gaussian solver'):
            for _ in range(n_steps):
                # Do loop iteration here
                Z = rk4_step(t, Z, h, F, M, A)
                t = t + h

                pass
        n = len(Z) // 2
        X, Y = Z[:n], Z[n:]
        return X