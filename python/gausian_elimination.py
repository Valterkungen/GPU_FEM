import numpy as np

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
    
    
    n = 3
    # Coefficient matrix augmented with the constant terms
    A = np.concatenate((A,b), axis=1)
    # Perform Gaussian elimination with partial pivot
    partial_pivot(A, n)
    x = back_substitute(A, n)
    return x



A = np.array([[3.0, 2.0, -4.0],
                [2.0, 3.0, 3.0],
                [5.0, -3, 1.0]])
b = np.array([[3.0],
              [15.0],
              [14.0]])

x = gaussian_elimination(A,b)
print(x.shape)