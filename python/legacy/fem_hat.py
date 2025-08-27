import numpy as np
def fem_hat(x, f, boundry = "Perodic"):
    N = len(x) - 1
    h = x[1] - x[0]
    #Initilize Mass, Stiffness matricies and load vector
    M = np.zeros((N+1, N+1))
    A = np.zeros((N+1, N+1))
    F = np.zeros((N+1))

    #Definition of basis
    basis = {"eval": {"phi0": lambda x: x, "phi1": lambda x: 1-x}, 
        "nabla": {"phi0_nabla": lambda x: 1, "phi1_nabla": lambda x: -1}}
    #Integration using simpson
    quad = [(0, 1/6), (0.5, 4/6), (1, 1/6)]

    for k in range(N):
        for base in basis.keys():   
            for j, phi_j in enumerate(basis[base].values()):
                for iter, weight in quad:
                    if base == "eval":
                        F[k + j] += h * weight * phi_j(iter) * f(x)
                    if base == "eval":
                        for i, phi_i in enumerate(basis[base].values()):    
                            M[k + j , k + i] += 1/h * weight * phi_j(iter) * phi_i(iter)
                    if base == "nabla":
                        for i, phi_i in enumerate(basis[base].values()):    
                            A[k + j , k + i] += 1/h * weight * phi_j(iter) * phi_i(iter)


    if boundry == "Periodic":
        A[0, 0] = A[N, N] = A[1, 1]
        M[0, 0] = M[N, N] = M[1, 1]
        A[0, -1] = A[0, 1]
        A[N, 0] = A[N, -2]
        M[0, -1] = M[0, 1]
        M[N, 0] = M[N, -2]
    if boundry == "Dirichlet":
        A[0, :] = A[N, :] = 0
        M[0, :] = M[N, :] = 0
        A[0, 0] = A[N, N] = 1
        M[0, 0] = M[N,N] = 1
        F[0] = F[N] = 0
    return M, A, F