import numpy as np
def fem_hat(x, f, boundry = "Periodic"):
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
    

    dim = len(basis)
    local_mass = np.zeros((dim, dim))
    local_stiffnes = np.zeros((dim, dim))


    for base in basis.keys():   
        for j, phi_j in enumerate(basis[base].values()):
            for iter, weight in quad:
                if base == "eval":
                    for i, phi_i in enumerate(basis[base].values()):    
                        local_mass[j , i] += weight * phi_j(iter) * phi_i(iter)
                if base == "nabla":
                    for i, phi_i in enumerate(basis[base].values()):    
                        local_stiffnes[j , i] += weight * phi_j(iter) * phi_i(iter)


    DOF = dim
    A2 = np.zeros((N + 1, N + 1))
    M2 = np.zeros((N + 1, N+1))
    h = 1
    for k in range(N):
        kl, ku = k, k + 2
        A2[kl:ku,kl:ku] += (1.0/h)*local_stiffnes
        M2[kl:ku,kl:ku] += (1.0/h)*local_mass

    A2[0, N] = A2[N, 0] = -1/h
    A2[0, 0] += 1/h
    A2[N, N] += 1/h

    M2[0, N] = M2[N, 0] = 1/6*h
    M2[0, 0] += 1/3*h
    M2[N, N] += 1/3*h
    print(M2)
    if boundry == "Dirichlet":
        A[0, :] = A[N, :] = 0
        M[0, :] = M[N, :] = 0
        A[0, 0] = A[N, N] = 1
        M[0, 0] = M[N,N] = 1
        F[0] = F[N] = 0
    # om inte Dirichlet? 
    return M, A, F


x = np.linspace(0, 1, 6)
fem_hat(x, lambda x: 0)