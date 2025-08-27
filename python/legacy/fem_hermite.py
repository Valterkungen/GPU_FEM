import numpy as np
def fem_hat(x, f, boundry = "Perodic"):

    #Definition of basis
    basis = {"eval": 
        {"H1": lambda x: (2*x**3 - 3*x**2 + 1), 
         "H2": lambda x: (x**3 - 2*x**2 + x), 
         "H3": lambda x: (-2*x**3 + 3*x**2), 
         "H4": lambda x: (x**3 - x**2)
         },
         "nabla": 
         {"H1_prime": lambda x: (6*x**2 - 6*x), 
          "H2_prime": lambda x: (3*x**2 - 4*x + 1),
          "H3_prime": lambda x: (6*x - 6*x**2),
          "H4_prime": lambda x: (3*x**2 - 2*x)
         } 
    }
    N = len(x) - 1
    N_iter = len(x)
    h = x[1] - x[0]
    nel = N-1
    ndof = N*2
    h = 1.0/nel
    steps = 1/N_iter
    iter = np.linspace(0, 1, N_iter + 1)
    #Initilize Mass, Stiffness matricies and load vector
    M = np.zeros((ndof, ndof))
    A = np.zeros((ndof, ndof))
    F = np.zeros(ndof)

    dim = 2*len(basis)
    local_mass = np.zeros((dim, dim))
    local_stiffnes = np.zeros((dim, dim))
    



    
    for k in range((nel + 1)//2):
        k = 4*k
        for base in basis.keys():   
            for j, Hj in enumerate(basis[base].values()):
                if base == "eval":
                    fx = Hj(iter) * f(iter)
                    F[k + j] = 2*h * (steps/3)*(fx[0] + 4*np.sum(fx[1:-1:2]) + 2*np.sum(fx[2:-2:2]) +fx[-1])
                    for i, Hi in enumerate(basis[base].values()):    
                            fx = Hj(iter)*Hi(iter)
                            local_mass[j , i] = (steps/3)*(fx[0] + 4*np.sum(fx[1:-1:2]) + 2*np.sum(fx[2:-2:2]) +fx[-1])
                if base == "nabla":
                        for i, Hi in enumerate(basis[base].values()):    
                            fx = Hj(iter)*Hi(iter)
                            local_stiffnes[j , i] = (steps/3)*(fx[0] + 4*np.sum(fx[1:-1:2]) + 2*np.sum(fx[2:-2:2]) +fx[-1])

    # Assembly loop

    print(local_mass)
    for k in range(nel):

        kl, ku = 2*k, 2*k + 4
        A[kl:ku,kl:ku] += (1.0/h)*local_stiffnes
        M[kl:ku,kl:ku] += (1.0/h)*local_mass

    if boundry == "Periodic":
        idr = ndof - 2
        A[0,:] = A[idr,:] = 0
        A[0,0] = A[idr,idr] = 1.0
        M[0,:] = M[idr,:] = 0
        M[0,0] = M[idr,idr] = 1.0
        F[0] = F[idr]

    if boundry == "Dirichlet":
        idr = ndof - 2
        A[0,:] = A[idr,:] = 0
        A[0,0] = A[idr,idr] = 1.0
        M[0,:] = M[idr,:] = 0
        M[0,0] = M[idr,idr] = 1.0
        F[0] = F[idr] = 0

    return M, A, F

def main():
    xvec = np.linspace(0,1,10)
    fem_hat(xvec,xvec,0)

if __name__ == '__main__':
     main()