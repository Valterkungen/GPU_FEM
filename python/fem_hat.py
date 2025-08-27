from interfaces import *

class FEMHatDiscretization(FEMDiscretization):
    def __init__(self, dimensions: int, boundary_condition: BoundaryCondition, x: np.ndarray | csr_matrix, f: callable):
        """
        Fem Discretizer with hat basis function. 

        ## Parameters
         - ``dimenions`` Number of dimenions, 1D or 2D. 
         - ``boundary_condition`` Boundary condition, dirichlet or periodic. 
         - ``x`` Grid vector/matrix. 
         - ``f`` Forcing function. 
        """
        self.dimensions = dimensions
        self.boundary_condition = boundary_condition
        self.x = x
        self.f = f
        self.M, self.A, self.F = None, None, None

        match dimensions:
            case 1:
                with Timer('1D hat discretization'):
                    self.M, self.A, self.F = self.__discretize1d(x, f)
            case 2:
                with Timer('2D hat discretization'):
                    self.M, self.A, self.F = self.__discretize2d(x, f)
            case _:
                raise ValueError(f'FEM Hat Discrtization for dimension number {dimensions} is not supported. ')

    
    def __discretize1d(self, x: np.ndarray, f: callable) -> None:
        # Do stuff here...
        
        N = len(x) - 1
        h = x[1] - x[0]
        
        # Initilize Mass, Stiffness matricies and load vector
        M = np.zeros((N+1, N+1))
        A = np.zeros((N+1, N+1))
        F = np.zeros((N+1))

        # Definition of basis
        basis = {"eval": {"phi0": lambda x: x, "phi1": lambda x: 1-x}, 
            "nabla": {"phi0_nabla": lambda x: 1, "phi1_nabla": lambda x: -1}}

        # Integration using simpson
        quad = [(0, 1/6), (0.5, 4/6), (1, 1/6)]

        # Pre determined local stiffness and mass 
        local_stiffnes = np.array([[1, -1], [-1, 1]])
        local_mass = np.array([[1/3, 1/6],[1/6, 1/3]])

        for k in range(N):
            for base in basis.keys():   
                for j, phi_j in enumerate(basis[base].values()):
                    for iter, weight in quad:
                        if base == "eval":
                            F[k + j] += h * weight * phi_j(iter) * f(x)
        
        
        for k in range(N):
            kl, ku = k, k + 2
            A[kl:ku,kl:ku] += (1.0/h) * local_stiffnes
            M[kl:ku,kl:ku] += (1.0/h) * local_mass

        if self.boundary_condition == BoundaryCondition.periodic:
            A[0, N] = A[N, 0] = -1/h
            A[0, 0] += 1/h
            A[N, N] += 1/h

            M[0, N] = M[N, 0] = 1/6*h
            M[0, 0] += 1/3*h
            M[N, N] += 1/3*h

        if self.boundary_condition == BoundaryCondition.dirichlet:
            A[0, :] = A[N, :] = 0
            A[0, 0] = A[N, N] = 1

            M[0, :] = M[N, :] = 0
            M[0, 0] = M[N,N] = 1

            F[0] = F[N] = 0
        return M, A, F
    
    def __discretize2d(self, x: csr_matrix, f: callable) -> None:
        # Do stuff here...

        if self.boundary_condition == BoundaryCondition.periodic:
            pass

        self.M = None   
        self.A = None
        self.F = None

    @property
    def matrices(self):
        return self.M, self.A, self.F


