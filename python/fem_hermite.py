from interfaces import *
import hermite_polynomials as hp
import hermitebases as base
import fractions
np.set_printoptions(formatter={'all':lambda x: str(fractions.Fraction(x).limit_denominator())})

def f(xvec):
    return np.exp(-xvec**2)*0

"""
def LoadProduct(M,dim,f):
    return F
"""

def GenMatrices(M,dim):

    Msub = np.zeros((dim,dim))
    Asub = np.zeros((dim,dim))
    
    for i in range(dim):

        coeffs1 = base.getPoly(i,dim,M)
        p1 = hp.Polynomial(coeffs1)
        p1prim = p1.derive()

        for j in range(dim):

            coeffs2 = base.getPoly(j,dim,M)
            p2 = hp.Polynomial(coeffs2)
            p2prim = p2.derive()

            Msub[i][j] = (p1*p2).integrate(0,1)
            Asub[i][j] = (p1prim*p2prim).integrate(0,1)

    return Msub, Asub

class FEMHermiteDiscretization(FEMDiscretization):
    def __init__(self, dimensions: int, boundary_condition: BoundaryCondition, x: np.ndarray | csr_matrix, f: callable):
        """
        Fem Discretizer with hermite basis function. 
        ## Parameters
         - ``dimenions`` Number of dimenions, 1D or 2D. 
         - ``boundary_condition`` Boundary condition, dirichlet or periodic. 
         - ``x`` Grid vector/matrix. 
         - ``f`` Forcing function. 
        """
        self.dimensions = dimensions
        self.boundary_condition = boundary_condition

        match dimensions:
            case 1:
                with Timer('1D hermite discretization'):
                    self.__discretize1d(x, f)
            case 2:
                with Timer('2D hermite discretization'):
                    self.__discretize2d(x, f)
            case _:
                raise ValueError(f'FEM Hat Discrtization for dimension number {dimensions} is not supported. ')
    
    def __discretize1d(self, x: np.ndarray, f: callable) -> None:
        # Do stuff here...

        if self.boundary_condition == BoundaryCondition.periodic:
            pass
        
        dim = 4
        # A is hermite basis matrix
        Mh = base.genA(dim)
        # invert A to get M where columns in M are hermite basis functions
        H = np.linalg.inv(Mh)
        # M = base.round_matrix(M,digits)

        Msub,Asub = GenMatrices(H,dim)
        print(Msub)
        print(Asub)

        #grid
        N = len(x)
        nel = N-1
        ndof = 2*N
        h = 1/nel

        M = np.zeros((ndof,ndof))
        A = np.zeros((ndof,ndof))
        F = np.zeros(ndof)

        # Assembly loop
        for k in range(nel):
            kl, ku = 2*k, (2*k + 4)
            M[kl:ku,kl:ku] += Msub*(1.0/h)
            A[kl:ku,kl:ku] += Asub*(1.0/h)

        BC = 'Periodic'
        if BC == "Periodic":
            idr = ndof - 2
            A[0,:] = A[idr,:] = 0
            A[0,0] = A[idr,idr] = 1.0
            M[0,:] = M[idr,:] = 0
            M[0,0] = M[idr,idr] = 1.0

        if BC == "Dirichlet":
            idr = ndof - 2
            A[0,:] = A[idr,:] = 0
            A[0,0] = A[idr,idr] = 1.0
            M[0,:] = M[idr,:] = 0
            M[0,0] = M[idr,idr] = 1.0

        self.M, self.A = GenMatrices(H, dim)
        self.F = None
    
    def __discretize2d(self, x: csr_matrix, f: callable) -> None:
        # Do stuff here...

        if self.boundary_condition == BoundaryCondition.periodic:
            pass

        self.M = None
        self.A = None
        self.F = None


