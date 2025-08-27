"""
Main script. 
"""

import sys
from colorprint import print_error, print_info, print_sucess
from interfaces import *

from fem_hat import FEMHatDiscretization
from fem_hermite import FEMHermiteDiscretization

from solver_gaussian import GaussianSolver
from solver_lu import LUSolver
from solver_cg import ConjugateGradientSolver

if sys.version_info < (3, 10):
    raise SystemError('Python version >= 3.10 is requried to run the program. ')

def main():
    # Input arguments
    if len(sys.argv) != 8:
        print_error('Incorrect number of arguments. ')
        print_info('Usage: python3 main.py <dim> <bc> <n steps> <dt> <grid size> <basis function> <solver>')
        exit(-1)

    dimensions = int(sys.argv[1])
    boundary_condition = sys.argv[2]
    n_steps = int(sys.argv[3])
    dt = float(sys.argv[4])
    grid_size = int(sys.argv[5])
    basis_function = sys.argv[6]
    solver = sys.argv[7]

    SelectedDiscretization: FEMDiscretization = None
    SelectedSolver: Solver = None

    match boundary_condition.lower():
        case 'dirichlet':
            boundary_condition = BoundaryCondition.dirichlet
        case 'periodic':
            boundary_condition = BoundaryCondition.periodic
        case _:
            print_error(f'Invalid boundary condition "{boundary_condition}"')
            print_info('Boundary conditions: dirichlet, periodic')
            exit(-1)
    
    match basis_function.lower():
        case 'hat':
            SelectedDiscretization = FEMHatDiscretization
        case 'hermite':
            SelectedDiscretization = FEMHermiteDiscretization
        case _:
            print_error(f'Invalid basis function "{basis_function}"')
            print_info('Basis functions: hat, hermite')
            exit(-1)

    match solver.lower():
        case 'gaussian':
            SelectedSolver = GaussianSolver
        case 'lu': 
            SelectedSolver = LUSolver
        case 'cg':
            SelectedSolver = ConjugateGradientSolver
        case _:
            print_error(f'Invalid solver "{solver}"')
            print_info('Basis functions: gaussian, lu, cg')
            exit(-1)
    
    # Initilize grid
    match dimensions:
        case 1:
            x = np.linspace(-1, 1, grid_size)
            x = 1/np.sqrt(2*np.pi)*np.exp(-x**2/2)
        case 2:
            x = np.fromfunction(
                lambda i, j: np.sqrt((2*i/grid_size-1)**2 + (2*j/grid_size-1)**2), 
                (grid_size, grid_size)
            )
            x = 1/np.sqrt(2*np.pi)*np.exp(-x**2/2)
        case _:
            print_error(f'Unsupported number of dimensions {dimensions}')
    f = lambda x: 0

    # Discretize
    discreization = SelectedDiscretization(dimensions, boundary_condition, x, f)

    # Solve
    solver = SelectedSolver(discreization)
    x = solver.step(n_steps, x, dt)

    # Export result
    pass

    # Print success
    print_sucess("Program finished successfully")
    
if __name__ == '__main__':
    main()
