from enum import Enum
from scipy.sparse import csr_matrix
import numpy as np
from typing import Callable
import sys
from time import perf_counter

if sys.version_info < (3, 10):
    raise SystemError('Python version >= 3.10 is requried to run the program. ')

BoundaryCondition = Enum('BoundaryConditions', ['dirichlet', 'periodic'])

class FEMDiscretization:
    def __init__(self, dimensions: int, boundary_condition: BoundaryCondition, x: np.ndarray | csr_matrix, f: Callable):
        raise NotImplementedError('Constructor not implemented. ')

    def __init_subclass__(cls):
        cls.dimensions: int = None
        cls.boundary_condition: BoundaryCondition = None

        cls.M: csr_matrix = None
        cls.A: csr_matrix = None
        cls.F: np.ndarray | csr_matrix = None

    @property
    def matrices(self):
        return self.M, self.A, self.F

class Solver:
    def __init__(self, discretization: FEMDiscretization):
        raise NotImplementedError('Constructor not implemented. ')

    def __init_subclass__(cls):
        cls.discretization = None
    
    def step(self, n_steps: int, x: np.ndarray | csr_matrix, dt: float) -> np.ndarray | csr_matrix:
        raise NotImplementedError('Step function not implemented. ')

class Timer:
    def __init__(self, name):
        self.name = name
    
    def __enter__(self):
        self.time = perf_counter()
    
    def __exit__(self, *errors):
        self.time = perf_counter() - self.time
        print(f'{self.name} took {self.time:.2f} s')
