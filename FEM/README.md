# Finite Element Method Theory and Examples

This directory contains theoretical background, derivations, and example implementations for the FEM solver project.

## Contents

### Jupyter Notebooks

- **`FEM.ipynb`**: Core finite element method theory
  - Weak form derivations
  - Basis function mathematics
  - Discretization techniques
  - Time integration methods

- **`Poisson_FEM.ipynb`**: Poisson equation examples
  - 1D and 2D implementations
  - Boundary condition handling
  - Convergence analysis
  - Visualization examples

### Implementation Files

- **`WaveEq.py`**: Wave equation solver implementation
- **`fem_hat.py`**: Hat function basis implementation
- **`fem_hermite.py`**: Hermite polynomial basis functions
- **`hermite_polynomials.py`**: Hermite polynomial utilities
- **`*.cu` files**: CUDA implementations for comparison

## Mathematical Background

The notebooks provide comprehensive coverage of:

1. **Weak Formulation**: Converting PDEs to weak forms suitable for FEM
2. **Basis Functions**: Mathematical properties of hat and Hermite functions
3. **Assembly Process**: Building global matrices from local elements
4. **Time Integration**: Runge-Kutta and other time-stepping schemes
5. **Boundary Conditions**: Dirichlet and periodic boundary implementations

## Usage

Start Jupyter and open the notebooks:

```bash
jupyter notebook FEM.ipynb
```

The notebooks are self-contained with explanations, mathematical derivations, and executable code examples.

## Theory References

The implementations follow standard FEM literature and include:
- Variational formulations
- Galerkin methods
- Finite element spaces
- Error analysis foundations
