# GPU Acceleration for Finite Element Method Solvers

*Undergraduate Research Project in Engineering Physics*

## Overview

This repository contains a comprehensive implementation and performance analysis of GPU-accelerated Finite Element Method (FEM) solvers. The project explores the computational advantages of CUDA-based parallel computing for solving partial differential equations using finite element methods.

### Key Features

- **Multiple Implementations**: Both Python and CUDA implementations for performance comparison
- **Various Basis Functions**: Support for hat functions and Hermite polynomials
- **Multiple Solvers**: Direct, LU decomposition, and Conjugate Gradient solvers
- **Performance Analysis**: Comprehensive benchmarking of CPU vs GPU performance
- **Flexible Architecture**: Modular design allowing easy extension and modification

## Repository Structure

```
├── cuda/                   # CUDA implementation
│   ├── main.cu            # Main CUDA application
│   ├── fem_hat.cu         # Hat function basis implementation
│   ├── fem_hermite.cu     # Hermite polynomial basis implementation
│   ├── solver_*.cu        # Various solver implementations
│   └── Makefile           # Build configuration
├── python/                # Python reference implementation
│   ├── main.py            # Main Python application
│   ├── fem_*.py           # FEM discretization modules
│   ├── solver_*.py        # Solver implementations
│   └── legacy/            # Legacy implementations
├── FEM/                   # Theoretical notes and Jupyter notebooks
│   ├── FEM.ipynb          # FEM theory and examples
│   └── Poisson_FEM.ipynb  # Poisson equation examples
├── Hermite_C/             # C++ Hermite implementation
├── Hermite_Cuda/          # CUDA Hermite implementation
└── docs/                  # Documentation
```

## Mathematical Background

This project focuses on solving the wave equation using finite element methods:

$$\frac{\partial^2 u}{\partial t^2} - \frac{\partial^2 u}{\partial x^2} = f(t), \quad x \in (0, L), \quad t \in (0, T)$$

With boundary conditions:
$$u(0, t) = u(L, t) = 0$$

And initial conditions:
$$u(x, 0) = u_0(x), \quad \frac{\partial u}{\partial t}(x, 0) = u_1(x)$$

## Building and Running

### Prerequisites

- **CUDA Development**: NVIDIA GPU with CUDA support, CUDA Toolkit (11.0+)
- **C++ Compilation**: GCC compiler
- **Python Environment**: Python 3.10+, NumPy, SciPy, Matplotlib

### CUDA Implementation

Navigate to the `cuda/` directory and build:

```bash
cd cuda/
make
```

Run the CUDA solver:

```bash
./app <plot> <n_steps> <grid_size> <dt> <blocks> <threads> <basis> <solver> <device>
```

**Parameters:**
- `plot`: Export data for visualization (yes/no)
- `n_steps`: Number of time steps
- `grid_size`: Number of grid points
- `dt`: Time step size
- `blocks`: Number of CUDA blocks
- `threads`: Threads per CUDA block
- `basis`: Basis function type (hat/hermite)
- `solver`: Solver type (direct/cg/lu)
- `device`: Computation device (cpu/gpu)

**Example:**
```bash
./app no 1000 256 0.01 32 256 hat lu gpu
```

### Python Implementation

Navigate to the `python/` directory and run:

```bash
cd python/
python3 main.py <dim> <bc> <n_steps> <dt> <grid_size> <basis> <solver>
```

**Parameters:**
- `dim`: Problem dimensions (1 or 2)
- `bc`: Boundary conditions (dirichlet/periodic)
- `n_steps`: Number of time steps
- `dt`: Time step size
- `grid_size`: Number of grid points
- `basis`: Basis function (hat/hermite)
- `solver`: Solver type (gaussian/lu/cg)

**Example:**
```bash
python3 main.py 1 dirichlet 1000 0.01 256 hat lu
```

## Performance Analysis

The project includes comprehensive performance analysis comparing:

- **CPU vs GPU**: Execution time comparison for different problem sizes
- **Solver Efficiency**: Performance of different linear algebra solvers
- **Basis Function Impact**: Computational cost of different basis functions
- **Memory Usage**: GPU memory utilization patterns

## Development Environment

For development setup instructions, see:
- [Development Environment Setup](docs/wslhost.md)
- [Remote Development Guide](docs/sshclient.md)

## Theoretical Background

The `FEM/` directory contains Jupyter notebooks with:
- Finite Element Method theory
- Weak form derivations
- Numerical implementation details
- Example problems and solutions

## Contributing

This project uses a modular architecture allowing easy extension:

1. **Adding New Basis Functions**: Implement in both `python/fem_*.py` and `cuda/fem_*.cu`
2. **Adding New Solvers**: Follow the solver interface in `interfaces.py` and `interfaces.cuh`
3. **Performance Testing**: Use the provided benchmarking framework

## License

This project is part of an undergraduate research thesis. Please contact the authors for usage permissions.

---

*This implementation demonstrates the significant performance advantages of GPU acceleration for finite element computations, with speedups of up to 10-100x depending on problem size and complexity.*
