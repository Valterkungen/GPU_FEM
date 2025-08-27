# Python Reference Implementation

This directory contains the Python reference implementation of the FEM solver, serving as a baseline for performance comparison and algorithm verification.

## Architecture

The Python implementation mirrors the CUDA structure for easy comparison:

- **`main.py`**: Main application with command-line interface
- **`fem_hat.py`**: Hat function basis implementation using NumPy
- **`fem_hermite.py`**: Hermite polynomial basis functions
- **`solver_*.py`**: Various linear algebra solvers (Gaussian, LU, CG)
- **`interfaces.py`**: Abstract base classes for extensibility
- **`colorprint.py`**: Colored terminal output utilities

## Dependencies

```bash
pip install numpy scipy matplotlib
```

## Usage

```bash
python3 main.py <dim> <bc> <n_steps> <dt> <grid_size> <basis> <solver>
```

### Example Runs

**1D Dirichlet problem:**
```bash
python3 main.py 1 dirichlet 1000 0.01 256 hat lu
```

**2D problem with periodic boundaries:**
```bash
python3 main.py 2 periodic 500 0.01 128 hermite cg
```

**Solver comparison:**
```bash
python3 main.py 1 dirichlet 1000 0.01 256 hat gaussian
python3 main.py 1 dirichlet 1000 0.01 256 hat lu
python3 main.py 1 dirichlet 1000 0.01 256 hat cg
```

## Performance Notes

- **NumPy Optimization**: Uses vectorized operations where possible
- **SciPy Integration**: Leverages optimized linear algebra routines
- **Memory Efficiency**: Careful memory management for large problems
- **Baseline Reference**: Provides verified results for CUDA comparison

## Legacy Code

The `legacy/` subdirectory contains earlier implementations and experimental code that may be useful for reference or alternative approaches.
