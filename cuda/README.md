# CUDA Implementation

This directory contains the GPU-accelerated implementation of the FEM solver using CUDA.

## Architecture

The CUDA implementation follows a modular design:

- **`main.cu`**: Main application entry point with command-line interface
- **`fem_hat.cu/cuh`**: Hat function basis implementation with CUDA kernels
- **`fem_hermite.cu/cuh`**: Hermite polynomial basis with GPU acceleration
- **`solver_*.cu/cuh`**: Various linear algebra solvers (LU, CG, Direct)
- **`linalg.cu/cuh`**: Core linear algebra operations optimized for GPU
- **`rk4.cu/cuh`**: Runge-Kutta 4th order time integration
- **`interfaces.cu/cuh`**: Abstract interfaces for extensibility

## Building

```bash
make
```

This creates the `app` executable.

## Usage

```bash
./app <plot> <n_steps> <grid_size> <dt> <blocks> <threads> <basis> <solver> <device>
```

### Example Runs

**Basic GPU computation:**
```bash
./app no 1000 256 0.01 32 256 hat lu gpu
```

**Performance comparison (CPU vs GPU):**
```bash
./app no 1000 256 0.01 32 256 hat lu cpu
./app no 1000 256 0.01 32 256 hat lu gpu
```

**Different basis functions:**
```bash
./app no 1000 256 0.01 32 256 hermite lu gpu
```

## Performance Tuning

- **Grid Size**: Larger grids show more GPU advantage
- **Blocks/Threads**: Optimize for your specific GPU architecture
- **Solver Choice**: LU decomposition typically fastest for dense systems
- **Basis Functions**: Hat functions are computationally simpler than Hermite

## CUDA Kernel Optimization

The implementation includes several GPU optimization techniques:
- Coalesced memory access patterns
- Shared memory utilization
- Thread block optimization
- Efficient reduction operations