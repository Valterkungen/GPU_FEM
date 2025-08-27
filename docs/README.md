# Documentation

This directory contains setup and development documentation for the FEM GPU acceleration project.

## Contents

- **`wslhost.md`**: Development environment setup guide
- **`sshclient.md`**: Remote development configuration

## Quick Start

For setting up a CUDA development environment:

1. Install CUDA Toolkit (11.0+)
2. Verify GPU support with `nvidia-smi`
3. Clone and build the project
4. Run performance comparisons

## Development Workflow

1. **Local Development**: Use the Python implementation for algorithm development
2. **GPU Testing**: Port to CUDA for performance optimization
3. **Benchmarking**: Compare CPU vs GPU performance across different problem sizes
4. **Profiling**: Use CUDA profiling tools for optimization

## Performance Guidelines

- Start with small problem sizes for debugging
- Use larger grids (1000+ points) to see GPU advantages
- Profile memory usage and kernel efficiency
- Compare different solver algorithms for your specific use case
