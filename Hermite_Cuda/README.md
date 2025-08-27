# CUDA Hermite Implementation

This directory contains a specialized CUDA implementation focusing on Hermite polynomial basis functions.

## Contents

- **`hermite_fem.cu`**: CUDA kernels for Hermite polynomial evaluation
- **`hermite_bases.cuh`**: CUDA header with Hermite basis definitions
- **`Makefile`**: Build configuration

## Building

```bash
make
```

## Features

- **GPU-Optimized Hermite Evaluation**: Efficient parallel computation of Hermite polynomials
- **Memory-Optimized Kernels**: Minimized memory transfers between host and device
- **High-Order Accuracy**: Support for higher-order Hermite polynomials

## Performance Characteristics

Hermite polynomials on GPU show significant advantages for:
- Large problem sizes (>1000 grid points)
- High-order polynomial evaluations
- Problems requiring smooth basis functions

## Integration

This implementation can be integrated with the main CUDA solver by selecting the Hermite basis option in the main application.
