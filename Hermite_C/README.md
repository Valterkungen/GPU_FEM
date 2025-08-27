# C++ Hermite Implementation

This directory contains a C++ implementation of Hermite polynomial basis functions for finite element methods.

## Contents

- **`hermite_fem.cpp`**: Main C++ implementation
- **`hermite_bases.h`**: Hermite polynomial basis function definitions
- **`Makefile`**: Build configuration

## Building

```bash
make
```

## Purpose

This implementation serves as:
- A reference implementation for Hermite polynomials
- Performance baseline for CPU computations
- Bridge between Python and CUDA implementations

## Mathematical Background

Hermite polynomials provide higher-order basis functions with:
- Better approximation properties than hat functions
- Continuous derivatives
- Orthogonality properties useful for certain problems

The implementation includes efficient evaluation of Hermite polynomials and their derivatives for FEM assembly processes.
