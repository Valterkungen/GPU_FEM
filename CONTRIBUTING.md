# Contributing to GPU-Accelerated FEM Solver

Thank you for your interest in contributing to this project! This guide will help you get started.

## Development Setup

1. **Clone the repository**
   ```bash
   git clone <repository-url>
   cd FEM
   ```

2. **Install dependencies**
   - CUDA Toolkit (11.0+)
   - Python 3.10+ with NumPy, SciPy, Matplotlib
   - GCC compiler

3. **Build and test**
   ```bash
   cd cuda && make
   cd ../python && python3 main.py 1 dirichlet 100 0.01 64 hat lu
   ```

## Project Structure

- **`cuda/`**: GPU implementations
- **`python/`**: Reference implementations
- **`FEM/`**: Theory and examples
- **`docs/`**: Documentation

## Adding New Features

### New Basis Functions
1. Implement in `python/fem_newbasis.py`
2. Add corresponding CUDA implementation in `cuda/fem_newbasis.cu`
3. Update interfaces in both implementations
4. Add tests and documentation

### New Solvers
1. Follow the solver interface pattern
2. Implement both CPU and GPU versions
3. Add performance benchmarks
4. Document algorithm complexity

## Code Style

- **Python**: Follow PEP 8
- **C++/CUDA**: Use consistent indentation and naming
- **Comments**: Document complex algorithms and GPU kernel optimizations
- **Headers**: Include brief descriptions of file purpose

## Testing

- Test both CPU and GPU implementations
- Verify numerical accuracy against reference solutions
- Include performance regression tests
- Test different problem sizes and configurations

## Performance Guidelines

- Profile GPU kernels using CUDA profiling tools
- Optimize memory access patterns
- Document performance characteristics
- Include benchmark results in pull requests

## Documentation

- Update relevant README files
- Add examples for new features
- Include mathematical background where applicable
- Update the main README if adding major features

## Pull Request Process

1. Create a feature branch from `main`
2. Make your changes with clear commit messages
3. Test thoroughly on different configurations
4. Update documentation
5. Submit pull request with detailed description

## Questions?

Feel free to open an issue for questions about:
- Implementation details
- Performance optimization
- Mathematical formulations
- CUDA programming techniques
