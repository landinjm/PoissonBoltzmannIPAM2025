"""
JuliaDraftSolver Python Wrapper

A Python package that provides a clean interface to the JuliaSolver Julia package
for solving 1D nonlinear Poisson equation.

The Julia environment is automatically set up when this package is imported.

Example usage:
    >>> from juliasolver_wrapper import NonlinearPoissonSolver
    >>> solver = NonlinearPoissonSolver()
    >>> x, u = solver.solve(n=31, domain=[0, 2], bcval=[0.5, 2.0], source=2, m=3)
"""
from .core import NonlinearPoissonSolver
