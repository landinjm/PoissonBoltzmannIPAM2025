"""
Core solver module for JuliaSolver Python wrapper.

This module provides the main NonlinearPoissonSolver class that wraps
the Julia implementation.
"""

import numpy as np

from typing import Tuple, List
from juliacall import Main as jl


class JuliaMPBSolver:
    """
    Python wrapper for the JuliaDraftSolver nonlinear Poisson equation solver.

    This class provides a clean Python interface to solve 1D nonlinear Poisson
    equations of the form: Δu^m = f
    """

    def __init__(self):
        """
        Initialize the NonlinearPoissonSolver.

        The Julia environment is automatically set up when the package is imported.
        """
        self._ensure_julia_ready()

    def _ensure_julia_ready(self):
        """Ensure Julia modules are loaded."""
        # Julia environment is already set up by package import
        # Just ensure the JuliaDraftSolver module is loaded
        jl.seval("using JuliaMPBSolver")

    def mpbpsolve(
        self,
        n: int = 21,
        domain: List[float] = [0, 10.0e-9],
        chargenumbers: List[float] = [-1, 1],
        bulkmolarity: float = 0.1,
        surfacecharge: List[float] = [0.16, -0.16],
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Solve the 1D nonlinear Poisson equation: Δu^m = f

        Args:
            n: Number of grid points
            domain: Domain boundaries [left, right]
            bcval: Dirichlet boundary values [left_bc, right_bc]
            source: Constant source term
            m: Nonlinearity exponent
            inival: Initial value for solver (if None, uses average of bcval)

        Returns:
            Tuple of (x_coordinates, solution_values) as numpy arrays

        Example:
            >>> solver = NonlinearPoissonSolver()
            >>> x, u = solver.solve(n=31, domain=[0, 2], bcval=[0.5, 2.0],
            ...                     source=2, m=3)
        """
        self._ensure_julia_ready()

        # Call the Julia function
        X, C0, Cp, Cm = jl.mpbpsolve(
            n=n,
            domain=domain,
            chargenumbers=chargenumbers,
            bulkmolarity=bulkmolarity,
            surfacecharge=surfacecharge,
        )

        # Convert to numpy arrays. Note the call to collect to make sure
        # we get an array rather than a Julia range.
        x_values = np.array(jl.collect(X), copy=True)
        c0_values = np.array(C0, copy=True)
        cp_values = np.array(Cp, copy=True)
        cm_values = np.array(Cm, copy=True)

        return x_values, c0_values, cp_values, cm_values

    def icmpbpsolve(
        self,
        n: int = 21,
        domain: List[float] = [0, 10.0e-9],
        chargenumbers: List[float] = [-1, 1],
        averagemolarity: float = 1,
        surfacecharge: List[float] = [0.16, -0.16],
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Solve the 1D nonlinear Poisson equation: Δu^m = f

        Args:
            n: Number of grid points
            domain: Domain boundaries [left, right]
            bcval: Dirichlet boundary values [left_bc, right_bc]
            source: Constant source term
            m: Nonlinearity exponent
            inival: Initial value for solver (if None, uses average of bcval)

        Returns:
            Tuple of (x_coordinates, solution_values) as numpy arrays

        Example:
            >>> solver = NonlinearPoissonSolver()
            >>> x, u = solver.solve(n=31, domain=[0, 2], bcval=[0.5, 2.0],
            ...                     source=2, m=3)
        """
        self._ensure_julia_ready()

        # Call the Julia function
        X, C0, Cp, Cm = jl.icmpbpsolve(
            n=n,
            domain=domain,
            chargenumbers=chargenumbers,
            averagemolarity=averagemolarity,
            surfacecharge=surfacecharge,
        )

        # Convert to numpy arrays
        x_values = np.array(X)
        c0_values = np.array(C0)
        cp_values = np.array(Cp)
        cm_values = np.array(Cm)

        return x_values, c0_values, cp_values, cm_values
