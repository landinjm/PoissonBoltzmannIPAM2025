"""
Core solver module for JuliaSolver Python wrapper.

This module provides the main NonlinearPoissonSolver class that wraps
the Julia implementation.
"""

import numpy as np

from typing import Tuple, List, Optional, Union
import sys

try:
    from juliacall import Main as jl
except ImportError:
    raise ImportError(
        "juliacall is required. Install with: pip install juliacall"
    )


class NonlinearPoissonSolver:
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
        jl.seval("using JuliaDraftSolver")
    
    def solve(
        self,
        n: int = 21,
        domain: List[float] = [0, 1],
        bcval: List[float] = [1, 1],
        source: float = 1,
        m: float = 2,
        inival: Optional[float] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
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
        
        if inival is None:
            inival = sum(bcval) / 2
        
        # Call the Julia function
        X, U = jl.solve_nonlinear_poisson(
            n=n,
            domain=domain,
            bcval=bcval,
            source=source,
            m=m,
            inival=inival
        )
        
        # Convert to numpy arrays
        x_values = np.array(X)
        u_values = np.array(U)
        
        return x_values, u_values
    
