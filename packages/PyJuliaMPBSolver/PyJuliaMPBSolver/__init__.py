import os
import sys
import warnings
from types import ModuleType
from typing import cast


def _import_juliacall():
    import juliacall


_import_juliacall()

from juliacall import Main as jl

jl = cast(ModuleType, jl)

jl.seval("using JuliaMPBSolver")
JuliaMPBSolver = jl.JuliaMPBSolver

jl.seval("using Pkg: Pkg")
Pkg = jl.Pkg

from .core import JuliaMPBSolver
