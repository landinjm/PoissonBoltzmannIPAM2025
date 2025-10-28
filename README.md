PoissonBoltzmannIPAM2025
=======================

Project to develop Modified Poisson-Boltzman solver for comparing with
molecular level simulations.

## Getting started
First, we'll clone the repository and instantiate the packages. **Make sure you have julia 1.11!**
```
git clone https://github.com/IPAM-ECH2025/PoissonBoltzmannIPAM2025.git
```
```
julia -e "using Pkg; Pkg.activate(); Pkg.instantiate()"
```
You can then open one of the notebooks or execute one of the scripts.

## Contents
```
.
├── assets                       # Some data used by the demo code 
│   └── FedorovKornyshev.png
├── LICENSE
├── notebooks
│   ├── demo-notebook.jl         # Notebook template
│   ├── pluto-equilibrium.jl
│   ├── pluto-simplecell-bsk.jl  # MPB demo (WIP) with Bazant-Storey-Kornyshev model
│   └── pluto-simplecell.jl      # Poisson-Boltzmann/Bikermann demo
├── packages
│   └── JuliaMPBSolver           # Main package
│       ├── Project.toml
│       └── src
│           └── JuliaMPBSolver.jl
├── Project.toml                 # Julia dependency description
├── README.md
└── scripts
    └── demo-script.jl           # Script template
```
