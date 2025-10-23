PoissonBoltzmannIPAM2025
=======================

Project to develop Modified Poisson-Boltzman solver for comparing with
molecular level simulations.

Some discussion how to work with this kind of repository are found in
[this Julia project HOWTO](https://j-fu.github.io/marginalia/julia/project-howto/).

## Initial contents
```
.
├── assets                       # Some data used by the demo code 
│   └── FedorovKornyshev.png
├── LICENSE
├── notebooks
│   ├── demo-notebook.jl         # Notebook template
│   ├── pluto-simplecell-bsk.jl  # MPB demo (WIP) with Bazant-Storey-Kornyshev model
│   └── pluto-simplecell.jl      # Poisson-Boltzmann/Bikermann demo
├── packages
│   └── JuliaMPBSolver           # Package for the code to be developed
│       ├── Project.toml
│       └── src
│           └── JuliaMPBSolver.jl
├── Project.toml                 # Julia dependency description
├── README.md
└── scripts
    └── demo-script.jl           # Script template
```
