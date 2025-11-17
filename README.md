PoissonBoltzmannIPAM2025
=======================

Project to develop Modified Poisson-Boltzmann solver for comparing with
molecular level simulations.

We are solving a class of Poisson-Boltzmann systems as described below. First, we have poisson equation for charge density.
$$
-\nabla \cdot (\varepsilon(\phi, y_\alpha) \nabla \phi) = \rho(\phi, p)\ in\ \Omega
$$
with mole fractions of species $\alpha$
$$
y_\alpha(\phi, p)=y^E_\alpha exp(-\frac{z_\alpha e(\phi-\phi^E)+v_\alpha(p-p^E)}{k_B T})
$$
and charge density
$$
\rho(\phi, p) = e \sum_{\alpha=0}^N z_\alpha \frac{y_\alpha (\phi, p)}{\sum_{\alpha=0}^N v_\alpha y_\alpha (\phi, p)}
$$
We must also solve for pressure to establish mechanical equilibrium with a simple pressure poisson.
$$
\nabla p = -q\nabla\phi \longrightarrow \Delta p = -\nabla \cdot (q\nabla\phi) \ in\ \Omega
$$

After simulating we get $q(V)$, $C(V)$, $c_\alpha(x)$, $\phi(x)$.

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
├── results                      # Subdirectory with calculation results
├── notebooks
│   ├── demo-notebook.jl         # Notebook template
│   ├── ICMPB.jl                 # Ion conserving MPB draft with different solver variants
│   ├── ICMPB-provide-csv.jl     # Ion conserving MPB solution with csv output of concentration
│   ├── pluto-equilibrium.jl     # MPB solver comparing pressure poisson with incompressibility
│   ├── pluto-simplecell-bsk.jl  # MPB demo (WIP) with Bazant-Storey-Kornyshev model
│   └── pluto-simplecell.jl      # Poisson-Boltzmann/Bikermann demo
├── src
│   └── PoissonBoltzmannIPAM2025.jl  # Module with shared code for project
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

## Workflow hints
For workflow hints, see
- https://modernjuliaworkflows.org/
- https://j-fu.github.io/marginalia/julia/
  
