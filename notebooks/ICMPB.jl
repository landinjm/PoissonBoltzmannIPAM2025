### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ a70cef7d-2a2f-4155-bdf3-fec9df94c63f
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    using PlutoUI, HypertextLiteral, UUIDs
    using LinearAlgebra
    using Interpolations
    using VoronoiFVM, ExtendableGrids
    using LaTeXStrings
    using LessUnitful, Unitful
    using PreallocationTools
    using LaTeXStrings
    using DoubleFloats
    using ForwardDiff
    using PythonPlot
    using NLsolve
end

# ╔═╡ 63edffa3-b937-4f60-ae9f-13fa02b1a8f5
TableOfContents(aside = true)

# ╔═╡ 47ab15c5-8805-40a2-a5e8-48f5b51045a5
md"""
# ICMPB: Ion conserving modified Poisson-Boltzmann
"""

# ╔═╡ 62e5cd7f-b929-4b44-8f47-a654c7ddeca5
md"""
- For the term, see e.g. "Ion-Conserving Modified Poisson–Boltzmann Theory Considering a Steric Effect in an Electrolyte", Journal of the Physical Society of Japan, December 15, 2016, Vol. 85, No. 12, DOI [10.7566/JPSJ.85.124006](https://doi.org/10.7566/JPSJ.85.124006)
- According to Google Scholar, there have been no studies so far of the situation with surface charges instead of applied voltages.
- This notebook studies the special case with equal ion sizes which results in a constant summary concentraion ``c̄=\sum_{i=0}^{N} c_i`` where ``c_0`` is the solvent concentation.
"""

# ╔═╡ ea15a97a-45d5-4076-989d-39024f532997
md"""
## MPB -- Modified Poisson-Boltzmann
"""

# ╔═╡ 3befe295-89cc-4da9-9dce-32563a22de81
md"""

Let ``\bar c`` be the summary concentration, and ``y_i`` be the solute molar fractions, such that for the species molar concentrations one has ``c_i=\bar c y_i``. Then the Poisson equation for ``N`` charged species states: 
```math
\begin{aligned}
      -\nabla ⋅ ((1+χ)ε_0 \nabla \phi) &=ρ= F\sum_{i=1}^n z_i \bar c y_i \quad\text(MPB)
\end{aligned}
```

The mole fractions ``y_i`` are calculated as

```math
\begin{aligned}
   y_i&= \frac{\frac{c_{i,ref}}{c_{0,ref}}\exp\left(z_i\phi \frac{F}{RT}\right)}{1+ f_{mod}\sum_{j=1}^n \frac{c_{j,ref}}{c_{0,ref}}\exp\left(z_i\phi \frac{F}{RT}\right)}
\end{aligned}
```
where ``f_{mod}=1`` indicates the Bikerman (modified Poisson Boltzmann) model and  ``f_{mod}=0`` leads to the unmodified Poisson-Boltzmann equation.
"""


# ╔═╡ e118211d-c341-47ad-8693-6e5fe5af2fff
md"""
We assume the the problem to be given in a domain ``Ω=(0,L)`` with boundary conditons
```math
\begin{aligned}
	(1+χ)\nabla \phi|_{z=0} &= -q\\
	(1+χ)\nabla \phi|_{z=L} &= q
\end{aligned}
```

Further parameters of the system are the reference concentrations ``c_{i,ref}`` ``(i=1\dots N)`` which allow to calculate the solvent reference concentration
``c_{0,ref}=c̄ - \sum_{i=1}^N c_{i,ref}``.
"""

# ╔═╡ 45503458-cd99-448b-b5ff-529bb858186f
md"""
The average  ion concentrations are given as 
```math
c_{i,avg}=\frac1{|\Omega|}\int_\Omega c_i d\omega
```
"""

# ╔═╡ 003cf25f-8930-43bb-82ab-d13933d12ecf
md"""
### Field dependent dielectric decrement

```math
χ=χ(E)= \frac{χ_S}{\sqrt{1+ a|E|^2}}
```
"""

# ╔═╡ 1266d252-56c0-4b19-8ac8-439eb6943866
md"""
## ICMPB: Ion conserving modified Poisson-Boltzmann
"""

# ╔═╡ 08eb4212-ae3e-4471-a22b-d1df5ed9a544
md"""
Then, the ion conserving modified Poisson-Boltzmann problem consists in solving
[MPB](#MPB-–-Modified-Poisson-Boltzmann) for ``ϕ`` and ``c_{i,ref}\; (i=1\dots N)`` for given values of ``q`` and ``c_{i,avg}`` with ``\sum_{i=1}^N z_ic_{i,avg}=0``.
"""

# ╔═╡ aa497ce6-2d0a-4b9a-8053-197136f97974
md"""
First numerical experiments suggested that this formulation misses uniqueness of the solution. Indeed, at least for ``f_{mod}=0`` it can be shown that for a given solution ``\big(ϕ, \{c_{i,ref}\}_{i=1\dots N}\big)``, and some ``δ``,  ``\big(ϕ+\delta, \{c^\delta_{i,ref}\}_{i=1\dots N}\big)`` is another solution, where ``c_{i,ref}^δ`` can be obtained from ``c_{i,ref}`` via a solving
a linear system of equations. 

"""

# ╔═╡ de71225f-e80d-4e03-b30a-78165475dec5
md"""
As a consequence, we need to add another condition to  setting of the ICMPB problem. There are several options, including
- Antisymmetry of the solution: ``ϕ(0)+ϕ(L)=0``
- Zero potential at domain centert: ``ϕ(\frac{L}{2}=0``
- Electroneutrality of the reference concentrations: ``\sum_{i=1}^N z_i c_{i,ref} = 0``
While all of them have been tested to work,  electroneutrality appears to be the most reasonable.
"""

# ╔═╡ 7ddbf1d9-1551-4e0e-a830-01f13dcc5bd5
md"""
## Implementation
"""

# ╔═╡ 9ba8a2a3-db8b-405b-92e1-38908b486f70
md"""
### Constants and units
"""

# ╔═╡ 67de5632-0645-41f0-a9fd-610e40bfc456
begin
    const χ_S = 78.49 - 1
    const F = ph"N_A" * ph"e"
    const K = ufac"K"
    const nm = ufac"nm"
    const m = ufac"m"
    const dm = ufac"dm"
    const V = ufac"V"
    const mol = ufac"mol"
    const T = (273.15 + 25) * ufac"K"
    const RT = ph"R" * T
    const c̄ = 55.508mol / dm^3 # summary molar concentration
    const ε_0 = ph"ε_0"
end

# ╔═╡ 832362ed-8e40-4c58-af87-9e72a80bb420
md"""
### Parameters
"""

# ╔═╡ 4cb3283e-d2b4-41e6-ab57-083ec02a6b14
begin
    const f_mod = true # model choice: 0: Boltzmann, 1: Bikerman
    const nref = 4 # grid refinement level
    const L = 10nm # computational domain size
    const n_e = 10 # number of electrons/nm^2 at interfaces, defines q
    const M_avg = 10 # average molarity
    const E_0 = 10V / nm # decrement parameter
    const a = 5.0 / E_0^2 # decrement parameter in χ(E)
    const z = [-1, 1.0] # species charge numbers
end;

# ╔═╡ b8bf15c2-a822-4294-a877-490a918103e6
const N = length(z)

# ╔═╡ e74073e5-62af-4cb3-8c03-134700d16714
const c_avg = fill(M_avg * mol / dm^3, N)

# ╔═╡ 6f08912b-f401-4d4b-96f1-a6404ca4027f
surfcharge(n) = n * ph"e" / ufac"nm^2"

# ╔═╡ cf2b2ef6-e407-4424-9004-aba5b6855bdd
const q = surfcharge(n_e) # surface charge

# ╔═╡ 5bbab594-84c0-4c61-8f1f-e55dac40a328
const c0_avg = c̄ - sum(c_avg) # solvent bulk molar concentration

# ╔═╡ d65b73d4-28f5-4cf9-9491-9fb545afe2ba
const l_debye = sqrt((1 + χ_S) * ε_0 * RT / (F^2 * c_avg[1])) |> u"nm"

# ╔═╡ 869a0bba-3cc8-459f-8dcd-d70c111df4d4
md"""
### Discretization grid
"""

# ╔═╡ ab73c7cd-d20b-42dd-a89e-e5258b8d2bd0
begin
    const hmin = 1.0e-1 * nm * 2.0^(-nref) # grid size at working electrode
    const hmax = 1.0 * nm * 2.0^(-nref) # grid size at bulk

    X0 = geomspace(0, L / 2, hmin, hmax)
    X1 = geomspace(L / 2, L, hmax, hmin)
    grid = simplexgrid(glue(X0, X1))
    bfacemask!(grid, [L / 2], [L / 2], 3, tol = 1.0e-10 * nm)
end

# ╔═╡ 1ad518f4-9968-4dcf-b566-ed6fd70f0474
const i3 = grid[BFaceNodes][1, 3] # Index  of grid midpoint

# ╔═╡ 2f189782-a7b4-4e32-86f8-ad42e3a5250f
md"""
### Data & callbacks
"""

# ╔═╡ 099c8ab7-1432-4ff1-9fe8-70e085430cfa
Base.@kwdef mutable struct PBData
    c_ref::Vector{Float64} = c_avg
    c_avg::Vector{Float64} = c_avg
    c̄::Float64 = c̄
    c0_ref::Float64 = c̄ - sum(c_ref)
    q::Float64 = q
    z::Vector{Float64} = z
    N::Int = length(z)
    a::Float64 = a
    E_0::Float64 = E_0
    f_mod::Bool = f_mod
    F::Float64 = F
    RT::Float64 = RT
    ε_0::Float64 = ε_0
    χ_S::Float64 = χ_S
    iϕ::Int = 1
    cache = DiffCache(ones(Float64, length(z)))
end

# ╔═╡ e3123f2a-9607-4b93-8568-dde3091a40ce
function molfractions!(y, ϕ, c_ref, c0_ref, data)
    (; z, N, f_mod, F, RT) = data
    for ic in 1:N
        y[ic] = exp(-z[ic] * ϕ * F / RT) * c_ref[ic] / c0_ref
    end
    denom = 1.0 / (one(ϕ) + f_mod * sum(y))
    for ic in 1:N
        y[ic] = y[ic] * denom
    end
    return nothing
end

# ╔═╡ 03ac2f7e-1d6b-419b-a4cb-9e4a531ab9f1
function concentrations(sol, data; c_ref = data.c_ref)
    (; c̄, iϕ, N) = data
    n = size(sol, 2)
    c = zeros(N, n)
    y = zeros(N)
    c0_ref = data.c̄ - sum(c_ref)
    for iz in 1:n
        molfractions!(y, sol[data.iϕ, iz], c_ref, c0_ref, data)
        for ic in 1:N
            c[ic, iz] = c̄ * y[ic]
        end
    end
    return c
end

# ╔═╡ d3bcbe56-761d-4897-9581-ffb71bba0a16
function avgconcentrations(sol, sys; data = data(sys), c_ref = data.c_ref)
    c = concentrations(sol, data; c_ref)
    (; N = data)
    cavg = zeros(N)
    nv = nodevolumes(sys)
    L = sum(nv)
    for ic in 1:N
        cavg[ic] = nv ⋅ c[ic, :] / L
    end
    return cavg
end

# ╔═╡ 0c914971-48a0-4fee-996f-90bb156f8552
function flux!(y, u, edge, data)
    (; χ_S, a, ε_0, iϕ) = data
    eins = one(eltype(u))
    h = edgelength(edge)
    E = (u[iϕ, 1] - u[iϕ, 2]) / h
    χ = χ_S / sqrt(eins + a * E^2)
    ε = (eins + χ) * ε_0
    y[iϕ] = ε * ((u[iϕ, 1] - u[iϕ, 2]))
    return nothing
end

# ╔═╡ 8adfa102-417a-4106-a387-0cc16788460a
function relpermittivity(sol, data; grid)
    (; χ_S, a, iϕ) = data
    X = grid[Coordinates][1, :]
    return χ_S ./ (a * ((sol[iϕ, 2:end] - sol[iϕ, 1:(end - 1)]) ./ (X[2:end] - X[1:(end - 1)])) .^ 2 .+ 1) .+ 1

end

# ╔═╡ b5ed07a7-d5a8-4eaf-9592-7a3445e06d26
function spacecharge!(y, ϕ, c_ref, c0_ref, data)
    (; N, F, c_ref, c0_ref, c̄, z) = data
    molfractions!(y, ϕ, c_ref, c0_ref, data)
    sumyz = zero(ϕ)
    for i in 1:N
        sumyz += z[i] * y[i]
    end
    return F * c̄ * sumyz
end

# ╔═╡ 041890d8-491c-4ca4-a80a-3e31e5b7a8ba
function reaction!(y, u, node, data)
    (; cache, c_ref, c0_ref, iϕ) = data
    tmp = get_tmp(data.cache, u)
    y[iϕ] = -spacecharge!(tmp, u[iϕ], c_ref, c0_ref, data)
    return nothing
end

# ╔═╡ 68570541-90b8-4e5a-90a2-322b65ecd510
function bcondition!(y, u, bnode, data)
    boundary_neumann!(y, u, bnode, species = data.iϕ, region = 2, value = -data.q)
    boundary_neumann!(y, u, bnode, species = data.iϕ, region = 1, value = data.q)
    return nothing
end

# ╔═╡ 111388ba-674e-453e-a3ea-a7b482d1e9d8
md"""
## MPB system with given c_ref
"""

# ╔═╡ fbc93d62-6472-47a0-b7d1-e666fe17f341
md"""
Given ``q``, ``c_{ref}`` solve  ([MPB](#MPB-–-Modified-Poisson-Boltzmann)).
"""

# ╔═╡ 5b34f6d3-f4d4-4f2a-8a0a-f9772f7818cc
function MPBSystem(grid, data)
    return VoronoiFVM.System(
        grid;
        data,
        reaction = reaction!,
        flux = flux!,
        bcondition = bcondition!,
        species = [1]
    )
end

# ╔═╡ 2bb0f3a9-abdf-4368-8ce4-2af9a8ad14d9
md"""
``M_{ref}``: 
$(@bind M1_ref confirm(PlutoUI.Slider(0.1:0.1:10, default=1, show_value=true),label="Go"))
``n_e``: $(@bind n1_e confirm(PlutoUI.Slider(1:1:20, default=10, show_value=true),label="Go"))
"""

# ╔═╡ e684c78f-4663-46ae-abd0-002a58985461
data1 = PBData(c_ref = fill(M1_ref * mol / dm^3, 2), q = surfcharge(n1_e));

# ╔═╡ 6e7dec28-af3f-48f4-bc9f-63852504d186
sys1 = MPBSystem(grid, data1);

# ╔═╡ 0f05cb3b-50f6-4298-b3c8-19dda854b203
sol1 = solve(sys1);

# ╔═╡ 56b5138c-6d79-4957-9bab-232c8cd10fa3
md"""
## ICMPB solution using outer iteration
"""

# ╔═╡ 1bd01480-309a-44bc-b130-8a5be25e71ec
function extcref(cref0, data)
    (; z, N) = data
    return push!(copy(cref0), -z[1:(N - 1)] ⋅ cref0 / z[N])
end

# ╔═╡ 3c6728cc-28ed-4829-a5f2-f0eed91282c9
lastsol = unknowns(sys1, inival = 0)

# ╔═╡ b7ccc508-71d2-4f84-be64-4de290c1e2c5
function clonedata(data0, c_ref)
    data = deepcopy(data0)
    data.c_ref .= c_ref
    data.c0_ref = c̄ - sum(data.c_ref)
    return data
end

# ╔═╡ 4589c222-126a-4995-b637-baf6ad68ba5d
cb = (0.01:0.01:1) * mol / dm^3

# ╔═╡ 770b30d3-0cec-41c7-bad8-4245a8e0631b
md"""
``M_{avg}``: 
$(@bind M2_avg confirm(PlutoUI.Slider(0.1:0.1:10, default=5, show_value=true),label="Go"))
``n_e``: $(@bind n2_e confirm(PlutoUI.Slider(1:1:20, default=10, show_value=true),label="Go"))
"""

# ╔═╡ 39554f93-18b7-4b87-a803-2777b2277b32
function cavg(c_ref; q = surfcharge(n2_e), sys = sys1, kwargs...)
    data = deepcopy(VoronoiFVM.data(sys))
    data = clonedata(data, extcref(c_ref, data))
    data.q = q
    (; c̄) = data
    state = VoronoiFVM.SystemState(sys; data)
    sol = solve!(state; inival = 0, kwargs...)
    lastsol .= sol
    ca = avgconcentrations(sol, sys; data)
    return [ca[1]], sol
end

# ╔═╡ 75e46b6c-0736-42af-990f-d5a7fc9fdd95
ca = [ cavg([c]; q = surfcharge(n_e), sys = sys1, verbose = "", damp_initial = 0.1)[1][1]  for c in cb]

# ╔═╡ 8bafe67d-ffd1-4b2b-90d5-0e7e364ce050
let
    PythonPlot.clf()
    fig, ax = pyplot.subplots(1, 1)
    fig.set_size_inches(8, 2)
    ax.grid()
    ax.plot(cb, ca)
    PythonPlot.gcf()
end

# ╔═╡ 8fd77ff6-ed4c-4c02-bc57-d6ec91b9fba2
c2_avg = [M2_avg * mol / dm^3]

# ╔═╡ 7cd7090b-08fe-4588-8179-f86e09f874b5
res = nlsolve(
    c_ref -> cavg(
        c_ref;
        verbose = "",
        damp_initial = 0.1
    )[1] - c2_avg,
    c2_avg * 0.1,
    ftol = 1.0e-14
)

# ╔═╡ fee4a70b-4da4-45ff-b3af-ef1177113e3f
c2_ref = extcref(res.zero, VoronoiFVM.data(sys1))

# ╔═╡ 16c77b3f-84ee-46ad-b828-7fb60d5eba57
sol2 = cavg(res.zero; sys = sys1, verbose = "")[2];

# ╔═╡ d84b0ad9-0c3c-48b2-b5b7-cfc7d532bdd9
md"""
## ICMPB solution using extended system 
"""

# ╔═╡ 258ff2ba-c6c1-4d59-8d46-64cecc1059e9
function ICMPBSystem(; data = PBData(), generic = VoronoiFVM.nofunc)
    sys = VoronoiFVM.System(
        grid;
        data,
        generic,
        flux = flux!,
        bcondition = bcondition!,
        unknown_storage = :sparse
    )
    enable_species!(sys, data.iϕ, [1])
    for ic in 1:(data.N - 1)
        enable_boundary_species!(sys, data.iϕ + ic, [3])
    end
    return sys
end

# ╔═╡ f7a0fd72-5493-44c0-ab73-e7639bb1cb8b
sys3_0 = ICMPBSystem();

# ╔═╡ fedf30ad-df53-4571-a727-3913e4fb48e1
const nv = nodevolumes(sys3_0);

# ╔═╡ a3a3b816-d7e7-4942-b92f-7028b20c4020
const idx = unknown_indices(unknowns(sys3_0));

# ╔═╡ ca58971d-515e-491b-b7b8-2ebed16409c4
idx[2, i3]

# ╔═╡ 4445180b-3e7b-4315-8d6e-37c02a9886eb
md"""
``M_{avg}``: 
$(@bind M3_avg confirm(PlutoUI.Slider(0.1:0.1:10, default=5, show_value=true),label="Go"))
``n_e``: $(@bind n3_e confirm(PlutoUI.Slider(1:1:20, default=10, show_value=true),label="Go"))
"""

# ╔═╡ 30a93ebe-43a6-4ce9-afe5-de8cb3cc9a5d
c3_avg = fill(M3_avg * mol / dm^3, 2)

# ╔═╡ 633ec843-2aaa-46b6-beda-cd1f2e0ab430
data3 = PBData(c_avg = c3_avg, c_ref = 0.5 * c3_avg, q = surfcharge(n3_e))

# ╔═╡ 37581680-dcfd-46f7-abe3-7b64338a5c07
function xreaction!(f, u, sys)
    data = data3 #VoronoiFVM.data(sys)
    (; cache, iϕ, N, z, F, c_avg, c̄) = data
    y = get_tmp(cache, u)
    c_ref = [u[idx[ic + iϕ, i3]] for ic in 1:(N - 1)]
    push!(c_ref, -z[1:(N - 1)] ⋅ c_ref / z[N])
    c0_ref = c̄ - sum(c_ref)
    L = sum(nv)

    for ic in 1:(N - 1)
        f[idx[ic + iϕ, i3]] = 0
    end

    for iv in 1:length(nv)
        molfractions!(y, u[idx[iϕ, iv]], c_ref, c0_ref, data)
        f[idx[iϕ, iv]] = 0
        for ic in 1:N
            f[idx[iϕ, iv]] -= y[ic] * c̄ * nv[iv] * z[ic] * F
        end
        for ic in 1:(N - 1)
            f[idx[ic + iϕ, i3]] += y[ic] * c̄ * nv[iv]
        end
    end
    for ic in 1:(N - 1)
        f[idx[ic + iϕ, i3]] = f[idx[ic + iϕ, i3]] - c_avg[ic] * L
    end
    return nothing
end


# ╔═╡ e5b3dd26-7df7-41cf-ad5f-887d95f73135
sys3 = ICMPBSystem(data = data3, generic = xreaction!)

# ╔═╡ 64f4c78a-22b3-4c46-a4ab-59cf0b511756
state3 = VoronoiFVM.SystemState(sys3; data = data3)

# ╔═╡ a72ee1be-e1ee-49de-a5d5-4a7699236cfa
state3.matrix

# ╔═╡ 73791e32-0cb7-4feb-93b1-2933ac662a0c
begin
    inival3 = unknowns(sys3, inival = 0.0)
    inival3[2:N, i3] .= c3_avg[1:(N - 1)] / 2
end;

# ╔═╡ 7d5cf0bf-515f-46a7-9970-742911e27974
sol3 = solve!(state3; inival = inival3, verbose = "n", damp_initial = 0.5)

# ╔═╡ 8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
html"""<hr>"""

# ╔═╡ 2e8d2fcf-c623-4142-a06e-a1a3bd02bf30
md"""
## Service functions
"""

# ╔═╡ fb1c3580-149a-41d6-97eb-e3ce76be331b
myround(x) = round(x, sigdigits = 4)

# ╔═╡ 3fa9d05a-a3f6-4e07-9a5d-f53ee2126cae
function plotsol(sol, sys; data = data(sys), grid = grid, c_ref = data.c_ref, size = (600, 400))
    PythonPlot.clf()
    fig, ax = pyplot.subplots(2, 1)
    fig.set_size_inches(size[1] / 100, size[2] / 100)
    ax1 = ax[0]
    ax2 = ax[1]
    ax1.grid()
    ax1r = ax1.twinx()
    X = grid[Coordinates][1, :]
    ε_r = relpermittivity(sol, data; grid)
    c = concentrations(sol, data; c_ref)
    c0 = -(sum(c, dims = 1) .- c̄)
    ax1.set_title("ϕ∈$(round.(Float64.(extrema(sol[1, :])), sigdigits = 3)), ε_r ∈$(round.(Float64.(extrema(ε_r)), sigdigits = 3))")
    ax1.plot(X / nm, sol[1, :], color = "green", linewidth = 2, label = "ϕ")
    ax1r.plot(X[1:(end - 1)] / nm, ε_r, color = "pink", linewidth = 3, label = L"ε_r")
    #   ax1.set_ylim(-10, 10)
    ax1.set_xlabel("z/nm")
    ax1.set_ylabel("ϕ/V")
    ax1.legend(loc = (0.1, 0.1))
    ax1r.legend(loc = (0.8, 0.1))
    ax1r.set_ylim(0, 80)

    ax2.grid()
    cavg = avgconcentrations(sol, sys; data, c_ref)
    cm, cp = cavg[1] / (mol / dm^3), cavg[2] / (mol / dm^3)
    crm, crp = c_ref[1] / (mol / dm^3), c_ref[2] / (mol / dm^3)


    ax2.set_title("M_ref=$(myround.((crm, crp))),  M_avg=$(myround.((cm, cp)))")
    ax2.set_xlabel("z/nm")
    ax2.set_ylabel("c/(mol/L)")
    ax2.set_ylim(0, 60)

    ax2.plot(X / nm, c[1, :] / (mol / dm^3), color = "blue", linewidth = 2, label = L"c^-")
    ax2.plot(X / nm, c[2, :] / (mol / dm^3), color = "red", linewidth = 2, label = L"c^+")
    ax2.plot(X / nm, c0[1, :] / (mol / dm^3), color = "green", linewidth = 2, label = L"c_{solvent}")
    ax2.legend(loc = (0.4, 0.1))


    tight_layout()
    return PythonPlot.gcf()
end

# ╔═╡ 97a6fd15-9e27-49ed-94f5-4b881d7a907c
plotsol(sol1, sys1; size = (600, 400))

# ╔═╡ c5a3e303-98a2-46a6-8cd9-8953248c5899
plotsol(sol2, sys1, data = data(sys1); c_ref = c2_ref, size = (600, 400))

# ╔═╡ 13cb4936-ff78-4064-85c5-a7d51c4963e3
plotsol(
    sol3, sys3;
    data = data3,
    c_ref = extcref(sol3[2:N, i3], data3),
    size = (600, 400)
)

# ╔═╡ 784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
html"""<style>.dont-panic{ display: none }</style>"""

# ╔═╡ afe4745f-f9f1-4e23-8735-cbec6fb79c41
begin
    function floataside(text::Markdown.MD; top = 1)
        uuid = uuid1()
        return @htl(
            """
            		<style>


            		@media (min-width: calc(700px + 30px + 300px)) {
            			aside.plutoui-aside-wrapper-$(uuid) {

            	color: var(--pluto-output-color);
            	position:fixed;
            	right: 1rem;
            	top: $(top)px;
            	width: 400px;
            	padding: 10px;
            	border: 3px solid rgba(0, 0, 0, 0.15);
            	border-radius: 10px;
            	box-shadow: 0 0 11px 0px #00000010;
            	/* That is, viewport minus top minus Live Docs */
            	max-height: calc(100vh - 5rem - 56px);
            	overflow: auto;
            	z-index: 40;
            	background-color: var(--main-bg-color);
            	transition: transform 300ms cubic-bezier(0.18, 0.89, 0.45, 1.12);

            			}
            			aside.plutoui-aside-wrapper > div {
            #				width: 300px;
            			}
            		}
            		</style>

            		<aside class="plutoui-aside-wrapper-$(uuid)">
            		<div>
            		$(text)
            		</div>
            		</aside>

            		"""
        )
    end
    floataside(stuff; kwargs...) = floataside(md"""$(stuff)"""; kwargs...)
end;


# ╔═╡ b8fd36a7-d8d1-45f7-b66e-df9132168bfc
# https://discourse.julialang.org/t/adding-a-restart-process-button-in-pluto/76812/5
restart_button() = html"""
<script>
	const button = document.createElement("button")

	button.addEventListener("click", () => {
		editor_state_set(old_state => ({
			notebook: {
				...old_state.notebook,
				process_status: "no_process",
			},
		})).then(() => {
			window.requestAnimationFrame(() => {
				document.querySelector("#process_status a").click()
			})
		})
	})
	button.innerText = "Restart notebook"

	return button
</script>
""";

# ╔═╡ Cell order:
# ╠═63edffa3-b937-4f60-ae9f-13fa02b1a8f5
# ╠═a70cef7d-2a2f-4155-bdf3-fec9df94c63f
# ╟─47ab15c5-8805-40a2-a5e8-48f5b51045a5
# ╟─62e5cd7f-b929-4b44-8f47-a654c7ddeca5
# ╟─ea15a97a-45d5-4076-989d-39024f532997
# ╟─3befe295-89cc-4da9-9dce-32563a22de81
# ╟─e118211d-c341-47ad-8693-6e5fe5af2fff
# ╟─45503458-cd99-448b-b5ff-529bb858186f
# ╟─003cf25f-8930-43bb-82ab-d13933d12ecf
# ╟─1266d252-56c0-4b19-8ac8-439eb6943866
# ╟─08eb4212-ae3e-4471-a22b-d1df5ed9a544
# ╟─aa497ce6-2d0a-4b9a-8053-197136f97974
# ╟─de71225f-e80d-4e03-b30a-78165475dec5
# ╟─7ddbf1d9-1551-4e0e-a830-01f13dcc5bd5
# ╟─9ba8a2a3-db8b-405b-92e1-38908b486f70
# ╠═67de5632-0645-41f0-a9fd-610e40bfc456
# ╟─832362ed-8e40-4c58-af87-9e72a80bb420
# ╠═4cb3283e-d2b4-41e6-ab57-083ec02a6b14
# ╠═b8bf15c2-a822-4294-a877-490a918103e6
# ╠═e74073e5-62af-4cb3-8c03-134700d16714
# ╠═6f08912b-f401-4d4b-96f1-a6404ca4027f
# ╠═cf2b2ef6-e407-4424-9004-aba5b6855bdd
# ╠═5bbab594-84c0-4c61-8f1f-e55dac40a328
# ╠═d65b73d4-28f5-4cf9-9491-9fb545afe2ba
# ╟─869a0bba-3cc8-459f-8dcd-d70c111df4d4
# ╠═ab73c7cd-d20b-42dd-a89e-e5258b8d2bd0
# ╠═1ad518f4-9968-4dcf-b566-ed6fd70f0474
# ╟─2f189782-a7b4-4e32-86f8-ad42e3a5250f
# ╠═099c8ab7-1432-4ff1-9fe8-70e085430cfa
# ╠═e3123f2a-9607-4b93-8568-dde3091a40ce
# ╠═03ac2f7e-1d6b-419b-a4cb-9e4a531ab9f1
# ╠═d3bcbe56-761d-4897-9581-ffb71bba0a16
# ╠═0c914971-48a0-4fee-996f-90bb156f8552
# ╠═8adfa102-417a-4106-a387-0cc16788460a
# ╠═b5ed07a7-d5a8-4eaf-9592-7a3445e06d26
# ╠═041890d8-491c-4ca4-a80a-3e31e5b7a8ba
# ╠═68570541-90b8-4e5a-90a2-322b65ecd510
# ╟─111388ba-674e-453e-a3ea-a7b482d1e9d8
# ╟─fbc93d62-6472-47a0-b7d1-e666fe17f341
# ╠═5b34f6d3-f4d4-4f2a-8a0a-f9772f7818cc
# ╟─2bb0f3a9-abdf-4368-8ce4-2af9a8ad14d9
# ╠═97a6fd15-9e27-49ed-94f5-4b881d7a907c
# ╠═e684c78f-4663-46ae-abd0-002a58985461
# ╠═6e7dec28-af3f-48f4-bc9f-63852504d186
# ╠═0f05cb3b-50f6-4298-b3c8-19dda854b203
# ╟─56b5138c-6d79-4957-9bab-232c8cd10fa3
# ╠═1bd01480-309a-44bc-b130-8a5be25e71ec
# ╠═3c6728cc-28ed-4829-a5f2-f0eed91282c9
# ╠═b7ccc508-71d2-4f84-be64-4de290c1e2c5
# ╠═39554f93-18b7-4b87-a803-2777b2277b32
# ╠═4589c222-126a-4995-b637-baf6ad68ba5d
# ╠═75e46b6c-0736-42af-990f-d5a7fc9fdd95
# ╠═8bafe67d-ffd1-4b2b-90d5-0e7e364ce050
# ╟─770b30d3-0cec-41c7-bad8-4245a8e0631b
# ╠═c5a3e303-98a2-46a6-8cd9-8953248c5899
# ╠═8fd77ff6-ed4c-4c02-bc57-d6ec91b9fba2
# ╠═7cd7090b-08fe-4588-8179-f86e09f874b5
# ╠═fee4a70b-4da4-45ff-b3af-ef1177113e3f
# ╠═16c77b3f-84ee-46ad-b828-7fb60d5eba57
# ╟─d84b0ad9-0c3c-48b2-b5b7-cfc7d532bdd9
# ╠═258ff2ba-c6c1-4d59-8d46-64cecc1059e9
# ╠═f7a0fd72-5493-44c0-ab73-e7639bb1cb8b
# ╠═fedf30ad-df53-4571-a727-3913e4fb48e1
# ╠═a3a3b816-d7e7-4942-b92f-7028b20c4020
# ╠═ca58971d-515e-491b-b7b8-2ebed16409c4
# ╠═37581680-dcfd-46f7-abe3-7b64338a5c07
# ╟─4445180b-3e7b-4315-8d6e-37c02a9886eb
# ╠═13cb4936-ff78-4064-85c5-a7d51c4963e3
# ╠═30a93ebe-43a6-4ce9-afe5-de8cb3cc9a5d
# ╠═633ec843-2aaa-46b6-beda-cd1f2e0ab430
# ╠═e5b3dd26-7df7-41cf-ad5f-887d95f73135
# ╠═64f4c78a-22b3-4c46-a4ab-59cf0b511756
# ╠═a72ee1be-e1ee-49de-a5d5-4a7699236cfa
# ╠═73791e32-0cb7-4feb-93b1-2933ac662a0c
# ╠═7d5cf0bf-515f-46a7-9970-742911e27974
# ╟─8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
# ╟─2e8d2fcf-c623-4142-a06e-a1a3bd02bf30
# ╠═fb1c3580-149a-41d6-97eb-e3ce76be331b
# ╠═3fa9d05a-a3f6-4e07-9a5d-f53ee2126cae
# ╟─784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
# ╟─afe4745f-f9f1-4e23-8735-cbec6fb79c41
# ╟─b8fd36a7-d8d1-45f7-b66e-df9132168bfc
