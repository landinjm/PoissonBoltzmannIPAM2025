### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
if isdefined(Main, :PlutoRunner)
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    using VoronoiFVM
    using ExtendableGrids
    using LinearAlgebra
    using LessUnitful
end

# ╔═╡ ef660f6f-9de3-4896-a65e-13c60df5de1e
md"""
## Ion conserving MPB solver
"""

# ╔═╡ 920b7d84-56c6-4958-aed9-fc67ba0c43f6
md"""
## 1. Intro

This code implements the model described in
[Müller, R., Fuhrmann, J., & Landstorfer, M. (2020). Modeling polycrystalline electrode-electrolyte interfaces: The differential capacitance. Journal of The Electrochemical Society, 167(10), 106512](https://iopscience.iop.org/article/10.1149/1945-7111/ab9cca/meta)


Equation numbers refer to the paper.

Concentrations are given in number densities.
"""

# ╔═╡ 87ac16f4-a4fc-4205-8fb9-e5459517e1b8
md"""
If not stated otherwise, all calculations and calculation results are in coherent SI units.
"""

# ╔═╡ 7d77ad32-3df6-4243-8bad-b8df4126e6ea
md"""
## 2. Model data
"""

# ╔═╡ 4cabef42-d9f9-43fe-988e-7b54462dc775
md"""
#### ICMPBData
"""

# ╔═╡ 0d825f88-cd67-4368-90b3-29f316b72e6e
begin
    """
        ICMPBData

    Data structure containing data for equilibrium calculations
    All data including molarity in SI basic units
    """
    Base.@kwdef mutable struct ICMPBData

        "Ion charge numbers."
        z::Vector{Int} = [-1, 1]

        "Number of ionic species"
        N::Int64 = length(z)

        "Ion solvation numbers"
        κ::Vector{Float64} = fill(10.0, N)

        "Bulk molarity"
        molarity::Float64 = 0.1 * ph"N_A" / ufac"dm^3"

        "Bulk ion number densities"
        n_E::Vector{Float64} = fill(molarity, N)

        "Average ion number densities"
        n_avg::Vector{Float64} = fill(molarity, N)

        "Surface charges"
        q::Vector{Float64} = [0, 0]


        "Solvent molarity"
        n0_ref::Float64 = 55.508 * ph"N_A" / ufac"dm^3"

        "Solvent molecule volume"
        v0::Float64 = 1 / n0_ref

        "Unsolvated ion volume"
        vu::Vector{Float64} = fill(1 / n0_ref, N)

        "Dielectric susceptibility"
        χ::Float64 = 78.49 - 1

        "Solvent molar fraction index"
        i0::Int = N + 1

        "Electric potential species index"
        iφ::Int = i0 + 1

        "Pressure species index"
        ip::Int = iφ + 1

        "Offset of n_E in species list"
        coffset::Int = ip

        "Reference pressure"
        p_ref::Float64 = 1.0e5 * ufac"Pa"

        "Pressure scaling nparameter"
        pscale::Float64 = 1.0 * ufac"GPa"

        "Concentration scaling parameter"
        cscale::Float64 = ph"N_A"

        "Reference voltage"
        E_ref::Float64 = 0.0 * ufac"V"

        "Temperature"
        T::Float64 = 298.15 * ufac"K"

        "Temperature times Boltzmann constant"
        kT::Float64 = ph"k_B" * T

        "Electron charge"
        e::Float64 = ph"e"

        "Vacuum permittivity"
        ε_0::Float64 = ph"ε_0"


        "Ion conservation"
        conserveions::Bool = false

        "node volumes" # we should be able to query this from the system
        nv::Vector{Float64} = Float64[]

    end
end

# ╔═╡ 30c6a176-935b-423f-9447-86f78746322f
md"""
#### debyelength(data)

```math
L_{Debye}=\sqrt{ \frac{(1+χ)ε_0k_BT}{e^2n_E}}
```
"""

# ╔═╡ a41c6c1f-ceb5-4590-a421-cae5078d167b
function L_Debye(data)
    return sqrt(
        (1 + data.χ) * data.ε_0 * ph"k_B" * data.T / (data.e^2 * data.n_E[1]),
    )
end;

# ╔═╡ f3049938-2637-401d-9411-4d7be07c19ca
md"""
#### set_molarity!(data,M)
"""

# ╔═╡ 5d6340c4-2ddd-429b-a60b-3de5570a7398
function set_molarity!(data::ICMPBData, M_E)
    n_E = M_E * ph"N_A" / ufac"dm^3"
    data.molarity = n_E
    data.n_E = fill(n_E, data.N)
    data.n_avg = fill(n_E, data.N)
    return data
end

# ╔═╡ a21545da-3b53-47af-b0c4-f253b37dc84f
md"""

#### dlcap0(data)
Double layer capacitance at $φ=0$
```math
C_{dl,0}=\sqrt{\frac{2(1+χ) ε_0e^2 n_E}{k_BT}}
```
"""

# ╔═╡ 1d22b09e-99c1-4026-9505-07bdffc98582
function dlcap0(data::ICMPBData)
    return sqrt(
        2 * (1 + data.χ) * ph"ε_0" * ph"e"^2 * data.n_E[1] / (ph"k_B" * data.T),
    )
end;

# ╔═╡ 5a210961-19fc-40be-a5f6-033a80f1414d
md"""
Check with Bard/Faulkner: the value must be $(22.8)μF/cm^2")
"""

# ╔═╡ fe704fb4-d07c-4591-b834-d6cf2f4f7075
# ╠═╡ skip_as_script = true
#=╠═╡
let
    data=ICMPBData()
    set_molarity!(data,0.01)
    data.χ=78.49-1
    cdl0=dlcap0(data)
    @assert cdl0 ≈ 22.84669184882525ufac"μF/cm^2"
end
  ╠═╡ =#

# ╔═╡ 5eca37ba-f858-45fb-a66a-3795327dfd18
md"""
## 3. Model equations
"""

# ╔═╡ a26cf11b-0ce1-4c1d-a64d-1917178ff676
md"""
### Mole fractions
Equilibrium expression for mole fractions (``α≥0``) (16)
```math
y_α(φ,p)=y_α^E\exp\left(\frac{-z_αe}{k_BT}(φ- φ^E)-\frac{v_α}{k_BT}(p-p^E)\right)
```
"""

# ╔═╡ cdd1d359-08fa-45a1-a857-e19f2adefcab
md"""
#### y_α(φ,p,α,data)

Ion molar fractions
"""

# ╔═╡ 188f67d8-2ae8-474c-8e58-68b8b4fde02e
function y_α(φ, p, α, data, ddata)
    η_φ = data.z[α] * data.e * (φ - data.E_ref)
    η_p = ddata.v[α] * (p * data.pscale - data.p_ref)
    return ddata.y_E[α] * exp(-(η_φ + η_p) / (data.kT))
end;

# ╔═╡ f70eed13-a6c2-4d54-9f30-113367afaf7d
md"""
#### y0(p,data)

Solvent molar fraction
"""

# ╔═╡ d7531d5f-fc2d-42b2-9cf9-6a737b0f0f8d
function y0(p, data, ddata)
    return ddata.y0_E * exp(-data.v0 * (p * data.pscale - data.p_ref) / (data.kT))
end;

# ╔═╡ f6f004a6-d71b-4813-a363-9f51dc37e42a
md"""
### Poisson equation
Poisson equation (32a)

```math
-∇⋅(1+χ)ε_0∇φ = q(φ,p)
```
"""

# ╔═╡ 3810cc88-07f1-4741-853f-331e71c87923
md"""
#### poisson_flux!(f,u,edge,data)

VoronoiFVM flux function for left hand side of Poisson equation
"""

# ╔═╡ 0e2d20a1-5f26-4263-9a91-3b40b2c2996a
function poisson_flux!(f, u, edge, data)
    return f[iφ] = (1.0 + data.χ) * data.ε_0 * (u[iφ, 1] - u[iφ, 2])
end;

# ╔═╡ 824c610b-6e5e-48a3-be37-19104f52d1d9
md"""
#### Space charge expression
"""

# ╔═╡ 2e11ce81-7d0d-498f-9ddd-7d4d836ab42f
md"""
Solvated ion volumes:
```math
	v_α=(1+κ_α)v_0
```
"""

# ╔═╡ b1e062c6-f245-4edc-aa02-871e2c776998
md"""
Incompressibility condition (14)
```math
\begin{aligned}
1&=∑\limits_α v_αn_α= n ∑\limits_α v_α y_α\\
n&=\frac1{∑\limits_α v_α y_α}
\end{aligned}
```
"""

# ╔═╡ c4cc940c-74aa-45f8-a2fa-6016d7c3c145
md"""
Space charge
```math
\begin{aligned}
q(φ,p)&=e∑\limits_α z_αn_α = ne∑\limits_α z_αy_α\\
      &=e\frac{∑\limits_α z_αy_α(\phi,p)}{∑\limits_α v_α y_α(\phi,p)}
\end{aligned}
```
"""

# ╔═╡ b07246b8-aec5-4161-8879-8cefb350aced
function spacecharge(u, data)
    (; iφ, ip, i0) = data
    y = u[i0]
    sumyz = zero(eltype(u))
    sumyv = data.v0 * y
    for α in 1:(data.N)
        y = u[α]
        sumyz += data.z[α] * y
        v = data.vu[α] + data.κ[α] * data.v0
        sumyv += v * y
    end
    return data.e * sumyz / sumyv
end

# ╔═╡ b41838bb-3d5b-499c-9eb5-137c252ae366
md"""
#### Sum of mole fractions
"""

# ╔═╡ a468f43a-aa20-45dc-9c21-77f5adf2d700
function local_ysum(u, data)
    (; iφ, ip, i0) = data

    sumy = u[i0]
    for α in 1:(data.N)
        sumy += u[α]
    end
    return sumy
end

# ╔═╡ 13fc2859-496e-4f6e-8b22-36d9d55768b8
md"""
#### Derived data

Update derived data in data record.

Calculate bulk mole fractions from incompressibiltiy:
```math
\begin{aligned}
∑\limits_αv_αn_α^E&=1\\
n_0^E&=\frac1{v_0}\left(1-∑\limits_{α>0}v_αn_α^E\right)\\
n^E&=\frac1{v_0}\left(1-∑\limits_{α>0}v_αn_α^E\right)+ ∑\limits_{α>0}n_α^E\\
   &=\frac1{v_0}\left(1-∑\limits_{α>0}(v_α-v_0)n_α^E\right)\\
   &=\frac1{v_0}\left(1-∑\limits_{α>0}((1+ κ_α)v_0-v_0)n_α^E\right)\\
   &=\frac1{v_0}\left(1-∑\limits_{α>0}κ_αv_0n_α^E\right)\\
   &=\frac1{v_0}-∑\limits_{α>0}κ_αn_α^E\\
y_α^E&=\frac{n_α^E}{n^E}
\end{aligned}
```
"""

# ╔═╡ 7c5e9686-9d66-4ff0-84a7-c7b22596ab57
begin
    struct DerivedData{T}
        "Effective ion volumes"
        v::Vector{Float64}
        "Bulk ion mole fractions"
        y_E::Vector{T}
        "Bulk solvent mole fraction"
        y0_E::T
    end

    function DerivedData(data::ICMPBData, n_E)
        (; κ, v0, vu, T) = data
        c0 = zero(eltype(n_E)) + 1 / v0
        barc = zero(eltype(n_E))
        v = vu + κ * v0
        N = length(κ)
        for α in 1:N
            barc += n_E[α]
            c0 -= n_E[α] * (1 + κ[α])
        end
        barc += c0
        y_E = n_E / barc
        y0_E = c0 / barc
        return DerivedData(v, y_E, y0_E)
    end

    function DerivedData(data::ICMPBData)
        return DerivedData(data, data.n_E)
    end

end

# ╔═╡ b1e333c0-cdaa-4242-b71d-b54ff71aef83
let
    data = ICMPBData()
    set_molarity!(data, 0.01)
    ddata = DerivedData(data)
    sumyz = 0.0
    sumyv = ddata.y0_E * data.v0
    sumy = ddata.y0_E
    for α in 1:data.N
        v = (1.0 + data.κ[α]) * data.v0
        sumyz += ddata.y_E[α] * data.z[α]
        sumyv += ddata.y_E[α] * v
        sumy += ddata.y_E[α]
    end
    @assert sumy ≈ 1.0
end

# ╔═╡ 55bd7b9a-a191-4a0b-9c6b-13733be5023e
md"""
#### c_num!(c,φ,p, data)
Calculate number concentration at discretization node
```math
	n_α=ny_α
```
"""

# ╔═╡ 3ceda3b1-bf1c-4126-b94f-2ee03e8dde99
function c_num!(c, y, data, ddata)
    (; i0) = data
    sumyv = data.v0 * y[i0]
    for α in 1:(data.N)
        c[α] = y[α]
        sumyv += y[α] * ddata.v[α]
    end
    return c ./= sumyv
end;

# ╔═╡ 97c5942c-8eb4-4b5c-8951-87ac0c9f396d
function c0_num(y, data, ddata)
    (; i0) = data

    y0 = y[i0]
    sumyv = data.v0 * y0
    for α in 1:(data.N)
        sumyv += y[α] * ddata.v[α]
    end
    return y0 / sumyv
end;

# ╔═╡ 0c54efd0-f279-4dc6-8b00-ba092dd13f44
md"""
#### calc_cnum(sol,sys)

Obtain ion number densities from system
"""

# ╔═╡ 800dfed8-9f29-4138-96f8-e8bf1f2f00e6
function calc_cnum(sol, sys)
    data = sys.physics.data
    i3 = sys.grid[BFaceNodes][3][1]
    if data.conserveions
        ddata = DerivedData(data, sol[(data.coffset + 1):end, i3])
    else
        ddata = DerivedData(data)
    end
    (; iφ, ip, N) = data
    grid = sys.grid
    nnodes = num_nodes(grid)
    conc = zeros(data.N, nnodes)
    for i in 1:nnodes
        @views c_num!(conc[:, i], sol[1:(N + 1), i], data, ddata)
    end
    return conc
end;

# ╔═╡ 24910762-7d56-446b-a758-d8e830fe9a09
function calc_c0num(sol, sys)
    data = sys.physics.data
    grid = sys.grid
    i3 = sys.grid[BFaceNodes][3][1]
    if data.conserveions
        ddata = DerivedData(data, sol[(data.coffset + 1):end, i3])
    else
        ddata = DerivedData(data)
    end
    (; iφ, ip, N) = data
    nnodes = num_nodes(grid)
    c0 = zeros(nnodes)
    for i in 1:nnodes
        @views c0[i] = c0_num(sol[1:(N + 1), i], data, ddata)
    end
    return c0
end;

# ╔═╡ 9fe3ca93-c051-426e-8b9a-cc59f59319ad
md"""
#### calc_cmol(sol,sys)

Obtain ion  molarities (molar densities in mol/L)  from system
"""

# ╔═╡ 2ee34d76-7238-46c2-94d1-a40d8b017af6
calc_cmol(sol, sys) = calc_cnum(sol, sys) / (ph"N_A" * ufac"mol/dm^3");

# ╔═╡ 79cc671b-ef6e-42da-8641-61e43f221cb1
calc_c0mol(sol, sys) = calc_c0num(sol, sys) / (ph"N_A" * ufac"mol/dm^3");

# ╔═╡ 7a607454-7b75-4313-920a-2dbdad258015
md"""
## 7. The pressure Poisson equation
"""

# ╔═╡ 9cb8324c-896f-40f8-baa8-b7d47a93e9f5
md"""
An alternative possibility to handle the pressure has been introduced in 

[J. Fuhrmann, “Comparison and numerical treatment of generalised Nernst–Planck models,” Computer Physics Communications, vol. 196, pp. 166–178, 2015.](https://dx.doi.org/10.1016/j.cpc.2015.06.004).

Starting with the momentum balance in mechanical equilibrium
```math
	\nabla p = -q\nabla \varphi
```
by taking the divergence on both sides of the equation, one derives the pressure Poisson problem
```math
\begin{aligned}
	-\Delta p &= \nabla\cdot q\nabla \varphi & \text{in}\; \Omega\\
      p&=p_{bulk} & \text{on}\; \Gamma_{bulk}\\
	(\nabla p + q\nabla \varphi)\cdot \vec n &=0 & \text{on}\; \partial\Omega\setminus\Gamma_{bulk}\\
\end{aligned}
```
"""

# ╔═╡ 003a5c0b-17c7-4407-ad23-21c0ac000fd4
md"""
The bulk Dirichlet boundary condition for the pressure is necessary to make the solution unique. It is reasonable to set the ``\varphi`` to a bulk value at ``\Gamma_{bulk}`` as well, and to calculate ``p_{bulk}`` from the molar fraction sum constraint.
"""

# ╔═╡ e1c13f1e-5b67-464b-967b-25e3a93e33d9
function reaction!(f, u, node, data)
    (; i0, iφ, ip, N) = data
    φ = u[iφ]
    p = u[ip]
    f[iφ] = -spacecharge(u, data)
    if !data.conserveions
        ddata = DerivedData(data)
        f[i0] = u[i0] - y0(p, data, ddata)
        for α in 1:N
            f[α] = u[α] - y_α(φ, p, α, data, ddata)
        end
    end
    return
end;

# ╔═╡ 64e47917-9c61-4d64-a6a1-c6e8c7b28c59
function poisson_and_p_flux!(f, u, edge, data)
    (; iφ, ip, N) = data
    f[iφ] = (1.0 + data.χ) * data.ε_0 * (u[iφ, 1] - u[iφ, 2])
    uu = zeros(eltype(u), N + 3)
    for i in 1:(N + 3)
        uu[i] = u[i, 1]
    end
    q1 = spacecharge(uu, data)
    for i in 1:(N + 3)
        uu[i] = u[i, 2]
    end
    q2 = spacecharge(uu, data)
    f[ip] =
        (u[ip, 1] - u[ip, 2]) + (u[iφ, 1] - u[iφ, 2]) * (q1 + q2) / (2 * data.pscale)
    return
end;

# ╔═╡ 743b9a7a-d6ac-4da0-8538-2045d965b547
function bcondition!(y, u, bnode, data)
    (; iφ, ip) = data
    boundary_neumann!(y, u, bnode, species = iφ, region = 2, value = data.q[2])
    boundary_neumann!(y, u, bnode, species = iφ, region = 1, value = data.q[1])
    boundary_dirichlet!(y, u, bnode, species = ip, region = 3, value = 0)
    return nothing
end

# ╔═╡ 0b646215-32db-4219-904b-f86f8861b46a
function apply_charge!(data::ICMPBData, q)
    data.q .= [q, - q]
    return data
end

# ╔═╡ 48670f54-d303-4c3a-a191-06e6592a2e0a
function ysum(sys, sol)
    data = sys.physics.data
    n = size(sol, 2)
    sumy = zeros(n)
    for i in 1:n
        sumy[i] = local_ysum(sol[:, i], data)
    end
    return sumy
end

# ╔═╡ b0a45e53-8b98-4e18-8b41-7f6d0bc1f76e
function VoronoiFVM.unknowns(sys, data::ICMPBData)
    (; i0, iφ, ip, coffset, N) = data
    u = unknowns(sys, inival = 0)
    i3 = sys.grid[BFaceNodes][3][1]
    for i in 1:N
        u[i, :] .= 0.1
        if data.conserveions
            u[i + coffset, i3] = data.n_E[i] / data.cscale
        end
    end
    u[i0, :] .= 1 - N * 0.1
    return u
end

# ╔═╡ dbccaa88-65d9-47ab-be78-83df64a6db24
function ionconservation!(f, u, sys, data)
    (; coffset, i0, iφ, ip, N, z, nv, n_avg) = data
    f .= 0
    i3 = sys.grid[BFaceNodes][3][1]
    idx = unknown_indices(unknowns(sys))
    n_E = [u[idx[coffset + i, i3]] * data.cscale for i in 1:N]
    ddata = DerivedData(data, n_E)
    y = zeros(eltype(u), N)
    for i in 1:num_nodes(sys.grid)
        f[idx[i0, i]] = u[idx[i0, i]] - y0(u[idx[ip, i]], data, ddata)
        for α in 1:N
            f[idx[α, i]] = u[idx[α, i]] - y_α(u[idx[iφ, i]], u[idx[ip, i]], α, data, ddata)
        end
    end
    really_conserve = true
    L = sum(nv)
    f[idx[coffset + N, i3]] = u[idx[coffset + N, i3]]
    for ic in 1:(N - 1)
        if really_conserve
            f[idx[ic + coffset, i3]] = -n_avg[ic] * L / data.cscale
        else
            f[idx[ic + coffset, i3]] = u[idx[ic + coffset, i3]] - data.n_E[ic] / data.cscale
        end
        f[idx[coffset + N, i3]] += z[ic] * u[idx[coffset + ic, i3]] / z[N]
    end
    if really_conserve
        for iv in 1:length(nv)
            uu = [u[idx[ic, iv]] for ic in 1:(N + 1)]
            c_num!(y, uu, data, ddata)
            for ic in 1:(N - 1)
                f[idx[ic + coffset, i3]] += y[ic] * nv[iv] / data.cscale
            end
        end
    end
    return nothing
end

# ╔═╡ 7bf3a130-3b47-428e-916f-4a0ec1237844
function ICMPBSystem(
        grid,
        data
    )

    if data.conserveions
        data.nv = ones(num_nodes(grid)) # trigger sparsity detector
        sys = VoronoiFVM.System(
            grid;
            data = data,
            flux = poisson_and_p_flux!,
            reaction = reaction!,
            bcondition = bcondition!,
            generic = ionconservation!,
            unknown_storage = :sparse
        )
    else
        sys = VoronoiFVM.System(
            grid;
            data = data,
            flux = poisson_and_p_flux!,
            reaction = reaction!,
            bcondition = bcondition!,
        )
    end

    for i in 1:data.N
        enable_species!(sys, i, [1])
    end

    enable_species!(sys, data.i0, [1])
    enable_species!(sys, data.iφ, [1])
    enable_species!(sys, data.ip, [1])

    if data.conserveions
        for ic in 1:data.N
            enable_boundary_species!(sys, data.coffset + ic, [3])
        end
        data.nv = nodevolumes(sys)
    end

    return sys
end;

# ╔═╡ 178b947f-3fef-44ed-9eca-fdb9916bc2b6
function qsweep(sys; qmax = 10, nsteps = 100, verbose = "", kwargs...)
    data = deepcopy(sys.physics.data)
    (; ip, iφ) = data
    apply_charge!(data, 0 * ph"e" / ufac"nm^2")
    state = VoronoiFVM.SystemState(sys; data)
    sol = solve!(state; inival = unknowns(sys, data), damp_initial = 0.1, verbose, kwargs...)

    volts = []
    Q = []
    for q in range(0, qmax, length = 50)
        apply_charge!(data, q * ph"e" / ufac"nm^2")
        sol = solve!(state; inival = sol, damp_initial = 0.1, verbose, kwargs...)
        push!(volts, (sol[iφ, end] - sol[iφ, 1]) / 2)
        # Division by  comes in because the voltage we get here is the difference
        # between the electrodes and not the difference between electrode and bulk
        # which corresponds to the other standard dlcap experiment
        push!(Q, q * ph"e" / ufac"nm^2")
    end
    dlcaps = -(Q[2:end] - Q[1:(end - 1)]) ./ (volts[2:end] - volts[1:(end - 1)])
    return volts[1:(end - 1)], dlcaps
end

# ╔═╡ fc84996b-02c0-4c16-8632-79f4e1900f78
function capscalc(sys, molarities; kwargs...)
    result = []
    for imol in 1:length(molarities)
        data = sys.physics.data
        set_molarity!(data, molarities[imol])

        t = @elapsed volts, caps = qsweep(sys; kwargs...)
        cdl0 = dlcap0(data)
        @info "elapsed=$(t)"
        push!(
            result,
            (
                voltages = volts,
                dlcaps = caps,
                cdl0 = cdl0,
                molarity = molarities[imol],
            ),
        )
    end
    return result
end

# ╔═╡ Cell order:
# ╠═ef660f6f-9de3-4896-a65e-13c60df5de1e
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╟─920b7d84-56c6-4958-aed9-fc67ba0c43f6
# ╟─87ac16f4-a4fc-4205-8fb9-e5459517e1b8
# ╟─7d77ad32-3df6-4243-8bad-b8df4126e6ea
# ╟─4cabef42-d9f9-43fe-988e-7b54462dc775
# ╠═0d825f88-cd67-4368-90b3-29f316b72e6e
# ╟─30c6a176-935b-423f-9447-86f78746322f
# ╠═a41c6c1f-ceb5-4590-a421-cae5078d167b
# ╟─f3049938-2637-401d-9411-4d7be07c19ca
# ╠═5d6340c4-2ddd-429b-a60b-3de5570a7398
# ╟─a21545da-3b53-47af-b0c4-f253b37dc84f
# ╠═1d22b09e-99c1-4026-9505-07bdffc98582
# ╟─5a210961-19fc-40be-a5f6-033a80f1414d
# ╠═fe704fb4-d07c-4591-b834-d6cf2f4f7075
# ╟─5eca37ba-f858-45fb-a66a-3795327dfd18
# ╟─a26cf11b-0ce1-4c1d-a64d-1917178ff676
# ╟─cdd1d359-08fa-45a1-a857-e19f2adefcab
# ╠═188f67d8-2ae8-474c-8e58-68b8b4fde02e
# ╟─f70eed13-a6c2-4d54-9f30-113367afaf7d
# ╠═d7531d5f-fc2d-42b2-9cf9-6a737b0f0f8d
# ╟─f6f004a6-d71b-4813-a363-9f51dc37e42a
# ╟─3810cc88-07f1-4741-853f-331e71c87923
# ╠═0e2d20a1-5f26-4263-9a91-3b40b2c2996a
# ╟─824c610b-6e5e-48a3-be37-19104f52d1d9
# ╟─2e11ce81-7d0d-498f-9ddd-7d4d836ab42f
# ╟─b1e062c6-f245-4edc-aa02-871e2c776998
# ╟─c4cc940c-74aa-45f8-a2fa-6016d7c3c145
# ╠═b07246b8-aec5-4161-8879-8cefb350aced
# ╟─b41838bb-3d5b-499c-9eb5-137c252ae366
# ╠═a468f43a-aa20-45dc-9c21-77f5adf2d700
# ╠═13fc2859-496e-4f6e-8b22-36d9d55768b8
# ╠═7c5e9686-9d66-4ff0-84a7-c7b22596ab57
# ╠═b1e333c0-cdaa-4242-b71d-b54ff71aef83
# ╟─55bd7b9a-a191-4a0b-9c6b-13733be5023e
# ╠═3ceda3b1-bf1c-4126-b94f-2ee03e8dde99
# ╠═97c5942c-8eb4-4b5c-8951-87ac0c9f396d
# ╟─0c54efd0-f279-4dc6-8b00-ba092dd13f44
# ╠═800dfed8-9f29-4138-96f8-e8bf1f2f00e6
# ╠═24910762-7d56-446b-a758-d8e830fe9a09
# ╟─9fe3ca93-c051-426e-8b9a-cc59f59319ad
# ╠═2ee34d76-7238-46c2-94d1-a40d8b017af6
# ╠═79cc671b-ef6e-42da-8641-61e43f221cb1
# ╟─7a607454-7b75-4313-920a-2dbdad258015
# ╟─9cb8324c-896f-40f8-baa8-b7d47a93e9f5
# ╟─003a5c0b-17c7-4407-ad23-21c0ac000fd4
# ╠═e1c13f1e-5b67-464b-967b-25e3a93e33d9
# ╠═64e47917-9c61-4d64-a6a1-c6e8c7b28c59
# ╠═743b9a7a-d6ac-4da0-8538-2045d965b547
# ╠═dbccaa88-65d9-47ab-be78-83df64a6db24
# ╠═0b646215-32db-4219-904b-f86f8861b46a
# ╠═7bf3a130-3b47-428e-916f-4a0ec1237844
# ╠═48670f54-d303-4c3a-a191-06e6592a2e0a
# ╠═b0a45e53-8b98-4e18-8b41-7f6d0bc1f76e
# ╠═178b947f-3fef-44ed-9eca-fdb9916bc2b6
# ╠═fc84996b-02c0-4c16-8632-79f4e1900f78
