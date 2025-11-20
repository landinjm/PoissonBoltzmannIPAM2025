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

const f_mod = true # model choice: 0: Boltzmann, 1: Bikerman
const nref = 4 # grid refinement level
const L = 10nm # computational domain size
const n_e = 10 # number of electrons/nm^2 at interfaces, defines q
const M_avg = 10 # average molarity
const E_0 = 10V / nm # decrement parameter
const a = 5.0 / E_0^2 # decrement parameter in χ(E)
const z = [-1, 1.0] # species charge numbers

const N = length(z)

const c_avg = fill(M_avg * mol / dm^3, N)

surfcharge(n) = n * ph"e" / ufac"nm^2"

const q = surfcharge(n_e) # surface charge

const c0_avg = c̄ - sum(c_avg) # solvent bulk molar concentration

const l_debye = sqrt((1 + χ_S) * ε_0 * RT / (F^2 * c_avg[1])) |> u"nm"

const hmin = 1.0e-1 * nm * 2.0^(-nref) # grid size at working electrode
const hmax = 1.0 * nm * 2.0^(-nref) # grid size at bulk

X0 = geomspace(0, L / 2, hmin, hmax)
X1 = geomspace(L / 2, L, hmax, hmin)
grid = simplexgrid(glue(X0, X1))
bfacemask!(grid, [L / 2], [L / 2], 3, tol = 1.0e-10 * nm)

const i3 = grid[BFaceNodes][1, 3] # Index  of grid midpoint

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

function relpermittivity(sol, data; grid)
    (; χ_S, a, iϕ) = data
    X = grid[Coordinates][1, :]
    return χ_S ./ (
        a *
            ((sol[iϕ, 2:end] - sol[iϕ, 1:(end - 1)]) ./ (X[2:end] - X[1:(end - 1)])) .^ 2 .+
            1
    ) .+ 1
end

function spacecharge!(y, ϕ, c_ref, c0_ref, data)
    (; N, F, c_ref, c0_ref, c̄, z) = data
    molfractions!(y, ϕ, c_ref, c0_ref, data)
    sumyz = zero(ϕ)
    for i in 1:N
        sumyz += z[i] * y[i]
    end
    return F * c̄ * sumyz
end

function reaction!(y, u, node, data)
    (; cache, c_ref, c0_ref, iϕ) = data
    tmp = get_tmp(data.cache, u)
    y[iϕ] = -spacecharge!(tmp, u[iϕ], c_ref, c0_ref, data)
    return nothing
end

function bcondition!(y, u, bnode, data)
    boundary_neumann!(y, u, bnode, species = data.iϕ, region = 2, value = -data.q)
    boundary_neumann!(y, u, bnode, species = data.iϕ, region = 1, value = data.q)
    return nothing
end

function MPBSystem(grid, data)
    return VoronoiFVM.System(
        grid;
        data,
        reaction = reaction!,
        flux = flux!,
        bcondition = bcondition!,
        species = [1],
    )
end

M1_ref = 1.0
n1_e = 10.0
M2_avg = 5.0
n2_e = 10.0
M3_avg = 5.0
n3_e = 10.0

data1 = PBData(c_ref = fill(M1_ref * mol / dm^3, 2), q = surfcharge(n1_e));

sys1 = MPBSystem(grid, data1);

sol1 = solve(sys1);

function extcref(cref0, data)
    (; z, N) = data
    return push!(copy(cref0), -z[1:(N - 1)] ⋅ cref0 / z[N])
end

lastsol = unknowns(sys1, inival = 0)

function clonedata(data0, c_ref)
    data = deepcopy(data0)
    data.c_ref .= c_ref
    data.c0_ref = c̄ - sum(data.c_ref)
    return data
end

cb = (0.01:0.01:1) * mol / dm^3

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

ca = [
    cavg(
            [c];
            q = surfcharge(n_e),
            sys = sys1,
            verbose = "",
            damp_initial = 0.1,
        )[1][1] for c in cb
]

let
    PythonPlot.clf()
    fig, ax = pyplot.subplots(1, 1)
    fig.set_size_inches(8, 2)
    ax.grid()
    ax.plot(cb, ca)
    PythonPlot.gcf()
end

c2_avg = [M2_avg * mol / dm^3]

res = nlsolve(
    c_ref -> cavg(
        c_ref;
        verbose = "",
        damp_initial = 0.1
    )[1] - c2_avg,
    c2_avg * 0.1,
    ftol = 1.0e-14
)

c2_ref = extcref(res.zero, VoronoiFVM.data(sys1))

sol2 = cavg(res.zero; sys = sys1, verbose = "")[2];

function ICMPBSystem(; data = PBData(), generic = VoronoiFVM.nofunc)
    sys = VoronoiFVM.System(
        grid;
        data,
        generic,
        flux = flux!,
        bcondition = bcondition!,
        unknown_storage = :sparse,
    )
    enable_species!(sys, data.iϕ, [1])
    for ic in 1:(data.N - 1)
        enable_boundary_species!(sys, data.iϕ + ic, [3])
    end
    return sys
end

sys3_0 = ICMPBSystem();

const nv = nodevolumes(sys3_0);

const idx = unknown_indices(unknowns(sys3_0));

idx[2, i3]

function xreaction!(f, u, sys, data)
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

c3_avg = fill(M3_avg * mol / dm^3, 2)

data3 = PBData(c_avg = c3_avg, c_ref = 0.5 * c3_avg, q = surfcharge(n3_e))

sys3 = ICMPBSystem(data = data3, generic = xreaction!)

state3 = VoronoiFVM.SystemState(sys3; data = data3)

state3.matrix

begin
    inival3 = unknowns(sys3, inival = 0.0)
    inival3[2:N, i3] .= c3_avg[1:(N - 1)] / 2
end;

sol3 = solve!(state3; inival = inival3, verbose = "n", damp_initial = 0.5)

myround(x) = round(x, sigdigits = 4)

function plotsol(
        sol,
        sys;
        data = data(sys),
        grid = grid,
        c_ref = data.c_ref,
        size = (600, 400),
        id = 0,
    )
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
    ax1.set_title(
        "ϕ∈$(round.(Float64.(extrema(sol[1, :])), sigdigits = 3)), ε_r ∈$(round.(Float64.(extrema(ε_r)), sigdigits = 3))",
    )
    ax1.plot(X / nm, sol[1, :], color = "green", linewidth = 2, label = "ϕ")
    ax1r.plot(
        X[1:(end - 1)] / nm,
        ε_r,
        color = "pink",
        linewidth = 3,
        label = L"ε_r",
    )
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

    ax2.plot(
        X / nm,
        c[1, :] / (mol / dm^3),
        color = "blue",
        linewidth = 2,
        label = L"c^-",
    )
    ax2.plot(
        X / nm,
        c[2, :] / (mol / dm^3),
        color = "red",
        linewidth = 2,
        label = L"c^+",
    )
    ax2.plot(
        X / nm,
        c0[1, :] / (mol / dm^3),
        color = "green",
        linewidth = 2,
        label = L"c_{solvent}",
    )
    ax2.legend(loc = (0.4, 0.1))

    tight_layout()
    savefig("ICMPB_$id.jpg", dpi = 300)
    return PythonPlot.gcf()
end

plotsol(sol1, sys1; size = (600, 400), id = 1)

plotsol(
    sol2,
    sys1,
    data = data(sys1);
    c_ref = c2_ref,
    size = (600, 400),
    id = 2,
)

plotsol(
    sol3,
    sys3;
    data = data3,
    c_ref = extcref(sol3[2:N, i3], data3),
    size = (600, 400),
    id = 3,
)
