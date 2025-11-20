using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))
using Revise
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
using DelimitedFiles

using DrWatson

using PoissonBoltzmannIPAM2025
using JuliaMPBSolver

const nel = 1.0 * JuliaMPBSolver.Units.el_surface_density # number of electrons/nm^2 at interfaces
const M_bulk = 2.0 # (bulk) molarity at center of domain
const E0 = 10JuliaMPBSolver.Units.V / JuliaMPBSolver.Units.nm # decrement parameter
const a = 5.0 / E0^2 # decrement parameter in χ(E)
const c̄ = 55.508JuliaMPBSolver.Units.M # summary molar concentration
const z = [-1, 1]
const c_bulk = [M_bulk / abs(z[1]), M_bulk / abs(z[2])] * JuliaMPBSolver.Units.M # bulk  concentrations

const F = ph"N_A" * ph"e"
const K = ufac"K"
const nm = ufac"nm"
const m = ufac"m"
const dm = ufac"dm"
const V = ufac"V"
const mol = ufac"mol"
const T = (273.15 + 25) * ufac"K"
const RT = ph"R" * T
const ε_0 = ph"ε_0"

const f_mod = true # model choice: 0: Boltzmann, 1: Bikerman
const L = 2.278592867nm # computational domain size

const N = length(z)

# Parameters
user_parameters = JuliaMPBSolver.Parameters.UserParameters(
    273.15 + 25 * JuliaMPBSolver.Units.K,
    78.49 - 1,
    0.0,
    0.0,
    z,
    c̄,
    c_bulk,
    nel,
    a,
    true,
    true,
    false,
    0.0,
)

computed_parameters =
    JuliaMPBSolver.Parameters.ComputedParameters(user_parameters)

grid_parameters = JuliaMPBSolver.Grid.UniformGrid(
    domain_size = L,
    refinement = 4,
    n_points = 20,
    use_offset = false,
)
grid = JuliaMPBSolver.Grid.create_full_cell(grid_parameters)

const i3 = grid[BFaceNodes][1, 3] # Index  of grid midpoint

Base.@kwdef mutable struct PBData
    c_ref::Vector{Float64} = c_bulk
    c_avg::Vector{Float64} = c_bulk
    c̄::Float64 = c̄
    c0_ref::Float64 = c̄ - sum(c_ref)
    q::Float64 = 0.0
    z::Vector{Float64} = z
    N::Int = length(z)
    a::Float64 = a
    f_mod::Bool = f_mod
    F::Float64 = F
    RT::Float64 = RT
    ε_0::Float64 = ε_0
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

function flux!(y, u, edge, data)
    (; a, ε_0, iϕ) = data
    eins = one(eltype(u))
    h = edgelength(edge)
    E = (u[iϕ, 1] - u[iϕ, 2]) / h
    χ = user_parameters.dielectric_susceptibility / sqrt(eins + a * E^2)
    ε = (eins + χ) * ε_0
    y[iϕ] = ε * ((u[iϕ, 1] - u[iϕ, 2]))
    return nothing
end

function extcref(cref0, data)
    (; z, N) = data
    return push!(copy(cref0), -z[1:(N - 1)] ⋅ cref0 / z[N])
end

function ICMPBSystem(; data = PBData(), generic = VoronoiFVM.nofunc)
    sys = VoronoiFVM.System(
        grid;
        data,
        generic,
        flux = flux!,
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

Z = grid[Coordinates] / nm

indata = (M = M_bulk, q = 1, L = L / nm, n = length(Z))

savename(indata)

c3_avg = fill(M_bulk * JuliaMPBSolver.Units.M, 2)

data3 = PBData(
    c_avg = c3_avg,
    c_ref = 0.5 * c3_avg,
    q = user_parameters.boundary_electron_density,
)

sys3 = ICMPBSystem(data = data3, generic = xreaction!)

JuliaMPBSolver.Equations.add_boundary_charge!(
    sys3,
    1,
    2,
    -user_parameters.boundary_electron_density,
)
JuliaMPBSolver.Equations.add_boundary_charge!(
    sys3,
    1,
    1,
    user_parameters.boundary_electron_density,
)

state3 = VoronoiFVM.SystemState(sys3; data = data3)

begin
    inival3 = unknowns(sys3, inival = 0.0)
    inival3[2:N, i3] .= c3_avg[1:(N - 1)] / 2
end;

sol3 = solve!(state3; inival = inival3, verbose = "n", damp_initial = 0.5)

conc =
    concentrations(sol3, data3, c_ref = extcref(sol3[2:N, i3], data3)) /
    (JuliaMPBSolver.Units.M)

c0 = c̄ / (JuliaMPBSolver.Units.M) .- sum(conc, dims = 1)

reference_data_filename = resultsdir(savename("icmpb", indata, "hdf5"))
print(reference_data_filename)

if !isfile(reference_data_filename)
    throw(ErrorException("Reference data not found"))
end

(x_ref, c_solvent_ref) =
    JuliaMPBSolver.DataOut.read_hdf5_data(reference_data_filename, "c_solvent")
(x_ref, c_anion_ref) =
    JuliaMPBSolver.DataOut.read_hdf5_data(reference_data_filename, "c_anion")
(x_ref, c_cation_ref) =
    JuliaMPBSolver.DataOut.read_hdf5_data(reference_data_filename, "c_cation")

same_x = x_ref == JuliaMPBSolver.Grid.get_coordinates(grid)
if !same_x
    throw(ErrorException("Grid coordinates mismatch"))
end
same_solvent_concentration = c_solvent_ref == c0
if !same_solvent_concentration
    throw(ErrorException("Solvent concentration mismatch"))
end
same_anion_concentration = c_anion_ref == conc[1, :]
if !same_anion_concentration
    throw(ErrorException("Anion concentration mismatch"))
end
same_cation_concentration = c_cation_ref == conc[2, :]
if !same_cation_concentration
    throw(ErrorException("Cation concentration mismatch"))
end
