using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

using VoronoiFVM
using ExtendableGrids
using LinearAlgebra
using NLsolve
using Unitful
using LessUnitful
using LessUnitful.MoreUnitful
using Test
using PythonPlot
using Colors
using JuliaMPBSolver

function L_Debye(data)
    return sqrt(
        (1 + data.χ) * data.ε_0 * ph"k_B" * data.T / (data.e^2 * data.n_E[1]),
    )
end

const mol = ph"N_A"

begin
    const iφ = 1
    const ip = 2
    const iA = 1
    const iC = 2
end

function y_α(φ, p, α, data)
    η_φ = data.z[α] * data.e * (φ - data.E_ref)
    η_p = data.v[α] * (p * data.pscale - data.p_ref)
    return data.y_E[α] * exp(-(η_φ + η_p) / (data.kT))
end

function y0(p, data)
    return data.y0_E * exp(-data.v0 * (p * data.pscale - data.p_ref) / (data.kT))
end

function poisson_flux!(f, u, edge, data)
    return f[iφ] = (1.0 + data.χ) * data.ε_0 * (u[iφ, 1] - u[iφ, 2])
end

function spacecharge(φ, p, data)
    y = y0(p, data)
    sumyz = zero(eltype(p))
    sumyv = data.v0 * y
    for α in 1:(data.N)
        y = y_α(φ, p, α, data)
        sumyz += data.z[α] * y
        sumyv += data.v[α] * y
    end
    return data.e * sumyz / sumyv
end

function ysum(φ, p, data)
    sumy = y0(p, data)
    for α in 1:(data.N)
        sumy += y_α(φ, p, α, data)
    end
    return sumy
end

function derived(κ, v0, vu, n_E, T)
    c0 = 1 / v0
    barc = 0.0
    v = vu + κ * v0
    N = length(κ)
    for α in 1:N
        barc += n_E[α]
        c0 -= n_E[α] * (1 + κ[α])
    end
    barc += c0
    y_E = n_E / barc
    y0_E = c0 / barc
    return (; v, y_E, y0_E)
end

begin
    """
        EquilibriumData

    Data structure containing data for equilibrium calculations
    """
    Base.@kwdef mutable struct EquilibriumData
        N::Int64 = 2                     # number of ionic species
        T::Float64 = 298.15 * ufac"K"        # temperature
        kT::Float64 = ph"k_B" * T             # temperature
        p_ref::Float64 = 1.0e5 * ufac"Pa"        # reference pressure
        pscale::Float64 = 1.0 * ufac"GPa"         # pressure scaling nparameter
        E_ref::Float64 = 0.0 * ufac"V"           # reference voltage
        n0_ref::Float64 = 55.508 * ph"N_A" / ufac"dm^3"  # solvent molarity
        v0::Float64 = 1 / n0_ref              # solvent molecule volume
        vu::Vector{Float64} = [1 / n0_ref, 1 / n0_ref]           # unsolvated ion volume
        χ::Float64 = 15                    # dielectric susceptibility
        z::Vector{Int} = [-1, 1]                # ion charge numbers
        κ::Vector{Int} = [10, 10]               # ion solvation numbers
        molarity::Float64 = 0.1 * ph"N_A" / ufac"dm^3"
        n_E::Vector{Float64} = [molarity, molarity]  # bulk ion number densities
        μ_e::Vector{Float64} = [0.0]             # grain facet electron chemical potential

        e::Float64 = ph"e"
        ε_0::Float64 = ph"ε_0"

        v::Vector{Float64} = derived(κ, v0, vu, n_E, T).v   # effective ion volumes
        y_E::Vector{Float64} = derived(κ, v0, vu, n_E, T).y_E # bulk ion mole fractions
        y0_E::Float64 = derived(κ, v0, vu, n_E, T).y0_E       # bulk solvent mole fraction
    end
end

function set_molarity!(data::EquilibriumData, M_E)
    n_E = M_E * ph"N_A" / ufac"dm^3"
    data.molarity = n_E
    return data.n_E = fill(n_E, data.N)
end

function dlcap0(data::EquilibriumData)
    return sqrt(
        2 * (1 + data.χ) * ph"ε_0" * ph"e"^2 * data.n_E[1] / (ph"k_B" * data.T),
    )
end

begin
    function update_derived!(data::EquilibriumData)
        (; κ, v0, vu, n_E, T) = data
        return data.v, data.y_E, data.y0_E = derived(κ, v0, vu, n_E, T)
    end
end

let
    data = EquilibriumData()
    set_molarity!(data, 0.01)
    update_derived!(data)
    sumyz = 0.0
    sumyv = data.y0_E * data.v0
    sumy = data.y0_E
    for α in 1:data.N
        v = (1.0 + data.κ[α]) * data.v0
        sumyz += data.y_E[α] * data.z[α]
        sumyv += data.y_E[α] * v
        sumy += data.y_E[α]
    end
    @assert sumy ≈ 1.0
end

φ_Σ(ifacet, data, E) = data.μ_e[ifacet] / data.e - (E - data.E_ref)

function apply_voltage!(sys, E)
    data = sys.physics.data
    nbc = num_bfaceregions(sys.grid)
    nfacets = length(data.μ_e)
    @assert nbc > nfacets
    for ifacet in 1:nfacets
        boundary_dirichlet!(sys, iφ, ifacet, φ_Σ(ifacet, data, E))
    end
    return sys
end

calc_φ(sol, sys) = sol[iφ, :]

calc_p(sol, sys) = sol[ip, :] * sys.physics.data.pscale

function c_num!(c, φ, p, data)
    y = y0(p, data)
    sumyv = data.v0 * y
    for α in 1:(data.N)
        c[α] = y_α(φ, p, α, data)
        sumyv += c[α] * data.v[α]
    end
    return c ./= sumyv
end

function c0_num!(c, φ, p, data)
    y = y0(p, data)
    sumyv = data.v0 * y
    for α in 1:(data.N)
        c[α] = y_α(φ, p, α, data)
        sumyv += c[α] * data.v[α]
    end
    return y / sumyv
end

function calc_cnum(sol, sys)
    data = sys.physics.data
    grid = sys.grid
    nnodes = num_nodes(grid)
    conc = zeros(data.N, nnodes)
    for i in 1:nnodes
        @views c_num!(conc[:, i], sol[iφ, i], sol[ip, i], data)
    end
    return conc
end

function calc_c0num(sol, sys)
    data = sys.physics.data
    grid = sys.grid
    nnodes = num_nodes(grid)
    c0 = zeros(nnodes)
    conc = zeros(data.N)
    for i in 1:nnodes
        @views c0[i] = c0_num!(conc, sol[iφ, i], sol[ip, i], data)
    end
    return c0
end

calc_cmol(sol, sys) = calc_cnum(sol, sys) / (ph"N_A" * ufac"mol/dm^3")

calc_c0mol(sol, sys) = calc_c0num(sol, sys) / (ph"N_A" * ufac"mol/dm^3")

function spacecharge!(f, u, node, data)
    φ = u[iφ]
    p = u[ip]
    return f[iφ] = -spacecharge(u[iφ], u[ip], data)
end

function poisson_and_p_flux!(f, u, edge, data)
    f[iφ] = (1.0 + data.χ) * data.ε_0 * (u[iφ, 1] - u[iφ, 2])
    q1 = spacecharge(u[iφ, 1], u[ip, 1], data)
    q2 = spacecharge(u[iφ, 2], u[ip, 2], data)
    return f[ip] =
        (u[ip, 1] - u[ip, 2]) +
        (u[iφ, 1] - u[iφ, 2]) * (q1 + q2) / (2 * data.pscale)
end

function ysum(sys, sol)
    data = sys.physics.data
    n = size(sol, 2)
    sumy = zeros(n)
    for i in 1:n
        sumy[i] = ysum(sol[iφ, i], sol[ip, i], data)
    end
    return sumy
end

function spacecharge_and_ysum!(f, u, node, data)
    φ = u[iφ]
    p = u[ip]
    f[iφ] = -spacecharge(φ, p, data)
    return f[ip] = log(ysum(φ, p, data)) # this behaves much better with Newton's method
end

function create_equilibrium_system(
        grid,
        data::EquilibriumData = EquilibriumData();
        Γ_bulk = 0,
    )
    update_derived!(data)
    sys = VoronoiFVM.System(
        grid;
        data = data,
        flux = poisson_flux!,
        reaction = spacecharge_and_ysum!,
        species = [iφ, ip],
    )
    if Γ_bulk > 0
        boundary_dirichlet!(sys, iφ, Γ_bulk, 0.0)
    end
    return apply_voltage!(sys, 0)
end

function calc_QBL(sol, sys)
    return VoronoiFVM.integrate(sys, spacecharge_and_ysum!, sol)[iφ, 1]
end

function dlcapsweep_equi(
        sys;
        vmax = 2 * ufac"V",
        nsteps = 21,
        δV = 1.0e-3 * ufac"V",
        molarity = nothing,
        verbose = false,
    )
    if !isnothing(molarity)
        error(
            "The molarity kwarg of dlcapsweep_equie has been removed. Pass the molarity information with set_molarity!.",
        )
    end

    data = sys.physics.data
    update_derived!(data)
    apply_voltage!(sys, 0)

    c = VoronoiFVM.NewtonControl()
    #	c.damp_growth=1.1
    c.verbose = verbose
    c.tol_round = 1.0e-10
    c.max_round = 3
    c.damp_initial = 0.01
    c.damp_growth = 2

    inival = solve(sys; inival = 0, control = c)
    vstep = vmax / (nsteps - 1)

    c.damp_initial = 1

    function rundlcap(dir)
        volts = zeros(0)
        caps = zeros(0)
        volt = 0.0
        sol = inival
        for iv in 1:nsteps
            apply_voltage!(sys, volt)
            c.damp_initial = 1
            sol = solve(sys; inival = sol, control = c)

            Q = calc_QBL(sol, sys)
            apply_voltage!(sys, volt + dir * δV)
            c.damp_initial = 1
            sol = solve(sys; inival = sol, control = c)
            Qδ = calc_QBL(sol, sys)
            push!(caps, (Q - Qδ) / (dir * δV))
            push!(volts, volt)
            volt += dir * vstep
        end
        return volts, caps
    end
    Vf, Cf = rundlcap(1)
    Vr, Cr = rundlcap(-1)
    return vcat(reverse(Vr), Vf), vcat(reverse(Cr), Cf)
end

function create_equilibrium_pp_system(
        grid,
        data::EquilibriumData = EquilibriumData();
        Γ_bulk = 0,
    )
    update_derived!(data)

    sys = VoronoiFVM.System(
        grid;
        data = data,
        flux = poisson_and_p_flux!,
        reaction = spacecharge!,
        species = [iφ, ip],
    )
    if Γ_bulk > 0
        logysum!(y, p) = y[1] = log(ysum(0.0, p[1], data))
        res = nlsolve(
            logysum!,
            [0.0];
            autodiff = :forward,
            method = :newton,
            xtol = 1.0e-10,
            ftol = 1.0e-20,
        )
        boundary_dirichlet!(sys, iφ, Γ_bulk, 0.0)
        @info res
        boundary_dirichlet!(sys, ip, Γ_bulk, res.zero[1])
    end
    return apply_voltage!(sys, 0)
end

begin
    SI(x) = Float64(Unitful.ustrip(Unitful.upreferred(1 * x)))
    const V = SI(Unitful.V)
    const eV = SI(Unitful.eV)
    const nm = SI(Unitful.nm)
    const cm = SI(Unitful.cm)
    const μF = SI(Unitful.μF)
end

L_Debye(EquilibriumData(molarity = 0.01mol / ufac"dm^3")) / nm

dlcap_exact = 0.22846691848825248

begin
    equidata = EquilibriumData()
    set_molarity!(equidata, 0.01)
    equidata.χ = 78.49 - 1
end

@test dlcap0(equidata) |> unitfactor ≈ dlcap_exact

begin
    Vmax = 2 * V

    L = 20nm
end

grid_parameters = JuliaMPBSolver.Grid.GeometricGrid(
    domain_size = L,
    refinement = 0,
    hmin = 0.05 * nm,
    hmax = 0.5 * nm,
    use_offset = false,
)
grid = JuliaMPBSolver.Grid.create_half_cell(grid_parameters)
X = JuliaMPBSolver.Grid.get_coordinates(grid)

sys_sy = create_equilibrium_system(grid, equidata)

sys_pp = create_equilibrium_pp_system(grid, equidata, Γ_bulk = 2)

inival = unknowns(sys_sy, inival = 0)

molarities = [0.001, 0.01, 0.1, 1]

function capscalc(sys, molarities)
    result = []
    for imol in 1:length(molarities)
        data = sys.physics.data
        set_molarity!(data, molarities[imol])
        t = @elapsed volts, caps = dlcapsweep_equi(sys, vmax = 1V, nsteps = 101)
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

result_sy = capscalc(sys_sy, molarities)

result_pp = capscalc(sys_pp, molarities)

function resultcompare(r1, r2; tol = 1.0e-3)
    for i in 1:length(r1)
        for f in fieldnames(typeof(r1[i]))
            if !isapprox(r1[i][f], r2[i][f]; rtol = tol)
                return false
            end
        end
    end
    return true
end

@test resultcompare(result_pp, result_sy; tol = 1.0e-2)

function capsplot(ax, result, title)
    hmol = 1 / length(result)
    ax.set_title(title)
    ax.grid()
    ax.set_xlabel("U/V")
    ax.set_ylabel(L"C_{dl}/(μF/cm^3)")
    for imol in 1:length(result)
        c = (1 - imol * hmol, 0, imol * hmol)

        ax.plot(
            result[imol].voltages,
            result[imol].dlcaps / (μF / cm^2),
            color = c,
            label = "$(result[imol].molarity)M",
        )
        ax.scatter([0], [result[imol].cdl0] / (μF / cm^2), color = c)
    end
    return
end

let
    res = (600, 200)
    fig = figure(1, figsize = (6, 3))
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    capsplot(ax1, result_sy, "Algebraic pressure")
    capsplot(ax2, result_pp, "Pressure Poisson")

    tight_layout()
    savefig("equilibrium.jpg", dpi = 300)
end
