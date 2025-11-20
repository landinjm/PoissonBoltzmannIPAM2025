using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

using LinearAlgebra
using Interpolations
using VoronoiFVM
using PythonPlot
using JuliaMPBSolver

const nel = 20.0 * JuliaMPBSolver.Units.el_surface_density # number of electrons/nm^2 at interfaces
const M_bulk = 1 # (bulk) molarity at center of domain
const E0 = 10JuliaMPBSolver.Units.V / JuliaMPBSolver.Units.nm # decrement parameter
const a = 5.0 / E0^2 # decrement parameter in χ(E)
const c̄ = 55.508JuliaMPBSolver.Units.M # summary molar concentration
const z = [-1, 1]
const c_bulk = [M_bulk / abs(z[1]), M_bulk / abs(z[2])] * JuliaMPBSolver.Units.M # bulk  concentrations

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
)

computed_parameters =
    JuliaMPBSolver.Parameters.ComputedParameters(user_parameters)

# Grid generation
grid_parameters = JuliaMPBSolver.Grid.GeometricGrid(
    domain_size = 10.0 * JuliaMPBSolver.Units.nm,
    refinement = 4,
    hmin = 1.0e-1 * JuliaMPBSolver.Units.nm,
    hmax = 1.0 * JuliaMPBSolver.Units.nm,
    use_offset = false,
)

solution, X, nv, ε_r =
    JuliaMPBSolver.Equations.create_and_run_full_cell_problem(
    grid_parameters,
    user_parameters,
    computed_parameters,
)

function bee!(y, ϕ)
    N = length(user_parameters.charge_numbers)
    for i in 1:N
        y[i] =
            JuliaMPBSolver.Units.thermal_energy(user_parameters.temperature) *
            log(c_bulk[i] / computed_parameters.bulk_solvent_concentration) /
            JuliaMPBSolver.Units.F - user_parameters.charge_numbers[i] * ϕ
    end
    return nothing
end

function bee(sol)
    n = size(sol, 2)
    N = length(user_parameters.charge_numbers)
    e = zeros(N, n)
    y = zeros(N)
    for i in 1:n
        bee!(y, sol[1, i])
        for j in 1:N
            e[j, i] = y[j]
        end
    end
    return e
end

function plotsol(sol; size = (600, 400))
    PythonPlot.clf()
    fig, ax = pyplot.subplots(2, 1)
    ax1 = ax[0]
    ax2 = ax[1]
    ax1.grid()
    ax1r = ax1.twinx()

    c = JuliaMPBSolver.Postprocess.compute_concentrations(
        sol[1, :],
        user_parameters,
        computed_parameters,
    )
    cm = c[1, :] ⋅ nv / (JuliaMPBSolver.Units.M) / grid_parameters.domain_size
    cp = c[2, :] ⋅ nv / (JuliaMPBSolver.Units.M) / grid_parameters.domain_size
    c0 = -(sum(c, dims = 1) .- c̄)
    e = bee(sol)
    ax1.set_title(
        "ϕ∈$(round.(Float64.(extrema(sol[1, :])), sigdigits = 3)), ε_r ∈$(round.(Float64.(extrema(ε_r)), sigdigits = 3))",
    )
    ax1.plot(
        X / JuliaMPBSolver.Units.nm,
        sol[1, :],
        color = "green",
        linewidth = 2,
        label = "ϕ",
    )
    ax1.plot(
        X / JuliaMPBSolver.Units.nm,
        e[1, :],
        color = "blue",
        linewidth = 2,
        label = L"ψ^+",
        linestyle = "dotted",
    )
    ax1.plot(
        X / JuliaMPBSolver.Units.nm,
        e[2, :],
        color = "red",
        linewidth = 2,
        linestyle = "dotted",
        label = L"ψ^+",
    )
    ax1r.plot(
        X[1:(end - 1)] / JuliaMPBSolver.Units.nm,
        ε_r,
        color = "pink",
        linewidth = 3,
        label = L"ε_r",
    )
    ax1.set_ylim(-10, 10)
    ax1.set_xlabel("z/nm")
    ax1.set_ylabel("ϕ/V")
    ax1.legend(loc = (0.1, 0.1))
    ax1r.legend(loc = (0.8, 0.1))
    ax1r.set_ylim(0, 80)

    ax2.grid()
    ax2.set_title("M_avg=$(round.((cm, cp), sigdigits = 3))")
    ax2.set_xlabel("z/nm")
    ax2.set_ylabel("c/(mol/L)")
    ax2.set_ylim(0, 60)

    ax2.plot(
        X / JuliaMPBSolver.Units.nm,
        c[1, :] / (JuliaMPBSolver.Units.M),
        color = "blue",
        linewidth = 2,
        label = L"c^-",
    )
    ax2.plot(
        X / JuliaMPBSolver.Units.nm,
        c[2, :] / (JuliaMPBSolver.Units.M),
        color = "red",
        linewidth = 2,
        label = L"c^+",
    )
    ax2.plot(
        X / JuliaMPBSolver.Units.nm,
        c0[1, :] / (JuliaMPBSolver.Units.M),
        color = "green",
        linewidth = 2,
        label = L"c_{solvent}",
    )
    ax2.legend(loc = (0.4, 0.1))

    tight_layout()
    savefig("simplecell.jpg", dpi = 300)
    return PythonPlot.gcf()
end

plotsol(solution)
