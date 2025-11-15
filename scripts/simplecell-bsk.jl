using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

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
using JuliaMPBSolver

begin
  const F = ph"N_A" * ph"e"
  const K = ufac"K"
  const nm = ufac"nm"
  const m = ufac"m"
  const dm = ufac"dm"
  const V = ufac"V"
  const mol = ufac"mol"
  const T = (273.15 + 25) * ufac"K"
  const RT = ph"R" * T
end

begin
  const floattype = Float64
  const L = 10.0nm # computational domain size
  const nel = 20 # number of electrons/nm^2
  const lc = 1 * nm # BSK parameter
  const M_bulk = 1 # (bulk) molarity at center of domain
  const E0 = floattype(10V / nm) # decrement parameter
  const a = 5.0 / E0^2 # decrement parameter in χ(E)
  const q = nel * ph"e" / ufac"nm^2" # surface charge
  const c̄ = 55.508mol / dm^3 # summary molar concentration
  const ε_0 = floattype(ph"ε_0")
end

const z = [-1, 1]
const c_bulk = [M_bulk, M_bulk] * mol / dm^3 # bulk  concentrations
const c0_bulk = c̄ - sum(c_bulk) # solvent bulk molar concentration

# Parameters
user_parameters = JuliaMPBSolver.Parameters.UserParameters(
  273.15 + 25 * JuliaMPBSolver.Units.K,
  78.49 - 1,
  0.0,
  0.0,
  z,
  c̄,
  c_bulk,
  c0_bulk,
  true,
)

begin
  const l_debye = sqrt(
    (1 + user_parameters.dielectric_susceptibility) * ε_0 * RT /
    (F^2 * c_bulk[1]),
  ) # Debye length
  const dlcap0 = sqrt(
    2 *
    (1 + user_parameters.dielectric_susceptibility) *
    ε_0 *
    F^2 *
    c_bulk[1] / RT,
  ) # Double layer capacitance at point of zero charge (0V)
end

grid_parameters = JuliaMPBSolver.Grid.GeometricGrid(
  domain_size = L,
  refinement = 4,
  hmin = 1.0e-1 * nm,
  hmax = 1.0 * nm,
  use_offset = false,
)
grid = JuliaMPBSolver.Grid.create_full_cell(grid_parameters)
X = JuliaMPBSolver.Grid.get_coordinates(grid)

const Y = DiffCache(ones(floattype, length(z))) # place for temporary data in callbacks

function flux!(y, u, edge, data)
  eins = one(eltype(u))
  h = floattype(edgelength(edge))
  E = (u[1, 1] - u[1, 2]) / h
  χ = user_parameters.dielectric_susceptibility / sqrt(eins + a * E^2)
  ε = (eins + χ) * ε_0
  y[1] = ε * ((u[1, 1] - u[1, 2]) - lc^2 * (u[2, 1] - u[2, 2]))
  y[2] = u[1, 1] - u[1, 2]
  return nothing
end

function reaction!(y, u, node, data)
  tmp = get_tmp(Y, u)
  y[1] =
    -JuliaMPBSolver.Postprocess.compute_spacecharge(tmp, u[1], user_parameters)
  y[2] = -u[2]
  return nothing
end

function bcondition!(y, u, bnode, data)
  boundary_neumann!(y, u, bnode, species = 1, region = 2, value = -q)
  boundary_neumann!(y, u, bnode, species = 1, region = 1, value = q)
  boundary_dirichlet!(y, u, bnode, species = 2, region = 3, value = 0)
  return nothing
end

pbsystem = VoronoiFVM.System(
  grid;
  reaction = reaction!,
  flux = flux!,
  bcondition = bcondition!,
  species = [1, 2],
  valuetype = floattype,
)

sol = solve(
  pbsystem,
  inival = 0.1,
  verbose = "n",
  damp_initial = 0.1,
  maxiters = 1000,
)

function bee!(y, ϕ)
  N = length(z)
  for i in 1:N
    y[i] = RT * log(c_bulk[i] / c0_bulk) / F - z[i] * ϕ
  end
  return nothing
end

function bee(sol)
  n = size(sol, 2)
  N = length(z)
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

nv = nodevolumes(pbsystem)

c =
  JuliaMPBSolver.Postprocess.compute_concentrations(sol[1, :], user_parameters)

ε_r =
  user_parameters.dielectric_susceptibility ./ (
    a *
    ((sol[1, 2:end] - sol[1, 1:(end-1)]) ./ (X[2:end] - X[1:(end-1)])) .^ 2 .+
    1
  ) .+ 1

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
  )
  cm = c[1, :] ⋅ nv / (mol / dm^3) / L
  cp = c[2, :] ⋅ nv / (mol / dm^3) / L
  c0 = -(sum(c, dims = 1) .- c̄)
  e = bee(sol)
  ax1.set_title(
    "ϕ∈$(round.(Float64.(extrema(sol[1, :])), sigdigits = 3)), ε_r ∈$(round.(Float64.(extrema(ε_r)), sigdigits = 3))",
  )
  ax1.plot(X / nm, sol[1, :], color = "green", linewidth = 2, label = "ϕ")
  ax1.plot(
    X / nm,
    e[1, :],
    color = "blue",
    linewidth = 2,
    label = L"ψ^+",
    linestyle = "dotted",
  )
  ax1.plot(
    X / nm,
    e[2, :],
    color = "red",
    linewidth = 2,
    linestyle = "dotted",
    label = L"ψ^+",
  )
  ax1r.plot(
    X[1:(end-1)] / nm,
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
  savefig("simplecell-bsk.jpg", dpi = 300)
  return PythonPlot.gcf()
end

plotsol(sol)
