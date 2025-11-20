### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 8871344f-1ded-4d21-86b6-8a61a7038aa1
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
    using PythonPlot
    using LaTeXStrings
    using DoubleFloats
    using ForwardDiff
end

# ╔═╡ d1503e2a-76c8-4b40-9650-9121c9b771e7
begin
    using Markdown
    using InteractiveUtils
end

# ╔═╡ d5d26fae-c502-11f0-8122-b54486d2dc8d
begin
    using NLsolve


    # --- Physical constants (SI) ---
    const e0 = 1.602176634e-19    # elementary charge (C)
    const kB = 1.380649e-23       # Boltzmann constant (J/K)
    const NA = 6.02214076e23      # Avogadro
    const eps0 = 8.8541878128e-12 # vacuum permittivity (F/m)
    const Rgas = kB * NA          # gas constant (J/(mol K))

    # --- Model parameters (use values from Dreyer et al. 2014) ---
    T = 293.75                     # temperature (K)
    kBT = kB * T
    RT = Rgas * T

    # electrolyte: symmetric 1:1 (A = anion, C = cation)
    zA = -1; zC = 1

    # solvation parameters (from Dreyer et al. paper example)
    κA = 15.0; κC = 15.0           # solvation numbers
    vR0 = 1.797e-5                 # L/mol -> m^3/mol (1.797e-2 L/mol = 1.797e-5 m^3/mol)
    vRA = vR0                      # in the example vRα = vR0
    vRC = vR0

    # bulk molarities (mol/m^3) — paper uses 0.5 mol/L = 500 mol/m^3
    nB_A = 0.5e3   # mol/m^3
    nB_C = 0.5e3
    nB_0 = 55.508e3  # solvent (mol/m^3) ~ 55.508 mol/L

    # compute bulk eyB (volume fractions used in from Dreyer et al. paper notation)
    en_bulk = nB_0 + nB_A + nB_C
    eyB_A = nB_A / en_bulk
    eyB_C = nB_C / en_bulk
    eyB_0 = nB_0 / en_bulk

    # permittivity (field-independent here; we can make it field-dependent later)
    eps_r_bulk = 80.0
    const εr = eps_r_bulk
    const ε = eps0 * εr

    # BSK length scale (choose 1 nm as in BSK's notebook): convert nm -> m
    lc = 1.0e-9

    # domain and grid (same idea as in BSK's notebook)
    L = 10.0e-9    # domain length (m) — half cell like you had before; adjust as needed
    npts = 201     # number of grid points
    X = range(0.0, stop = L, length = npts) |> collect
    dx = X[2] - X[1]

    # boundary values
    phi_left = 0.4      # applied potential at electrode (V) —
    phi_right = 0.0      # bulk potential (V)
    p_bulk = 1.01325e5   # bulk pressure (Pa)

    # helper indices for unknown vector ordering
    # unknowns per node: [phi, psi, p, yA, yC]
    Nvar = 5
    N = npts
    tot = Nvar * N
    # index helper
    idx(i, v) = (i - 1) * Nvar + v   # v in 1..5

end


# ╔═╡ 371f9c44-444d-4d10-a0d3-79befef1698d
begin
    # Precompute some constants used in exponentials
    const alphaA = (vRA + κA * vR0) / (kBT)  # factor multiplying (p-pB) in exponent for anion
    const alphaC = (vRC + κC * vR0) / (kBT)

    # Bulk potentials used in eq. (27)
    phiB = phi_right
    pB = p_bulk

    # q(φ,p) according to eq. (28) specialized to binary electrolyte (A,C)
    function compute_q(phi, p, yA, yC)
        # ey0 from incompressibility (ey0 = 1 - yA - yC)
        y0 = 1.0 - yA - yC
        # numerator: e0 * sum zα eyα
        number = e0 * (zA * yA + zC * yC)
        # denominator: vR0*ey0 + sum (κ_α vR0 + vRα) eyα
        denom = vR0 * y0 + (κA * vR0 + vRA) * yA + (κC * vR0 + vRC) * yC
        # convert to charge density (C per m^3)
        return number / denom
    end

    # algebraic expressions for y from eq. (27)
    function y_from_phi_p(phi, p)
        # eyα(x) = eyBα * exp[- zα e0/(kBT) (phi-phiB) - (vRα + κα vR0)/(kBT) (p-pB)]
        expA = exp(- zA * e0 / (kBT) * (phi - phiB) - alphaA * (p - pB))
        expC = exp(- zC * e0 / (kBT) * (phi - phiB) - alphaC * (p - pB))
        yA = eyB_A * expA
        yC = eyB_C * expC
        # renormalization is not strictly necessary here because eq. (27) uses eyB and pressure,
        # but to keep numeric stability we optionally renormalize so ey0 >= 0:
        s = yA + yC
        # ensure ey0 = 1 - s (incompressibility implicit). If s > 0.9999, we allow it (?).
        return yA, yC
    end

    # second derivative helper (interior central, boundaries use one-sided)
    function d2_center(u, i)
        if i == 1
            # one-sided second derivative (forward)
            return (u[1] - 2.0 * u[2] + u[3]) / dx^2
        elseif i == N
            # backward difference
            return (u[N - 2] - 2.0 * u[N - 1] + u[N]) / dx^2
        else
            return (u[i - 1] - 2.0 * u[i] + u[i + 1]) / dx^2
        end
    end

    # central derivative first derivative (for momentum eq)
    function d1_center(u, i)
        if i == 1
            return (u[2] - u[1]) / dx
        elseif i == N
            return (u[N] - u[N - 1]) / dx
        else
            return (u[i + 1] - u[i - 1]) / (2.0 * dx)
        end
    end

    # Global residual function for NLsolve (autodiff-friendly)
    function residual!(F, U)
        # U is vector length tot
        # unpack into arrays
        phi = @view U[idx.(1:N, 1)]
        psi = @view U[idx.(1:N, 2)]
        p = @view U[idx.(1:N, 3)]
        yA = @view U[idx.(1:N, 4)]
        yC = @view U[idx.(1:N, 5)]

        # output F likewise
        # enforce boundary conditions and interior equations

        # 1) phi PDE (BSK): -ε d2φ + ε lc^2 d2ψ - q(φ,p) = 0  (we implement at interior nodes)
        for i in 1:N
            if i == 1
                # Dirichlet at electrode: phi = phi_left
                F[idx(i, 1)] = phi[i] - phi_left
            elseif i == N
                # bulk potential
                F[idx(i, 1)] = phi[i] - phi_right
            else
                d2phi = d2_center(phi, i)
                d2psi = d2_center(psi, i)
                # compute q at node i using current yA,yC
                qi = compute_q(phi[i], p[i], yA[i], yC[i])
                F[idx(i, 1)] = -ε * d2phi + ε * lc^2 * d2psi - qi
            end
        end

        # 2) psi relation: -d2φ + ψ = 0  (ψ = -Δφ) (impose on all nodes; use same boundary handling)
        for i in 1:N
            d2phi = d2_center(phi, i)
            F[idx(i, 2)] = -d2phi + psi[i]
        end

        # 3) momentum eq: ∂p/∂x + q ∂φ/∂x = 0
        for i in 1:N
            if i == N
                # Dirichlet pressure at bulk
                F[idx(i, 3)] = p[i] - pB
            else
                # use central/forward derivative for dp/dx and dφ/dx
                dpdx = d1_center(p, i)
                dphidx = d1_center(phi, i)
                qi = compute_q(phi[i], p[i], yA[i], yC[i])
                F[idx(i, 3)] = dpdx + qi * dphidx
            end
        end

        # 4) algebraic relations for y's (mass-balance → eq. (27))
        # enforce yA - eyB_A*exp[...] = 0, and similarly for yC
        for i in 1:N
            # compute predicted y from phi,p
            yA_pred, yC_pred = y_from_phi_p(phi[i], p[i])
            F[idx(i, 4)] = yA[i] - yA_pred
            F[idx(i, 5)] = yC[i] - yC_pred
        end
        return
    end
end


# ╔═╡ e4ecf2ab-0c9b-4680-8727-c8179a264df0
begin
    # Initial guess: use small potential ramp, zero psi, p = pB, and y's from bulk
    U0 = zeros(tot)
    for i in 1:N
        xi = X[i]
        # linear ramp potential from phi_left to phi_right
        phi_guess = phi_left * (1.0 - xi / L) + phi_right * (xi / L)
        U0[idx(i, 1)] = phi_guess
        U0[idx(i, 2)] = 0.0
        U0[idx(i, 3)] = pB
        # y from eq. (27) at bulk (approx)
        yAi, yCi = y_from_phi_p(phi_guess, pB)
        U0[idx(i, 4)] = yAi
        U0[idx(i, 5)] = yCi
    end

    # call NLsolve with autodiff Jacobian (ForwardDiff inside)
    println("Starting Newton solve...")
    sol = nlsolve(residual!, U0; method = :newton, autodiff = :forward, ftol = 1.0e-9, xtol = 1.0e-9, maxiters = 200)

    if sol.converged
        println("Converged in ", sol.iterations, " iterations.")
    else
        println("Solver did not converge; status: ", sol)
    end

    U = sol.zero
    phi = [U[idx(i, 1)] for i in 1:N]
    psi = [U[idx(i, 2)] for i in 1:N]
    p = [U[idx(i, 3)] for i in 1:N]
    yA = [U[idx(i, 4)] for i in 1:N]
    yC = [U[idx(i, 5)] for i in 1:N]
    y0 = [1.0 - yA[i] - yC[i] for i in 1:N]
    q = [compute_q(phi[i], p[i], yA[i], yC[i]) for i in 1:N]

    # Quick plots
    figure(figsize = (8, 10))

    subplot(4, 1, 1)
    plot(X .* 1.0e9, phi, "-k", lw = 2)
    ylabel("ϕ (V)")
    title("Potential")

    subplot(4, 1, 2)
    plot(X .* 1.0e9, p ./ (1.0e5), "-b", lw = 2)   # in bar-ish
    ylabel("p (bar)")
    title("Pressure")

    subplot(4, 1, 3)
    plot(X .* 1.0e9, yA, "-b", label = "y_A")
    plot(X .* 1.0e9, yC, "-r", label = "y_C")
    plot(X .* 1.0e9, y0, "-g", label = "y_0")
    legend()
    ylabel("ey (fraction)")
    title("Volume fractions (solvated ion fractions & free solvent)")

    subplot(4, 1, 4)
    plot(X .* 1.0e9, q, "-m", lw = 2)
    xlabel("x (nm)")
    ylabel("q (C/m^3)")
    title("Space charge density q")

    tight_layout()
    savefig("dreyer_mixture_solution.png", dpi = 300)
    println("Saved dreyer_mixture_solution.png")
end


# ╔═╡ Cell order:
# ╠═d1503e2a-76c8-4b40-9650-9121c9b771e7
# ╠═8871344f-1ded-4d21-86b6-8a61a7038aa1
# ╠═d5d26fae-c502-11f0-8122-b54486d2dc8d
# ╠═371f9c44-444d-4d10-a0d3-79befef1698d
# ╠═e4ecf2ab-0c9b-4680-8727-c8179a264df0
