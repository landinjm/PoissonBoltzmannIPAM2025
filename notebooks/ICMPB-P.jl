### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    using Revise
    using PlutoUI
    using VoronoiFVM
    using ExtendableGrids
    using LinearAlgebra
    using NLsolve
    using LessUnitful
    using LessUnitful.MoreUnitful
    using Test
    using PythonPlot
    using Colors
    using JuliaMPBSolver.ICMPBP: ICMPBData, ICMPBSystem, L_Debye, set_molarity!, dlcap0, update_derived!, apply_charge!, ysum
end

# ╔═╡ ef660f6f-9de3-4896-a65e-13c60df5de1e
md"""
# ICMPB with pressure, solvation and surface charges
"""

# ╔═╡ 4082c3d3-b728-4bcc-b480-cdee41d9ab99
# ╠═╡ skip_as_script = true
#=╠═╡
TableOfContents(title="",depth=5)
  ╠═╡ =#

# ╔═╡ 920b7d84-56c6-4958-aed9-fc67ba0c43f6
md"""
## 1. Intro

This code implements the model described in
[Müller, R., Fuhrmann, J., & Landstorfer, M. (2020). Modeling polycrystalline electrode-electrolyte interfaces: The differential capacitance. Journal of The Electrochemical Society, 167(10), 106512](https://iopscience.iop.org/article/10.1149/1945-7111/ab9cca/meta)


Equation numbers refer to the paper.


Concentrations are given in number densities, and calculations are done in mole fractions.
"""

# ╔═╡ 87ac16f4-a4fc-4205-8fb9-e5459517e1b8
md"""
If not stated otherwise, all calculations and calculation results are in coherent SI units.
"""

# ╔═╡ f4facb34-1f4a-432d-8a1e-30299e542bcd
begin
    const nm = ufac"nm"
    const V = ufac"V"
    const cm = ufac"cm"
    const μF = ufac"μF"
end

# ╔═╡ a3f23fe8-3b83-440a-8f4f-c4fedef5615b
L_Debye(ICMPBData(molarity = 0.01 * ph"N_A" / ufac"dm^3")) / nm

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
    cdl0=dlcap0(data)/ufac"μF/cm^2"
    @test cdl0 ≈ 22.84669184882525
end
  ╠═╡ =#

# ╔═╡ b1e333c0-cdaa-4242-b71d-b54ff71aef83
let
    data = ICMPBData()
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
    @test sumyz ≈ 0.0
    @test sumy ≈ 1.0
end

# ╔═╡ 97c5942c-8eb4-4b5c-8951-87ac0c9f396d
function c0_num!(c, φ, p, data)
    y = y0(p, data)
    sumyv = data.v0 * y
    for α in 1:(data.N)
        c[α] = y_α(φ, p, α, data)
        sumyv += c[α] * data.v[α]
    end
    return y / sumyv
end;

# ╔═╡ 00464966-2b1e-455c-a3a1-2af61c6649b7
dlcap_exact = 0.22846691848825248

# ╔═╡ 05334798-a072-41ae-b23e-f884baadb071
begin
    equidata = ICMPBData()
    set_molarity!(equidata, 0.01)
end

# ╔═╡ ddb3e60b-8571-465f-acf3-2403fb884363
@test dlcap0(equidata) ≈ dlcap_exact

# ╔═╡ a629e8a1-b1d7-42d8-8c17-43475785218e
begin

    L = 10nm
    n = 101

    X = range(0, L, length = n)

    grid = ExtendableGrids.simplexgrid(X)
    bfacemask!(grid, [L / 2], [L / 2], 3, tol = 1.0e-10 * nm)

end

# ╔═╡ 31a1f686-f0b6-430a-83af-187df411b293
sys_pp = ICMPBSystem(grid, equidata)

# ╔═╡ 14ac1c80-cae5-42f1-b0d3-33aa5bba4de6
begin
    sol0 = solve(sys_pp, verbose = "n")
    ysum(sys_pp, sol0)
end

# ╔═╡ efb12e12-825b-4dfd-aa10-c6afb304b6bf
ph"e" / ufac"nm^2"

# ╔═╡ 6f037b32-e2a8-4693-b46c-952d6b140e8e
begin
    apply_charge!(sys_pp, 2 * ph"e" / ufac"nm^2")
    sol1 = solve(sys_pp, verbose = "n", damp_initial = 0.1)
    ysum(sys_pp, sol1)
end

# ╔═╡ 1c0145d5-76b1-48c1-8852-de1a2668285a
molarities = [0.001, 0.01, 0.1, 1]

# ╔═╡ b1a69fe9-a3bd-4e52-95c3-efaa2d5f44c3
function qsweep(sys; qmax = 10, nsteps = 100)
    data = sys.physics.data
    (; ip, iφ) = data
    apply_charge!(sys, 0 * ph"e" / ufac"nm^2")
    sol = solve(sys, damp_initial = 0.1)

    volts = []
    Q = []
    for q in range(0, qmax, length = 50)
        @info q
        apply_charge!(sys, q * ph"e" / ufac"nm^2")
        sol = solve(sys, inival = sol)
        push!(volts, (sol[iφ, end] - sol[iφ, 1]) / 2)
        # Division by  comes in because the voltage we get here is the difference
        # between the electrodes and not the difference between electrode and bulk
        # which corresponds to the other standard dlcap experiment
        push!(Q, q * ph"e" / ufac"nm^2")
    end
    dlcaps = -(Q[2:end] - Q[1:(end - 1)]) ./ (volts[2:end] - volts[1:(end - 1)])
    return volts[1:(end - 1)], dlcaps
end

# ╔═╡ f1c33101-00e6-4af9-9e68-6cdf5fe92b59
qsweep(sys_pp)

# ╔═╡ 70e1a34b-9041-4151-91aa-4dd7907a5b13
function capscalc(sys, molarities)
    result = []
    for imol in 1:length(molarities)
        data = sys.physics.data
        set_molarity!(data, molarities[imol])
        update_derived!(data)

        t = @elapsed volts, caps = qsweep(sys; qmax = 5, nsteps = 100000)
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

# ╔═╡ a7f2692e-a15f-47b7-8486-8948ce7ab3f7
result_pp = capscalc(sys_pp, molarities)

# ╔═╡ e114ec0d-13d3-4455-b1c9-d1c5d76671d9
md"""
#### Pressure poisson problem
"""

# ╔═╡ 0b6f33b9-41d4-48fd-8026-8a3bddcc1989
md"""
#### Result plot

Compare with Fig 4.2 of [Fuhrmann (2015)](https://dx.doi.org/10.1016/j.cpc.2015.06.004) 

We see that for small voltages, the finite difference approach to calculate
things is quite inaccurate, this could be improved by a better control of the charge sweep. Otherwise, this result looks good.
"""

# ╔═╡ 000e87b5-cd1b-4b23-99be-6b7006502312
function capsplot(ax, result, title)
    hmol = 1 / length(result)
    ax.set_title(title)
    ax.grid()
    ax.set_xlabel("U/V")
    ax.set_ylabel(L"C_{dl}/(μF/cm^3)")

    for imol in 1:length(result)
        c = (1 - imol * hmol, 0, imol * hmol)

        ax.plot(
            result[imol].voltages, result[imol].dlcaps / (μF / cm^2),
            color = c, label = "$(result[imol].molarity)M"
        )
        ax.scatter(
            [0], [result[imol].cdl0] / (μF / cm^2),
            color = c
        )
    end
    return
end


# ╔═╡ 1859db0c-1c1e-4a78-9d4a-e2ec45c3cffc
let
    res = (600, 200)
    fig = figure(1, figsize = (6, 3))
    ax1 = fig.add_subplot(1, 1, 1)
    capsplot(ax1, result_pp, "Pressure poisson")
    tight_layout()
    matplotlib.pyplot.close()
    fig
end

# ╔═╡ Cell order:
# ╟─ef660f6f-9de3-4896-a65e-13c60df5de1e
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╟─4082c3d3-b728-4bcc-b480-cdee41d9ab99
# ╟─920b7d84-56c6-4958-aed9-fc67ba0c43f6
# ╟─87ac16f4-a4fc-4205-8fb9-e5459517e1b8
# ╠═f4facb34-1f4a-432d-8a1e-30299e542bcd
# ╠═a3f23fe8-3b83-440a-8f4f-c4fedef5615b
# ╟─5a210961-19fc-40be-a5f6-033a80f1414d
# ╠═fe704fb4-d07c-4591-b834-d6cf2f4f7075
# ╠═b1e333c0-cdaa-4242-b71d-b54ff71aef83
# ╟─97c5942c-8eb4-4b5c-8951-87ac0c9f396d
# ╠═00464966-2b1e-455c-a3a1-2af61c6649b7
# ╠═05334798-a072-41ae-b23e-f884baadb071
# ╠═ddb3e60b-8571-465f-acf3-2403fb884363
# ╠═a629e8a1-b1d7-42d8-8c17-43475785218e
# ╠═31a1f686-f0b6-430a-83af-187df411b293
# ╠═14ac1c80-cae5-42f1-b0d3-33aa5bba4de6
# ╠═efb12e12-825b-4dfd-aa10-c6afb304b6bf
# ╠═6f037b32-e2a8-4693-b46c-952d6b140e8e
# ╠═1c0145d5-76b1-48c1-8852-de1a2668285a
# ╠═b1a69fe9-a3bd-4e52-95c3-efaa2d5f44c3
# ╠═f1c33101-00e6-4af9-9e68-6cdf5fe92b59
# ╠═70e1a34b-9041-4151-91aa-4dd7907a5b13
# ╠═a7f2692e-a15f-47b7-8486-8948ce7ab3f7
# ╟─e114ec0d-13d3-4455-b1c9-d1c5d76671d9
# ╟─0b6f33b9-41d4-48fd-8026-8a3bddcc1989
# ╟─1859db0c-1c1e-4a78-9d4a-e2ec45c3cffc
# ╠═000e87b5-cd1b-4b23-99be-6b7006502312
