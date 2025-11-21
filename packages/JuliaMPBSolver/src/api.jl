using ExtendableGrids
using LessUnitful
using VoronoiFVM
function mpbpsolve(;
        n = 21,
        chargenumbers = [-1, 1],
        domain = [0, 1.0e-10],
        surfacecharge = [0.16, -0.16],
        bulkmolarity = 0.1
    )
    X = range(domain..., length = n)
    grid = ExtendableGrids.simplexgrid(X)
    xmid = (domain[1] + domain[2]) / 2
    bfacemask!(grid, [xmid], [xmid], 3, tol = 1.0e-10 * ufac"nm")
    data = ICMPBP.ICMPBData(; z = chargenumbers, q = surfacecharge)
    ICMPBP.set_molarity!(data, bulkmolarity)
    sys = ICMPBP.ICMPBSystem(grid, data)
    sol = solve(sys; inival = 0.01, damp_initial = 0.05, verbose = "n")
    c0 = ICMPBP.calc_c0mol(sol, sys)
    c = ICMPBP.calc_cmol(sol, sys)
    return X, c0, c[1, :], c[2, :]
end

function icmpbpsolve(;
        n = 21,
        chargenumbers = [-1, 1],
        domain = [0, 1.0e-10],
        surfacecharge = [0.16, -0.16],
        averagemolarity = 1
    )
    @info averagemolarity
    if iseven(n)
        n = n + 1
    end
    X = range(domain..., length = n)
    grid = ExtendableGrids.simplexgrid(X)
    xmid = (domain[1] + domain[2]) / 2
    bfacemask!(grid, [xmid], [xmid], 3, tol = 1.0e-10 * ufac"nm")
    data = ICMPBP.ICMPBData(; z = chargenumbers, q = surfacecharge, conserveions = true)
    ICMPBP.set_molarity!(data, averagemolarity)
    sys = ICMPBP.ICMPBSystem(grid, data)
    sol = solve(sys; inival = unknowns(sys, data), damp_initial = 0.05, verbose = "n")
    c0 = ICMPBP.calc_c0mol(sol, sys)
    c = ICMPBP.calc_cmol(sol, sys)
    return X, c0, c[1, :], c[2, :]
end
