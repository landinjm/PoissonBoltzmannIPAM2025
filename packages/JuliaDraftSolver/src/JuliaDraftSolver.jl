module JuliaDraftSolver

using VoronoiFVM
using ExtendableGrids

"""
    solve_nonlinear_poisson(;
                            n=21,
                            domain=[0,1],
                            bcval=[1,1],
                            source=1,
                            m=2,
                            inival=sum(bcval)/2)

Solve ``\\Delta u^m = f``.

Keyword arguments:
- `n`: Number of unknowns
- `domain`: Domain boundaries
- `bcval`: Dirichlet boundary values
- `source`: constant source term
- `m`: nonlinearity exponent 
- `inival`: initial value

Returns: vector of grid nodes `X` and vector of solution values `U`.
"""
function solve_nonlinear_poisson(;
        n = 21,
        domain = [0, 1],
        bcval = [1, 1],
        source = 1,
        m = 2,
        inival = sum(bcval) / 2
    )
    mybcval = Float64[bcval...]

    X = collect(range(domain..., length = n))
    grid = simplexgrid(X)

    function flux(y, u, edge, data)
        return y[1] = u[1, 1]^m - u[1, 2]^m
    end

    function bcondition(y, u, bnode, data)
        boundary_dirichlet!(y, u, bnode; region = 1, value = mybcval[1])
        return boundary_dirichlet!(y, u, bnode; region = 2, value = mybcval[2])
    end

    function src(y, node, data)
        return y[1] = source
    end

    sys = VoronoiFVM.System(grid; inival, flux, bcondition, source = src, species = [1])
    sol = solve(sys; inival)
    return X, Vector(sol[1, :])
end


export solve_nonlinear_poisson
end # module JuliaSolver
