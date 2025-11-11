using JuliaMPBSolver
using Test

@testset "JuliaMPBSolver" begin
  @testset "Units" begin
    @test JuliaMPBSolver.Units.RT(0) == 0
  end
  @testset "Grids" begin
    grid = JuliaMPBSolver.Grid.create_full_cell(
      JuliaMPBSolver.Grid.GeometricGrid(10.0f0, 0, 1.0f0, 1.0f0, false),
    )
    @test JuliaMPBSolver.Grid.is_equivalent(grid, 1, 11, 10, 3)

    grid = JuliaMPBSolver.Grid.create_half_cell(
      JuliaMPBSolver.Grid.GeometricGrid(10.0f0, 0, 1.0f0, 1.0f0, false),
    )
    @test JuliaMPBSolver.Grid.is_equivalent(grid, 1, 11, 10, 2)
  end
end