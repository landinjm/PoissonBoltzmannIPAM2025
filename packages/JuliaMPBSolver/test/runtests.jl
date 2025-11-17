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
  @testset "DataOut" begin
    grid = JuliaMPBSolver.Grid.create_full_cell(
      JuliaMPBSolver.Grid.GeometricGrid(10.0f0, 0, 1.0f0, 1.0f0, false),
    )
    coordinates = JuliaMPBSolver.Grid.get_coordinates(grid)
    solution = rand(length(coordinates))
    filename = "test_output.h5"
    JuliaMPBSolver.DataOut.write_hdf5_data!(filename, grid, solution)
    @test isfile(filename)
    JuliaMPBSolver.DataOut.read_hdf5_data(filename) == (coordinates, solution)
    rm(filename)

    filename = "test_output.csv"
    JuliaMPBSolver.DataOut.write_csv_data!(filename, solution)
    @test isfile(filename)
    rm(filename)
  end
end
