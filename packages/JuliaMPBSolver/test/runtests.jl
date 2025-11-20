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
    @testset "mpbsolve" begin
        X, C0, CP, CM = mpbpsolve(n = 5)

        X_ref = 0.0:2.5e-11:1.0e-10
        C0_ref = [53.07684745281815, 53.239200914598854, 53.292594052218185, 53.239200914598854, 53.07684745281816]
        CP_ref = [0.15690926770216584, 0.12582063392563822, 0.10070027035371927, 0.08043382838355795, 0.06410460022345751]
        CM_ref = [0.06410460022345683, 0.08043382838355714, 0.10070027035371823, 0.1258206339256369, 0.15690926770216423]
        @test X ≈ X_ref
        @test C0 ≈ C0_ref
        @test CP ≈ CP_ref
        @test CM ≈ CM_ref
    end
end
