module DataOut

using HDF5
using CSV
using DataFrames
using ExtendableGrids

using ..Grid

function write_hdf5_data!(
        base_filename::String,
        grid::ExtendableGrid,
        solution,
        solution_name::String = "solution",
    )
    @assert dim_grid(grid) == 1

    file_stream = h5open(base_filename, "cw")

    grid_coordinates = Grid.get_coordinates(grid)

    @assert length(grid_coordinates) == length(solution)

    if haskey(file_stream, "x")
        existing_x = read(file_stream["x"])
        if existing_x != grid_coordinates
            close(file_stream)
            error(
                "Grid coordinates mismatch: existing x-coordinates differ from current grid",
            )
        end
    else
        file_stream["x"] = grid_coordinates
    end

    file_stream[solution_name] = solution

    close(file_stream)
    return nothing
end

function read_hdf5_data(filename::String, solution_name::String = "solution")
    return h5open(filename, "r") do file_stream
        grid_coordinates = read(file_stream["x"])
        solution = read(file_stream[solution_name])
        return grid_coordinates, solution
    end
end

function write_csv_data!(filename::String, solution)
    df = DataFrame(value = solution)
    CSV.write(filename, df)

    return nothing
end

export write_data!, read_hdf5_data

end
