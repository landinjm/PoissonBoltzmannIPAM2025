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

  file_stream["x"] = grid_coordinates
  file_stream[solution_name] = solution

  close(file_stream)
  return nothing
end

function read_hdf5_data(filename::String, solution_name::String = "solution")
  file_stream = h5open(filename, "r")

  grid_coordinates = file_stream["x"][:]
  solution = file_stream[solution_name][:]

  close(file_stream)

  return grid_coordinates, solution
end

function write_csv_data!(filename::String, solution)
  df = DataFrame(value = solution)
  CSV.write(filename, df)

  return nothing
end

export write_data!, read_hdf5_data

end
