module Grid

using ExtendableGrids

struct GeometricGrid
  domain_size::Float32
  refinement::Int32
  hmin::Float32
  hmax::Float32
  use_offset::Bool
end

GeometricGrid() = GeometricGrid(10.0f0, 1, 0.1f0, 1.0f0, false)

function create_grid(grid::GeometricGrid)
  # Compute a minimum and maximum cell size from the provided values and
  # refinement. This is a little weird to me because the GeometricGrid
  # hmin and hmax should be obeyed here (logically). However, then
  # the refinement is of little use. 
  # TODO: IDK Fix this logic
  local_hmin = grid.hmin / 2.0^grid.refinement
  local_hmax = grid.hmax / 2.0^grid.refinement

  # Create a little offset 
  offset = 0.0f0
  if grid.use_offset
    offset = 1.0e-3 * grid.domain_size
  end

  # Create to geometric spacings and glue them together. This way
  # the coarsest cells are in the middle
  x_left = geomspace(offset, grid.domain_size / 2.0, local_hmin, local_hmax)
  x_right =
    geomspace(grid.domain_size / 2.0, grid.domain_size, local_hmax, local_hmin)
  x = glue(x_left, x_right)

  # Create the simplex grid
  x = simplexgrid(x)

  # Add a face to the grid in the middle
  # TODO: Refactor this as a function. We don't always need this additional 
  # face to impose boundary conditions
  bfacemask!(
    x,
    [grid.domain_size / 2],
    [grid.domain_size / 2],
    3,
    tol = 1.0e-2 * local_hmin,
  )

  return x
end

function is_equivalent(
  grid::ExtendableGrid,
  dim::Int,
  n_nodes::Int,
  n_cells::Int,
  n_boundary_faces::Int,
)
  return dim_grid(grid) == dim &&
         num_nodes(grid) == n_nodes &&
         num_cells(grid) == n_cells &&
         num_bfaces(grid) == n_boundary_faces
end

export GeometricGrid, create_grid, is_equivalent

end