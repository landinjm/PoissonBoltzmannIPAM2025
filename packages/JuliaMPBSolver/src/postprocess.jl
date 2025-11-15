module Postprocess

using ..Parameters
using ..Units

function compute_mole_fractions!(
  mole_fractions::AbstractVector,
  electric_potential,
  user_parameters::UserParameters,
)
  n_species = length(user_parameters.charge_numbers)

  @assert length(mole_fractions) == n_species

  for i in 1:n_species
    mole_fractions[i] =
      exp(
        -user_parameters.charge_numbers[i] * electric_potential * F /
        RT(user_parameters.temperature),
      ) * user_parameters.bulk_ion_concentrations[i] /
      user_parameters.bulk_solvent_concentration
  end

  if user_parameters.use_bikerman_model
    mole_fractions ./= one(electric_potential) + sum(mole_fractions)
  end

  return nothing
end

function compute_concentrations!(
  concentrations::Matrix{<:AbstractFloat},
  electric_potential,
  user_parameters::UserParameters,
)
  n_species = length(user_parameters.charge_numbers)
  n_nodes = length(electric_potential)

  @assert size(concentrations, 1) == n_species
  @assert size(concentrations, 2) == n_nodes

  point_wise_concentrations = zeros(n_species)

  for i in 1:n_nodes
    compute_mole_fractions!(
      point_wise_concentrations,
      electric_potential[i],
      user_parameters,
    )
    for j in 1:n_species
      concentrations[j, i] =
        user_parameters.total_concentration * point_wise_concentrations[j]
    end
  end

  return nothing
end

function compute_concentrations(
  electric_potential,
  user_parameters::UserParameters,
)
  n_species = length(user_parameters.charge_numbers)
  n_nodes = length(electric_potential)

  concentrations = zeros(n_species, n_nodes)

  compute_concentrations!(concentrations, electric_potential, user_parameters)

  return concentrations
end

function compute_spacecharge(
  mole_fractions,
  electric_potential,
  user_parameters::UserParameters,
)
  # TODO: This function is slow because of the dynamic allocation of arrays.
  # This won't fly with 2D and 3D problems. We need to pre-allocate the arrays
  # and pass them as arguments.
  n_species = length(user_parameters.charge_numbers)

  compute_mole_fractions!(mole_fractions, electric_potential, user_parameters)

  space_charge = zero(electric_potential)

  for i in 1:n_species
    space_charge +=
      Units.F *
      user_parameters.total_concentration *
      user_parameters.charge_numbers[i] *
      mole_fractions[i]
  end

  return space_charge
end

function compute_gradient(y::AbstractVector, x::AbstractVector)
  # Compute and return the gradient with a forward euler finite difference scheme.
  # Note that this function does not take into account boundary conditions and returns
  # an array of length n - 1, where n is the length of x. 
  @assert length(y) == length(x)

  return (y[2:end] - y[1:(end-1)]) ./ (x[2:end] - x[1:(end-1)])
end

export compute_mole_fractions!,
  compute_concentrations!, compute_concentrations, compute_spacecharge

end
