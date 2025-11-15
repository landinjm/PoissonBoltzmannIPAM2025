module Parameters

struct UserParameters{F<:AbstractFloat,I<:Integer}
  temperature::F
  dielectric_susceptibility::F
  reference_pressure::F
  reference_electric_potential::F
  charge_numbers::AbstractVector{<:I}
  total_concentration::F
  bulk_ion_concentrations::AbstractVector{<:F}
  bulk_solvent_concentration::F
  use_bikerman_model::Bool
end

function UserParameters(;
  temperature::AbstractFloat,
  dielectric_susceptibility::AbstractFloat,
  reference_pressure::AbstractFloat,
  reference_electric_potential::AbstractFloat,
  charge_numbers::AbstractVector{<:Integer},
  total_concentration::AbstractFloat,
  bulk_ion_concentrations::AbstractVector{<:AbstractFloat},
  bulk_solvent_concentration::AbstractFloat,
  use_bikerman_model::Bool,
)
  @assert dot(bulk_ion_concentrations, charge_numbers) == 0
  return UserParameters(
    temperature,
    dielectric_susceptibility,
    reference_pressure,
    reference_electric_potential,
    charge_numbers,
    total_concentration,
    bulk_ion_concentrations,
    bulk_solvent_concentration,
    use_bikerman_model,
  )
end

export UserParameters

end
