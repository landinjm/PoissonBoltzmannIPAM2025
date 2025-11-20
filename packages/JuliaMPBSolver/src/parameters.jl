module Parameters

using ..Units

struct UserParameters{F <: AbstractFloat, I <: Integer}
    temperature::F
    dielectric_susceptibility::F
    reference_pressure::F
    reference_electric_potential::F
    charge_numbers::AbstractVector{<:I}
    total_concentration::F
    bulk_ion_concentrations::AbstractVector{<:F}
    boundary_electron_density::F
    electric_susceptibility_decrement_parameter::F
    use_bikerman_model::Bool
end

struct ComputedParameters{F <: AbstractFloat}
    debye_length::F
    double_layer_capacitance::F
    bulk_solvent_concentration::F
end

function UserParameters(;
        temperature::AbstractFloat,
        dielectric_susceptibility::AbstractFloat,
        reference_pressure::AbstractFloat,
        reference_electric_potential::AbstractFloat,
        charge_numbers::AbstractVector{<:Integer},
        total_concentration::AbstractFloat,
        bulk_ion_concentrations::AbstractVector{<:AbstractFloat},
        boundary_electron_density::AbstractFloat,
        electric_susceptibility_decrement_parameter::AbstractFloat,
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
        boundary_electron_density,
        electric_susceptibility_decrement_parameter,
        use_bikerman_model,
    )
end

function ComputedParameters(user_parameters::UserParameters)
    debye_length = sqrt(
        (1 + user_parameters.dielectric_susceptibility) *
            Units.vacuum_permittivity *
            Units.thermal_energy(user_parameters.temperature) /
            (Units.F^2 * user_parameters.bulk_ion_concentrations[1]),
    )
    double_layer_capacitance = sqrt(
        2 *
            (1 + user_parameters.dielectric_susceptibility) *
            Units.vacuum_permittivity *
            Units.F^2 *
            user_parameters.bulk_ion_concentrations[1] /
            Units.thermal_energy(user_parameters.temperature),
    )
    bulk_solvent_concentration =
        user_parameters.total_concentration -
        sum(user_parameters.bulk_ion_concentrations)
    return ComputedParameters(
        debye_length,
        double_layer_capacitance,
        bulk_solvent_concentration,
    )
end

export UserParameters, ComputedParameters

end
