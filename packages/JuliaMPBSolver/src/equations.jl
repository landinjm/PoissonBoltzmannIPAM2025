module equations

#placeholder params for now TODO: Implement params module elsewhere
struct PBParams{T<:Real}
  z::Vector{T}
  c_bulk::Vector{T}
  c0_bulk::Vector{T}
  cbar::T
  T::T
  F::T
end

"""
    molefractions!(phi, z, c_bulk, c0_bulk, Temp; f_model=false) 

Compute mole fraction at phi (potential), z (charges), c_bulk (bulk conc.),
c0_bulk (reference) and Temp (temperature)
f_model=false -> Poisson Boltzmann
f_model=true -> Bikerman

Adapting these from the simplecell bsk code for the PB equations implemented
in the notebook; have to create other modified equations from the discussion
for pressure etc
"""
function χ_from_E() #susceptibility calculation from applied Field
  #TODO: implementation
end

function molefractions!(
  phi,
  z::AbstractVector,
  c_bulk::AbstractVector,
  c0_bulk::AbstractVector,
  T;
  f_model::Bool = false,
)
  @assert length(z) == length(c_bulk) == length(c0_bulk)
  η = (z .* (phi .* F)) ./ (R * T)
  common_numerator = exp.((z .* (phi .* F)) ./ (R * T)) .* (c_bulk ./ c0_bulk)

  if !f_model
    model_denom = one(eltype(phi))  # 1 for Poisson-Boltzmann
  else
    # sum over species
    model_denom = one(eltype(phi)) .+ sum(common_numerator)  # for Bikerman
  end

  return common_numerator ./ model_denom
end

function spacecharge!(
  y::AbstractVector{T},
  ϕ::T,
  params::PBParams;
  f_model::Bool = false,
) where {T<:Real}
  # Call molefractions with unpacked parameters #TODO: check
  y .= molefractions!(
    ϕ,
    params.z,
    params.c_bulk,
    params.c0_bulk,
    params.T;
    f_model = f_model,
  )
  sumyz = zero(ϕ)
  for i in 1:length(params.z)
    sumyz += params.z[i] * y[i]
  end
  return params.F * params.cbar * sumyz
end

function pb_rhs(phi, z, c_bulk, cbar, c0_bulk, T)
  let y = similar(z)
    params = PBParams(z, c_bulk, c0_bulk, cbar, T, F)
    spacecharge!(y, phi, params; f_model = false)
  end
end

function bikerman_rhs(phi, z, c_bulk, cbar, c0_bulk, T)
  let y = similar(z)
    params = PBParams(z, c_bulk, c0_bulk, cbar, T, F)
    spacecharge!(y, phi, params; f_model = true)
  end
end

export χ_from_E, molefractions!, pb_rhs, bikerman_rhs

end
