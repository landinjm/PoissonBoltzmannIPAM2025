module Units

using LessUnitful, Unitful

#
# A little foreward about LessUnitful.jl. If you decide to look up 
# the package you can find information there. In short, it's a small
# wrapper around Unitful.jl where the units are "unitless". This is
# an important distinction because our FVM and linear algebra 
# implementations cannot be made unit-aware. 
#

# Faraday constant
const F = ph"N_A" * ph"e"

# Temperature unit
const K = ufac"K"

# Distance units
const nm = ufac"nm"
const m = ufac"m"
const dm = ufac"dm"
const angstrom = ufac"Å"

# Voltage unit
const V = ufac"V"

# Moles
const mol = ufac"mol"

# Thermal energy
thermal_energy(temperature) = ph"R" * temperature
RT(x) = thermal_energy(x)

# Export all constants and functions
export F, K, nm, m, dm, angstrom, V, mol, thermal_energy, RT

end
