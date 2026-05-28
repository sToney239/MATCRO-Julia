# constants.jl - Shared constants (from Fortran param.inc)

const seconds_per_day = 86400.0 # how many seconds in a day




# Physics
const karman_constant = 0.4        # von Karman constant
const M_H2O = 0.018015      # molecular weight of water [kg/mol]
const R_air = 287.04         # gas constant of dry air [J/(kg·K)]
const R_vap = 461.5          # gas constant of water vapor [J/(kg·K)]
const ε_v = R_air / R_vap
const es0 = 611.0           # reference saturation vapor pressure [Pa]
const L_vaporization = 2.5e6  # latent heat of vaporization [J/kg]
const L_melt = 3.4e5        # latent heat of melting [J/kg]
const T_ice = 273.15        # ice point temperature [K]
const g0 = 9.8              # gravitational acceleration [m/s²]
const ρ_water = 1000.0      # density of water [kg/m³]

# ============ shared functions ============
# saturation_vapor_pressure: saturation specific humidity
function saturation_vapor_pressure(T::Float64, P::Float64)::Float64
    sgn = T >= T_ice ? 1.0 : -1.0
    return ε_v * es0 / P * exp((L_vaporization + L_melt / 2 * (1 - sgn)) / R_vap * (1 / T_ice - 1 / T))
end
