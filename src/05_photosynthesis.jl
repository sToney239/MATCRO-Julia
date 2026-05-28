# PHSYN - Canopy Photosynthesis Wrapper
# Dispatches to C3 or C4 leaf-level photosynthesis for sunlit/shade leaves,
# then aggregates by LAI. Mirrors Fortran SUB_PHSYN.f90 structure.
# Paper https://doi.org/10.5194/egusphere-2025-1885

include("05_1_photosynthesis_C3.jl")
include("05_2_photosynthesis_C4.jl")

# ============ main function ============
function calc_photosynthesis(;
    Qp_sunlit::Float64,           # Absorbed PAR per unit sunlit leaf area [W/m²]
    Qp_shade::Float64,            # Absorbed PAR per unit shade leaf area [W/m²]
    Vmax25_sunlit::Float64,       # Vmax25 per unit sunlit leaf area [mol/m²/s]
    Vmax25_shade::Float64,        # Vmax25 per unit shade leaf area [mol/m²/s]
    LAI_sunlit::Float64,          # Sunlit leaf area index [m²(leaf)/m²]
    LAI_shade::Float64,           # Shade leaf area index [m²(leaf)/m²]
    leaf_temperature::Float64,    # Leaf temperature (= air temperature) [K]
    wind_speed::Float64,          # Wind speed at 2m [m/s]
    specific_humidity::Float64,   # Specific humidity [kg/kg]
    pressure::Float64,            # Surface pressure [Pa]
    co2_ppm::Float64,             # Atmospheric CO2 [ppm]
    water_stress::Float64,        # Water stress factor [-]
    crop_height::Float64,         # Crop height [m]
    EFFCON::Float64,              # Quantum efficiency [mol/mol]
    atheta::Float64,              # Collatz coupling parameter (Rubisco vs RuBP)
    btheta::Float64,              # Collatz coupling parameter (co-limitation vs TPU/PEP)
    m_H2O::Float64 = 4.0,        # Ball-Berry slope (H2O)
    b_H2O::Float64 = 0.04,       # Ball-Berry intercept (H2O)
    crop_name::String = "Soybeans")  # Crop name: "Rice","Wheat","Soybeans","Maize"

    LAI_total = LAI_sunlit + LAI_shade

    if LAI_total <= 0.0
        return (gpp=0.0, rsp=0.0, tsp=0.0)
    end

    # Select leaf-level function based on crop type
    photosynthesis_function = crop_name == "Maize" ? leaf_photosynthesis_c4 : leaf_photosynthesis_c3

    # ----- Sunlit leaves -----
    # both C3 & C4 leaves us Vmax25_sunlit
    if LAI_sunlit > 0.0
        r_sun = photosynthesis_function(;
            leaf_temperature=leaf_temperature, wind_speed=wind_speed,
            specific_humidity=specific_humidity, pressure=pressure,
            co2_ppm=co2_ppm, water_stress=water_stress,
            crop_height=crop_height, Vmax25=Vmax25_sunlit, Qp=Qp_sunlit,
            EFFCON=EFFCON, atheta=atheta, btheta=btheta,
            m_H2O=m_H2O, b_H2O=b_H2O)
    else
        r_sun = (gpp=0.0, rsp=0.0, tsp=0.0)
    end

    # ----- Shade leaves -----
    # C3: shade leaves use Vmax25_shade
    # C4: shade leaves use Vmax25_sunlit
    Vmax25_shade_eff = crop_name == "Maize" ? Vmax25_sunlit : Vmax25_shade
    if LAI_shade > 0.0
        r_shade = photosynthesis_function(;
            leaf_temperature=leaf_temperature, wind_speed=wind_speed,
            specific_humidity=specific_humidity, pressure=pressure,
            co2_ppm=co2_ppm, water_stress=water_stress,
            crop_height=crop_height, Vmax25=Vmax25_shade_eff, Qp=Qp_shade,
            EFFCON=EFFCON, atheta=atheta, btheta=btheta,
            m_H2O=m_H2O, b_H2O=b_H2O)
    else
        r_shade = (gpp=0.0, rsp=0.0, tsp=0.0)
    end

    # ----- Aggregate by LAI -----
    GPP = r_sun.gpp * LAI_sunlit + r_shade.gpp * LAI_shade  # [mol/m²(ground)/s]
    RSP = r_sun.rsp * LAI_sunlit + r_shade.rsp * LAI_shade
    TSP = r_sun.tsp * LAI_sunlit + r_shade.tsp * LAI_shade  # [kg/m²(ground)/s]

    return (gpp=GPP, rsp=RSP, tsp=TSP)
end
