# LEAF_PHSYN_C4 - Maize (C4) Leaf Photosynthesis
# Paper https://doi.org/10.5194/egusphere-2025-1885


# ============ main function ============
function leaf_photosynthesis_c4(;
    leaf_temperature::Float64,      # Leaf temperature [K]
    wind_speed::Float64,             # Wind speed at 2m [m/s]
    specific_humidity::Float64,      # Specific humidity [kg/kg]
    pressure::Float64,               # Surface pressure [Pa]
    co2_ppm::Float64,                # Atmospheric CO2 [ppm]
    water_stress::Float64,           # Water stress factor [-]
    crop_height::Float64,            # Crop height [m]
    Vmax25::Float64,                 # Maximum Rubisco capacity at 25C [mol/m2/s]
    Qp::Float64,                     # Leaf absorbed PAR [W/m2]
    EFFCON::Float64,                 # Quantum efficiency [mol/mol]
    atheta::Float64,                 # Collatz coupling parameter (limit_Rubisco, limit_RuBP)
    btheta::Float64,                 # Collatz coupling parameter (limit_both, limit_PEP)
    m_H2O::Float64 = 4.0,            # Ball-Berry slope (H2O)
    b_H2O::Float64 = 0.04)           # Ball-Berry intercept (H2O)

    # ===== 1. Environmental variables =====
    CO2_atmosphere = co2_ppm / 1e6           # [mol/mol]
    H2O_internal = saturation_vapor_pressure(leaf_temperature, pressure) / ε_v  # saturated mole fraction
    H2O_ambient = specific_humidity / ε_v            # actual mole fraction
    Rh = min(1.0, specific_humidity / saturation_vapor_pressure(leaf_temperature, pressure))    # relative humidity [-]

    # ===== 2. Aerodynamic resistance (FAO56 P20) =====
    d0 = crop_height * 2.0 / 3.0                    # zero plane displacement
    z0m = crop_height * 0.123                       # roughness length for momentum
    z0h = 0.1 * z0m                         # roughness length for heat
    Rb = log((2.0 - d0) / z0m) * log((2.0 - d0) / z0h) / (karman_constant^2 * wind_speed)  # [s/m]
    GB_leaf = 1.0 / Rb                     # leaf boundary conductance [m/s]
    GB_H2O = GB_leaf * pressure / (leaf_temperature * R_vap * M_H2O)  # [mol(H2O)/m2/s]
    GB_CO2 = GB_H2O / 1.4                  # [mol(CO2)/m2/s]

    # ===== 3. Dark respiration (eq.3, Bonan 2011) =====
    tpem = 2.0^((leaf_temperature - 298.15) / 10)
    RSP_leaf = 0.025 * Vmax25 * tpem / (1 + exp(1.3 * (leaf_temperature - 328.15)))  # [mol/m2/s]

    # ===== 4. Actual carboxylation rate (eq.9-11) =====
    factor_high_temp = 1 + exp(0.3 * (leaf_temperature - 313.15))
    factor_low_temp = 1 + exp(0.2 * (288.15 - leaf_temperature))
    Vm = water_stress * Vmax25 * tpem / (factor_high_temp * factor_low_temp)

    # ===== 5. Rubisco / RUBP limited photosynthesis (eq.6-7) =====
    limit_Rubisco_Ac = Vm   # Rubisco-limited photosynthesis [mol/m2/s]
    limit_RuBP_Aj = EFFCON * 4.6 * Qp * 1e-6  # convert W/m2 to mol/m2/s  # [mol/m2/s]
    # limit_Rubisco_Ac + limit_RuBP_Aj -> limit_both_Ai (eq.4)
    square_root = max((limit_RuBP_Aj + limit_Rubisco_Ac)^2 - 4 * atheta * limit_RuBP_Aj * limit_Rubisco_Ac, 0.0)
    limit_both_Ai = ((limit_RuBP_Aj + limit_Rubisco_Ac) - sqrt(square_root)) / (2 * atheta)  # [mol/m2/s]

    # ===== 6. CO2 limitation quadratic equation solving (eq.20) =====
    kp = Vmax25 > 0 ? Vmax25 * 20000.0 * tpem : 0.7  # Collatz 1992

    G1_CO2 = m_H2O / 1.6
    G0_CO2 = b_H2O / 1.6

    aa = GB_CO2^2 * G1_CO2 * Rh - GB_CO2 * G0_CO2 - kp * (G0_CO2 - GB_CO2 * G1_CO2 * Rh + GB_CO2)
    bb = CO2_atmosphere * GB_CO2^2 * G0_CO2 - GB_CO2 * G0_CO2 * RSP_leaf + GB_CO2^2 * G1_CO2 * Rh * RSP_leaf -
         kp * CO2_atmosphere * (GB_CO2^2 * G1_CO2 * Rh - 2 * GB_CO2 * G0_CO2 - GB_CO2^2)
    cc = CO2_atmosphere * GB_CO2^2 * G0_CO2 * (RSP_leaf - kp * CO2_atmosphere)
    square_root = max(bb^2 - 4 * aa * cc, 0.0)

    # three analytical solutions
    NetAssim_solution_1 = (-bb - sqrt(square_root)) / (2 * aa)
    NetAssim_solution_2 = (-bb + sqrt(square_root)) / (2 * aa)
    NetAssim_solution_3 = (kp * CO2_atmosphere - RSP_leaf) / (1 + kp * (1 / GB_CO2 + 1 / G0_CO2))

    # select valid solution (last one wins, following Fortran logic)
    NetAssim_leaf = 0.0
    for (NetAssim_solution, cond_ge_0) in [(NetAssim_solution_1, true), (NetAssim_solution_2, true), (NetAssim_solution_3, false)]
        valid_init = cond_ge_0 ? (NetAssim_solution >= 0) : (NetAssim_solution < 0)
        if valid_init
            CO2_surface = CO2_atmosphere - NetAssim_solution / GB_CO2
            GS_CO2_candidates = G0_CO2 + G1_CO2 * Rh * NetAssim_solution / CO2_surface
            CO2_intercell_candidates = (NetAssim_solution + RSP_leaf) / kp
            if GS_CO2_candidates > 0 && CO2_intercell_candidates > 0
                NetAssim_leaf = NetAssim_solution
            end
        end
    end
    limit_PEP_Ap = NetAssim_leaf + RSP_leaf  # PEP-limited photosynthesis (eq.2)

    # limit_both_Ai + limit_PEP_Ap -> GPP (eq.5)
    square_root = max((limit_both_Ai + limit_PEP_Ap)^2 - 4 * btheta * limit_both_Ai * limit_PEP_Ap, 0.0)
    GPP_leaf = max(((limit_PEP_Ap + limit_both_Ai) - sqrt(square_root)) / (2 * btheta), 0.0)  # [mol/m2/s]

    # ===== 7. Stomatal conductance (eq.2) =====
    NetAssim_leaf_final = GPP_leaf - RSP_leaf
    GS_CO2 = NetAssim_leaf_final > 0 ? G1_CO2 * Rh * NetAssim_leaf_final / (CO2_atmosphere - NetAssim_leaf_final / GB_CO2) + G0_CO2 : G0_CO2
    GS_H2O = GS_CO2 * 1.6

    # ===== 8. Transpiration =====
    H2O_surface = min(H2O_internal, (GB_H2O * H2O_ambient + GS_H2O * H2O_internal) / (GB_H2O + GS_H2O))
    TSP_leaf = GS_H2O * (H2O_internal - H2O_surface) * M_H2O  # [kg/m2/s]

    return (gpp=GPP_leaf, rsp=RSP_leaf, tsp=TSP_leaf)
end