# LEAF_PHSYN_C3 - C3 Crop Leaf Photosynthesis (Rice, Wheat, Soybeans)
# Based on: Masutomi et al. (2016), Bernacchi et al. (2001, 2003), Farquhar et al. (1980)
# CO2 down-regulation: R_JV factor scales Jmax with CO2
# Paper https://doi.org/10.5194/egusphere-2025-1885

# ============ main function ============
function leaf_photosynthesis_c3(;
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
    atheta::Float64,                 # Collatz coupling parameter (Rubisco vs RuBP)
    btheta::Float64,                 # Collatz coupling parameter (co-limitation vs TPU)
    m_H2O::Float64 = 4.0,            # Ball-Berry slope (H2O)
    b_H2O::Float64 = 0.04)           # Ball-Berry intercept (H2O)

    # ===== 1. Absorbed PAR =====
    PAR_absorbed = 4.6 * Qp  # [μmol/m2/s]

    # ===== 2. Enzyme kinetics (Bernacchi et al. 2001, 2003) =====
    Kc = exp(38.05 - 79430.0 / (leaf_temperature * R_vap * M_H2O))     # [μmol/mol]
    Ko = exp(20.30 - 36380.0 / (leaf_temperature * R_vap * M_H2O)) / 1000.0  # [mol/mol]
    Gamma_star = exp(19.02 - 37830.0 / (leaf_temperature * R_vap * M_H2O))   # [μmol/mol]

    # ===== 3. Vcmax and Jmax (temperature-dependent) =====
    Vm = Vmax25 * 1e6 * exp(26.35 - 65330.0 / (leaf_temperature * R_vap * M_H2O)) * water_stress  # [μmol/m2/s]
    R_JV = 1.67 * (0.941 + 1.32e-4 * co2_ppm) / (0.941 + 1.32e-4 * 368.87)   # CO2 down-regulation
    Jm = 1.67 * Vmax25 * 1e6 * exp(17.70 - 43900.0 / (leaf_temperature * R_vap * M_H2O)) * R_JV   # [μmol/m2/s]

    # ===== 4. Dark respiration =====
    RSP_leaf = 0.015 * Vmax25 * 1e6 * exp(18.72 - 46390.0 / (leaf_temperature * R_vap * M_H2O))  # [μmol/m2/s]

    # ===== 5. Kc*(1 + O2/Ko) - Michaelis-Menten composite =====
    PO2 = 20900.0  # Atmospheric O2 partial pressure [Pa]
    RRKK = Kc * (1.0 + (PO2 / pressure) / Ko)  # [μmol/mol]

    # ===== 6. Aerodynamic resistance (FAO56 P20) =====
    d0 = crop_height * 2.0 / 3.0                    # zero plane displacement
    z0m = crop_height * 0.123                       # roughness length for momentum
    z0h = 0.1 * z0m                         # roughness length for heat
    Rb = log((2.0 - d0) / z0m) * log((2.0 - d0) / z0h) / (karman_constant^2 * wind_speed)  # [s/m]
    GB_leaf = 1.0 / Rb                     # leaf boundary conductance [m/s]
    GB_H2O = GB_leaf * pressure / (leaf_temperature * R_vap * M_H2O)  # [mol(H2O)/m2/s]
    GB_CO2 = GB_H2O / 1.4                  # [mol(CO2)/m2/s]

    # ===== 7. Electron transport rate (hyperbolic light response) =====
    e1 = 0.7
    e2 = -(PAR_absorbed * EFFCON + Jm)
    e3 = PAR_absorbed * EFFCON * Jm
    square_root = max(e2^2 - 4.0 * e1 * e3, 0.0)
    Je = (-e2 - sqrt(square_root)) / (2.0 * e1)  # [μmol/m2/s]

    # ===== 8. Environmental variables =====
    CO2_atmosphere = co2_ppm  # [μmol/mol = ppm]
    H2O_internal = saturation_vapor_pressure(leaf_temperature, pressure) / ε_v  # saturated mole fraction
    H2O_ambient = specific_humidity / ε_v            # actual mole fraction
    Rh = min(1.0, specific_humidity / saturation_vapor_pressure(leaf_temperature, pressure))    # relative humidity [-]

    # ===== 9. Solve for Rubisco-limited (Ac) and RuBP-limited (Aj) rates =====
    # Baldocchi (1994) approach: solve coupled A-ci-stomata equations
    aa = [Vm, Je]
    bb = [RRKK, 8.0 * Gamma_star]
    dd = [Gamma_star, Gamma_star]
    ee = [1.0, 4.0]

    G1_CO2 = m_H2O / 1.6
    G0_CO2 = b_H2O / 1.6

    limit_Rubisco_Ac = 0.0  # Rubisco-limited GPP [μmol/m2/s]
    limit_RuBP_Aj = 0.0  # RuBP-limited GPP [μmol/m2/s]

    for i in 1:2
        α = GB_CO2 * CO2_atmosphere
        β = G1_CO2 * GB_CO2 * Rh - G0_CO2
        γ = aa[i] * dd[i] + bb[i] * RSP_leaf
        ζ = aa[i] - ee[i] * RSP_leaf
        ω = CO2_atmosphere * GB_CO2 * ζ - γ * GB_CO2
        ψ = ee[i] * CO2_atmosphere * GB_CO2 + ζ + bb[i] * GB_CO2

        a1 = ee[i] * β - ee[i] * GB_CO2
        a2 = ee[i] * G0_CO2 * α - β * ψ + ee[i] * GB_CO2 * α + GB_CO2 * ζ
        a3 = -G0_CO2 * α * ψ + β * ω - GB_CO2 * α * ζ
        a4 = G0_CO2 * α * ω

        # Normalize to monic cubic: x^3 + A*x^2 + B*x + C = 0
        A = a2 / a1
        B = a3 / a1
        C = a4 / a1

        # Cardano's method (matching Fortran complex arithmetic)
        p = B - A^2 / 3.0
        q = 2.0 / 27.0 * A^3 - A * B / 3.0 + C

        D1 = q^2 / 4.0 + p^3 / 27.0

        # Three cubic roots via complex cube roots
        if D1 >= 0.0
            sqrt_D1 = sqrt(D1)
            u_real = -q * 0.5 + sqrt_D1
            v_real = -q * 0.5 - sqrt_D1
            u = cbrt(u_real)
            v = -p / (3.0 * u)

            X1 = u + v
            # Two complex roots (not needed, Fortran picks real root first)
            An1 = X1 - A / 3.0
            An2 = NaN  # complex root, skip
            An3 = NaN  # complex root, skip
        else
            # Three real roots
            sqrt_neg_D1 = sqrt(-D1)
            # u = complex(-q/2, sqrt_neg_D1), v = conj(u)
            u_r, u_i = -q * 0.5, sqrt_neg_D1
            r = sqrt(u_r^2 + u_i^2)
            θ = atan(u_i, u_r)

            cbrt_r = cbrt(r)
            X1 = 2.0 * cbrt_r * cos(θ / 3.0)
            X2 = 2.0 * cbrt_r * cos((θ + 2π) / 3.0)
            X3 = 2.0 * cbrt_r * cos((θ + 4π) / 3.0)

            An1 = X1 - A / 3.0
            An2 = X2 - A / 3.0
            An3 = X3 - A / 3.0
        end

        # Quadratic solutions (An < 0 branch)
        b1 = ee[i] * G0_CO2 + ee[i] * GB_CO2
        b2 = -G0_CO2 * ψ - GB_CO2 * ζ
        b3 = G0_CO2 * ω

        D2 = b2^2 - 4.0 * b1 * b3
        An4 = D2 >= 0.0 ? (-b2 + sqrt(D2)) / (2.0 * b1) : NaN
        An5 = D2 >= 0.0 ? (-b2 - sqrt(D2)) / (2.0 * b1) : NaN
        An6 = 0.0

        # Select valid solution (following Fortran logic exactly)
        solutions = [(An1, true), (An2, true), (An3, true), (An4, false), (An5, false), (An6, false)]
        found = false

        for (An_val, require_positive) in solutions
            isnan(An_val) && continue

            Cs = CO2_atmosphere - An_val / GB_CO2

            if require_positive
                # J=1..3: use Ball-Berry for stomatal conductance
                Cs == 0.0 && continue
                Gsc = An_val * G1_CO2 * Rh / Cs + G0_CO2
            else
                # J=4..6: stomata at minimum conductance
                Gsc = G0_CO2
            end

            Ci = Cs - An_val / Gsc

            if require_positive
                if An_val > 0.0 && Gsc > 0.0 && Ci > 0.0
                    if i == 1
                        limit_Rubisco_Ac = An_val + RSP_leaf
                    else
                        limit_RuBP_Aj = An_val + RSP_leaf
                    end
                    found = true
                    break
                end
            else
                if An_val < 0.0 && Gsc > 0.0 && Ci > 0.0
                    if i == 1
                        limit_Rubisco_Ac = An_val + RSP_leaf
                    else
                        limit_RuBP_Aj = An_val + RSP_leaf
                    end
                    found = true
                    break
                end
            end
        end

        if !found
            error("leaf_photosynthesis_c3: no valid solution for limitation case $i")
        end
    end

    # ===== 10. TPU-limited rate =====
    limit_TPU_Ap = Vm * 0.5  # [μmol/m2/s]

    # ===== 11. Co-limitation (Farquhar et al.) =====
    # Rubisco + RuBP co-limitation
    square_root = max((limit_RuBP_Aj + limit_Rubisco_Ac)^2 - 4.0 * atheta * limit_RuBP_Aj * limit_Rubisco_Ac, 0.0)
    limit_both_Ai = ((limit_RuBP_Aj + limit_Rubisco_Ac) - sqrt(square_root)) / (2.0 * atheta)

    # Final GPP with TPU limitation
    square_root = max((limit_both_Ai + limit_TPU_Ap)^2 - 4.0 * btheta * limit_both_Ai * limit_TPU_Ap, 0.0)
    GPP_leaf = ((limit_TPU_Ap + limit_both_Ai) - sqrt(square_root)) / (2.0 * btheta)  # [μmol/m2/s]

    # ===== 12. Net assimilation =====
    NetAssim_leaf = GPP_leaf - RSP_leaf  # [μmol/m2/s]

    # ===== 13. Stomatal conductance =====
    Cs = CO2_atmosphere - NetAssim_leaf / GB_CO2  # [μmol/mol]

    if NetAssim_leaf > 0.0
        GS_CO2 = G1_CO2 * Rh * NetAssim_leaf / Cs + G0_CO2  # [mol/m2/s for CO2]
    else
        GS_CO2 = G0_CO2  # [mol/m2/s for CO2]
    end
    GS_H2O = GS_CO2 * 1.6  # [mol/m2/s for H2O]

    # ===== 14. Transpiration =====
    H2O_surface = min(H2O_internal, (GB_H2O * H2O_ambient + GS_H2O * H2O_internal) / (GB_H2O + GS_H2O))
    TSP_leaf = GS_H2O * (H2O_internal - H2O_surface) * M_H2O  # [kg/m2/s]

    # Convert to mol/m2/s (matching Fortran GPPLF_D = GPPLF/1000000)
    return (gpp=GPP_leaf / 1e6, rsp=RSP_leaf / 1e6, tsp=TSP_leaf)
end