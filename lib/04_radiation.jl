# RAD - Two-stream canopy radiation model
# Paper: Masutomi et al. (2016) https://doi.org/10.5194/gmd-9-4133-2016
# Radiation transfer: Goudriaan and van Laar (1994)


# ============ constants (use Julia built-in π) ============
const f_leaf_orientation = 0.5       # leaf orientation factor F [-] (spherical distribution, Sect.3.2)
const ALBS = 0.1       # soil albedo r_g [-] (Eq.A7-A8)
const DF   = 1.66      # diffusivity factor d_f = sec(53°) (Kondoh 1994, Sect.3.2)
const FLFE = 1.0       # effective leaf PARam (0.75 in original MATSIRO)

# ============ main function ============
function calc_radiation(;
    leaf_nitrogen::Float64,  # specific leaf nitrogen [g N/m²(leaf)]
    kn::Float64,             # nitrogen extinction coefficient [-]
    shortwave_radiation::Float64,  # downward shortwave radiation R_s^d(0) [W/m²]
    LAI::Float64,            # leaf area index L [m²(leaf)/m²]
    RLFv::Float64,           # leaf PAR reflectance r_1 [-]
    TLFv::Float64,           # leaf PAR transmissivity t_1 [-]
    RLFn::Float64,           # leaf NIR reflectance r_2 [-]
    TLFn::Float64,           # leaf NIR transmissivity t_2 [-]
    lat::Float64,            # latitude [degree]
    doy::Int,                # day of year D_oy
    hour::Float64,           # hour of day
    crop_name::String,       # crop name ("Rice","Wheat","Soybean","Maize")
    development_stage::Float64)  # development stage [-]

    # ===== 1. Solar geometry =====
    # Appendix B
    cos_θ = sin_solar_elevation(doy, hour, lat)  # cos(zenith) = sin(solar elevation)

    # ===== 2. Leaf optical properties =====
    TLF = [TLFv, TLFn]  # [PAR, NIR]
    RLF = [RLFv, RLFn]

    if LAI <= 0.0
        return (PAR_abs_sunlit_leaf=0.0, PAR_abs_shade_leaf=0.0,
                Vmax_sunlit_leaf=0.0, Vmax_shade_leaf=0.0,
                LAI_sunlit=0.0, LAI_shade=0.0)
    end

    RFE = zeros(Float64, 2)  # effective leaf reflectivity
    TFE = zeros(Float64, 2)  # effective leaf transmissivity
    SRF = zeros(Float64, 2)  #  part of the variable when calculating a_i in eq.C1
    K_diffuse = zeros(Float64, 2)  # a_i in paper, extinction coefficient for scattered radiation

    # eq.C1
    for i in 1:2
        RFE[i] = FLFE * RLF[i] + (1.0 - FLFE) * TLF[i]
        TFE[i] = FLFE * TLF[i] + (1.0 - FLFE) * RLF[i]
        SRF[i] = sqrt((1.0 - TFE[i])^2 - RFE[i]^2)
        K_diffuse[i] = f_leaf_orientation * DF * SRF[i]
    end

    # ===== 3. Canopy Vmax top (crop-specific leaf_nitrogen→Vmax) =====
    # Vmax_top: maximum Rubisco capacity at the canopy top
    # Masutomi et al. (2016)
    if crop_name == "Rice" || crop_name == "Wheat"
        Vmax_top = min(87.04 * (leaf_nitrogen - 0.487), 138.77)             # [μmol/m²/s]
    elseif crop_name == "Soybean"
        Vmax_top = min(max(-18.516 * leaf_nitrogen^2 + 114.33 * leaf_nitrogen - 73.336, 0.0), 103.0)
    elseif crop_name == "Maize"
        if development_stage < 0.52
            # Vegetative: Vos et al. (2005)
            Vmax_top = 45.1 * (2.0 / (1.0 + exp(-2.9 * (leaf_nitrogen - 0.25))) - 1.0)
        elseif development_stage < 1.01
            # Reproductive (after silking): Drouet and Bonhomme (2004)
            Vmax_top = 40.2 * (2.0 / (1.0 + exp(-1.41 * (leaf_nitrogen - 0.43))) - 1.0)
        else
            error("Maize development_stage=$development_stage out of range (must be < 1.01)")
        end
    else
        error("Unknown crop: $crop_name. Use Rice, Wheat, Soybean, or Maize")
    end

    # Canopy-average Vmax (CLM scheme)
    Vmax = (Vmax_top / 1e6) * (1.0 - exp(-kn * LAI)) / kn # [mol/m²/s]

    # ===== 4. PAR / NIR split =====

    # ? Eqs.15-16: PAR = NIR = 0.5*R_s^d(0) (assuming equal split)
    PAR = 0.5 * shortwave_radiation  # [W/m²]
    NIR = shortwave_radiation - PAR  # [W/m²]
    if NIR < 0.0
        NIR = 0.0
        shortwave_radiation = PAR
    end

    # ===== 5. Direct/diffuse PARtitioning (Goudriaan & van Laar 1994) =====
    # Eq.19: 
    # sc: extraterrestrial radiation R_ex
    radiation_extraterrestrial = 1370.0 * (1.0 + 0.033 * cos(2.0 * π * Float64(doy) / 365.0))  # [W/m²]

    if cos_θ > 0.0
        # ===== 5a. Sun above horizon =====
        sec_θ = 1.0 / cos_θ   # sec(theta)

        # Eq.18: atmospheric transmissivity tau_atm
        transmissivity_atmosphere = shortwave_radiation / (radiation_extraterrestrial * cos_θ)

        # Eq.17: fraction of scattered radiation to total radiation
        if transmissivity_atmosphere <= 0.22
            frac_rad = 1.0
        elseif transmissivity_atmosphere <= 0.35
            frac_rad = 1.0 - 6.4 * (transmissivity_atmosphere - 0.22)^2
        else
            frac_rad = 1.47 - 1.66 * transmissivity_atmosphere
        end
        frac_rad = max(frac_rad, 0.15 + 0.85 * (1.0 - exp(-0.1 / cos_θ)))
        
        # Eqs.15-16: 
        # S_i^d(0): the radiant flux density for downward scattered at top
        dif = [PAR * frac_rad, NIR * frac_rad] 
        # D_i^d(0): the radiant flux density for downward direct at top
        dir = [PAR - dif[1], NIR - dif[2]]  

        # ===== 6. Two-stream radiation transfer coefficients =====
        eap = zeros(Float64, 2)
        ean = zeros(Float64, 2)
        coeff_A1 = zeros(Float64, 2)  # A_{1,i}
        coeff_A2 = zeros(Float64, 2)  # A_{2,i}
        coeff_A3 = zeros(Float64, 2)  # A_{3,i}
        coeff_C1 = zeros(Float64, 2)  # C_{1,i}
        coeff_C2 = zeros(Float64, 2)  # C_{2,i}
        coeff_C3 = zeros(Float64, 2)  # C_{3,i}
        coeff_C4 = zeros(Float64, 2)  # C_{4,i}

        for i in 1:2
            # Eq.9 context: exp(±a_i*L)
            eap[i] = exp(K_diffuse[i] * LAI)
            ean[i] = exp(-K_diffuse[i] * LAI)

            # Eq.12: direct beam attenuated through canopy
            rad_direct_extinct = dir[i] * exp(-f_leaf_orientation * LAI * sec_θ)       # D_i^d(L)

            # Eq.C6: A_{1,i} = (1-t_i + sqrt(...)) / r_i
            coeff_A1[i] = (1.0 - TFE[i] + SRF[i]) / RFE[i]
            # Eq.C7: A_{2,i} = (1-t_i - sqrt(...)) / r_i
            coeff_A2[i] = (1.0 - TFE[i] - SRF[i]) / RFE[i]
            # Eq.C8 A_{3,i} = (A1-r_g)*exp(aL) - (A2-r_g)*exp(-aL)
            coeff_A3[i] = (coeff_A1[i] - ALBS) * eap[i] - (coeff_A2[i] - ALBS) * ean[i]

            # Eq.C4: C_{3,i}
            denom = DF^2 * ((1.0 - TFE[i])^2 - RFE[i]^2) - sec_θ^2
            coeff_C3[i] = sec_θ * (TFE[i] * sec_θ + DF * TFE[i] * (1.0 - TFE[i]) + DF * RFE[i]^2) / denom
            # Eq.C5: C_{4,i}
            coeff_C4[i] = RFE[i] * (DF - sec_θ) * sec_θ / denom

            # Eq.C2: C_{1,i}
            coeff_C1[i] = (-(coeff_A2[i] - ALBS) * ean[i] * (dif[i] - coeff_C3[i] * dir[i]) +
                           (coeff_C3[i] * ALBS + ALBS - coeff_C4[i]) * rad_direct_extinct) / coeff_A3[i]
            # Eq.C3: C_{2,i}
            coeff_C2[i] = ((coeff_A1[i] - ALBS) * eap[i] * (dif[i] - coeff_C3[i] * dir[i]) -
                           (coeff_C3[i] * ALBS + ALBS - coeff_C4[i]) * rad_direct_extinct) / coeff_A3[i]
        end

        # ===== 7. Sunlit / shaded LAI =====
        # ? Eq.65 context: LAI_sunlit = integral of exp(-F*sec(theta)*l) dl
        LAI_sunlit = (1.0 - exp(-f_leaf_orientation * sec_θ * LAI)) / (f_leaf_orientation * sec_θ)
        LAI_shade  = LAI - LAI_sunlit

        # ===== 8. Absorbed PAR by canopy components =====
        # Direct PAR absorbed by sunlit leaves (from Eq.12)
        PAR_abs_sunlit_direct = dir[1] * (1.0 - exp(-f_leaf_orientation * sec_θ * LAI))  # [W/m²]

        # Diffuse PAR absorbed by entire canopy (from Eq.13 integral)
        PAR_abs_diffuse_total = coeff_C1[1] * (1.0 - coeff_A1[1]) * (1.0 - exp(K_diffuse[1] * LAI)) +
                                coeff_C2[1] * (1.0 - coeff_A2[1]) * (1.0 - exp(-K_diffuse[1] * LAI)) +
                                (coeff_C3[1] - coeff_C4[1]) * dir[1] * (1.0 - exp(-f_leaf_orientation * sec_θ * LAI))
        PAR_abs_diffuse_total = max(PAR_abs_diffuse_total, 0.0)

        # Diffuse PAR absorbed by sunlit leaves (from Eq.13 weighted by exp(-F*sec*l))
        k_s = f_leaf_orientation * sec_θ  # direct beam extinction coefficient
        PAR_abs_sunlit_diffuse = K_diffuse[1] * coeff_C1[1] * (1.0 - coeff_A1[1]) /
                                 (K_diffuse[1] - k_s) * (1.0 - exp((K_diffuse[1] - k_s) * LAI)) +
                                 K_diffuse[1] * coeff_C2[1] * (1.0 - coeff_A2[1]) /
                                 (K_diffuse[1] + k_s) * (1.0 - exp(-(K_diffuse[1] + k_s) * LAI)) +
                                 (coeff_C3[1] - coeff_C4[1]) * dir[1] / 2.0 * (1.0 - exp(-2.0 * k_s * LAI))
        PAR_abs_sunlit_diffuse = max(PAR_abs_sunlit_diffuse, 0.0)

        # Diffuse PAR absorbed by shaded leaves
        PAR_abs_shade_diffuse = max(PAR_abs_diffuse_total - PAR_abs_sunlit_diffuse, 0.0)

        # Total absorbed PAR per leaf area
        PAR_abs_sunlit_total = PAR_abs_sunlit_direct + PAR_abs_sunlit_diffuse
        PAR_abs_shade_total = PAR_abs_shade_diffuse

        PAR_abs_sunlit_leaf = LAI_sunlit > 0.0 ? PAR_abs_sunlit_total / LAI_sunlit : 0.0
        PAR_abs_shade_leaf  = LAI_shade  > 0.0 ? PAR_abs_shade_total / LAI_shade  : 0.0

        # ===== 9. Vmax PARtitioning (sunlit / shade) =====
        # CLM scheme: Vmax_sunlit integrates over sunlit LAI with (kn+F*sec(theta)) extinction
        Vmax_sunlit = (Vmax_top / 1e6) * (1.0 - exp(-(kn + f_leaf_orientation * sec_θ) * LAI)) / (kn + f_leaf_orientation * sec_θ)
        Vmax_sunlit = min(Vmax, Vmax_sunlit)
        Vmax_shade  = Vmax - Vmax_sunlit

        Vmax_sunlit_leaf = LAI_sunlit > 0.0 ? Vmax_sunlit / LAI_sunlit : 0.0
        Vmax_shade_leaf  = LAI_shade  > 0.0 ? Vmax_shade  / LAI_shade  : 0.0

    else
        # ===== 5b. Diffuse only (sun below horizon) =====
        frac_rad = 1.0
        dif = [PAR * frac_rad, NIR * frac_rad]

        eap = zeros(Float64, 2)
        ean = zeros(Float64, 2)
        coeff_A1 = zeros(Float64, 2)
        coeff_A2 = zeros(Float64, 2)
        coeff_A3 = zeros(Float64, 2)
        coeff_C1 = zeros(Float64, 2)
        coeff_C2 = zeros(Float64, 2)

        for i in 1:2
            # temporary variables used in following formulas
            eap[i] = exp(K_diffuse[i] * LAI)
            ean[i] = exp(-K_diffuse[i] * LAI)
            # eq. C6-C8
            coeff_A1[i] = (1.0 - TFE[i] + SRF[i]) / RFE[i]
            coeff_A2[i] = (1.0 - TFE[i] - SRF[i]) / RFE[i]
            coeff_A3[i] = (coeff_A1[i] - ALBS) * eap[i] - (coeff_A2[i] - ALBS) * ean[i]
            # eq. C2-C3
            coeff_C1[i] = (-(coeff_A2[i] - ALBS) * ean[i] * dif[i]) / coeff_A3[i]
            coeff_C2[i] = ((coeff_A1[i] - ALBS) * eap[i] * dif[i]) / coeff_A3[i]
        end

        # Diffuse PAR only (no sunlit leaves)
        PAR_abs_diffuse_total = coeff_C1[1] * (1.0 - coeff_A1[1]) * (1.0 - exp(K_diffuse[1] * LAI)) +
                                coeff_C2[1] * (1.0 - coeff_A2[1]) * (1.0 - exp(-K_diffuse[1] * LAI))
        PAR_abs_diffuse_total = max(PAR_abs_diffuse_total, 0.0)

        LAI_sunlit = 0.0
        LAI_shade  = LAI

        PAR_abs_sunlit_leaf = 0.0
        PAR_abs_shade_leaf  = LAI_shade > 0.0 ? PAR_abs_diffuse_total / LAI_shade : 0.0
        Vmax_sunlit_leaf = 0.0
        Vmax_shade_leaf  = LAI_shade > 0.0 ? Vmax / LAI_shade : 0.0
    end

    return (PAR_abs_sunlit_leaf=PAR_abs_sunlit_leaf,
            PAR_abs_shade_leaf=PAR_abs_shade_leaf,
            Vmax_sunlit_leaf=Vmax_sunlit_leaf,
            Vmax_shade_leaf=Vmax_shade_leaf,
            LAI_sunlit=LAI_sunlit,
            LAI_shade=LAI_shade)
end
