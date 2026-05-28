# SOIL - 5-layer Richards equation soil water model
# Paper: Masutomi et al. (2016) https://doi.org/10.5194/gmd-9-4133-2016
# Soil water: Richards equation with Campbell (1985) / Clapp-Hornberger (1978)
# Evaporation: Sellers (1996), FAO56 aerodynamic resistance


# ============ soil-specific constants ============
const n_layer = 5                                      # number of soil layers
const layer_width = [0.05, 0.20, 0.75, 1.00, 2.00]       # layer thickness [m]
const ZA  = 2.0                                    # reference height [m]

# Soil texture parameters (13 USDA textures, Table 6)
#  1:heavy clay  2:silty clay  3:clay  4:silty clay loam  5:clay loam
#  6:silt  7:silt loam  8:sandy clay  9:loam  10:sandy clay loam
#  11:sandy loam  12:loamy sand  13:sand
const ψ_sat_T = [-3.7, -3.4, -3.7, -3.3, -2.6, -2.1, -2.1, -2.9, -1.1, -2.8, -1.5, -0.9, -0.7]   # [×1e6 Pa]
const B_frac_T = [7.6, 7.9, 7.6, 6.6, 5.2, 4.7, 4.7, 6.0, 4.5, 4.0, 3.1, 2.1, 1.7]                 # B (Clapp-Hornberger)
const K_sat_T = [1.7e-5, 2.5e-5, 1.7e-5, 4.2e-5, 6.4e-5, 1.9e-4, 1.9e-4, 3.3e-5, 3.7e-4, 1.2e-4, 7.2e-4, 1.7e-3, 5.8e-3]  # K_s [m/s]
# const BLKS_T = [1330., 1260., 1330., 1300., 1390., 1380., 1380., 1470., 1430., 1500., 1460., 1430., 1430.]  # ρ_s [kg/m³]
const porosity_T = [0.50, 0.52, 0.50, 0.51, 0.48, 0.48, 0.48, 0.44, 0.46, 0.43, 0.45, 0.46, 0.46]   # w_sat [-]
const wilting_point_T = [0.30, 0.27, 0.30, 0.22, 0.22, 0.06, 0.11, 0.25, 0.14, 0.17, 0.08, 0.05, 0.05]   # w_wlt [-]
const field_capacity_T   = [0.42, 0.41, 0.42, 0.38, 0.36, 0.30, 0.31, 0.36, 0.28, 0.27, 0.18, 0.12, 0.10]   # field capacity [-]

# ============ main function ============
function calc_soil_water(;
    layer_water::Vector{Float64},    # soil water content per layer [m³/m³]
    transpiration::Float64,          # transpiration E_t [W/m²]
    W2SF::Float64,                   # water to soil surface F_c [kg/m²/s]
    depth_root::Float64,             # rooting depth z_rt [m]
    is_irrigated::Int,               # irrigation flag (0=rainfed, 1=flooded paddy)
    Δt::Int,                         # time step Δt [s]
    soil_type_i::Int,                # soil texture index (1-13)
    temperature::Float64,            # air temperature T [K]
    pressure::Float64,               # air pressure P_a [Pa]
    wind_speed::Float64,             # wind speed at 2m [m/s]
    specific_humidity::Float64,      # specific humidity Q [kg/kg]
    crop_height::Float64,            # crop height [m]
    is_planted::Int,                 # plant flag (0=no plant, 1/2=plant)
    crop_name::String)               # crop name

    # ===== 1. Soil texture parameters (Table 6) =====
    ψ_sat = Float64(ψ_sat_T[soil_type_i])  # ψ_s: saturated water potential (code units, Table 6)
    B_frac = Float64(B_frac_T[soil_type_i])  # B: a parameter that determines the relationship of hydraulic conductivity or water potentials between saturated and unsaturated soils
    K_sat = K_sat_T[soil_type_i]           # K_sat: saturated hydraulic conductivity [m/s]
    porosity = porosity_T[soil_type_i]           # w_sat: porosity [-]
    wilting_point = wilting_point_T[soil_type_i]           # w_wlt: wilting point [-]
    field_capacity = field_capacity_T[soil_type_i]             # field capacity [-]

    # θ = θ_saturation * (ψ_air_dryed / ψ_sat)^(-1.0 / B_frac)
    # absolute_potential_total_soil_ =-3e7 [Pa] (measured value, P107 in Saishindojyogaku)
    # ψ_air_dryed = absolute_potential_total_soil_ / ρ_water
    # θ_saturation = porosity
    # ADSW: Air-dried soil water
    # when soil water reach such threshold, no more evaporation
    ADSW = porosity * ((-3e7 / ρ_water) / ψ_sat)^(-1.0 / B_frac)


    # ===== 2. Layer width / depth / soil water content =====
    layer_width_ave = [(layer_width[i] + layer_width[i+1]) * 0.5 for i in 1:n_layer-1]
    depth = zeros(Float64, n_layer + 1)   # depth[1]=0, depth[i+1]=depth to bottom of layer i
    for i in 1:n_layer
        depth[i+1] = i == 1 ? layer_width[1] : depth[i] + layer_width[i]
    end

    # Clamp soil water
    for i in 1:n_layer
        layer_water[i] = min(layer_water[i], porosity)
    end


    # ===== 3. Root distribution (Eq.59) =====
    lowest_root_layer_num = 0  # Index of the lowest layer where root exists
    f_r = zeros(Float64, n_layer)

    if depth_root > 0.0
        for i in 1:n_layer
            if depth_root < depth[i+1]
                lowest_root_layer_num = i
                break
            end
        end

        for i in 1:n_layer
            d_top = depth[i]
            d_bot = depth[i+1]
            if i < lowest_root_layer_num
                # Eq.59: integrate f_r(z) = ∫ (3/2)(z_rt-z²)/z_rt³ over [d_top, d_bot]
                # the above formula is the derivative version 
                # the cumulative version should be
                # 1.5 / z_rt^3 * (z_rt^2 * z - 1/3 * z^3)
                # f_r[i] should be the rood distribution in layer i
                f_r[i] = 1.5 / depth_root^3 * (depth_root^2 * d_bot - d_bot^3 / 3 - (depth_root^2 * d_top - d_top^3 / 3))
            elseif i == lowest_root_layer_num
                # Root tip in this layer: integrate from d_top to depth_root
                f_r[i] = 1.5 / depth_root^3 * (2.0 * depth_root^3 / 3 - (depth_root^2 * d_top - d_top^3 / 3))
            else
                f_r[i] = 0.0
            end
        end
    end

    # ===== 4. Transpiration extraction (Eq.58) =====
    # S_s(z) = (E_t / ρ_w) * f_r(z), applied layer-by-layer
    evp_deficit = 0.0  # evapotranspiration deficit [W/m²]

    if lowest_root_layer_num > 0
        for i in 1:lowest_root_layer_num
            extract = f_r[i] * transpiration / L_vaporization * Float64(Δt) / ρ_water / layer_width[i]
            if (layer_water[i] - wilting_point) > extract
                layer_water[i] -= extract
            else
                evp_deficit += f_r[i] * transpiration - (layer_water[i] - wilting_point) * L_vaporization / Float64(Δt) * ρ_water * layer_width[i]
                layer_water[i] = wilting_point
            end
        end

        # Redistribute deficit across layers
        if evp_deficit > 0.0
            for i in 1:lowest_root_layer_num
                extract_deficit = evp_deficit / L_vaporization * Float64(Δt) / ρ_water / layer_width[i]
                if (layer_water[i] - wilting_point) > extract_deficit
                    layer_water[i] -= extract_deficit
                    evp_deficit = 0.0
                    break
                else
                    evp_deficit -= (layer_water[i] - wilting_point) * L_vaporization / Float64(Δt) * ρ_water * layer_width[i]
                    layer_water[i] = wilting_point
                end
            end
        end
    end

    # ===== 5. Soil evaporation (calculation mainly on top soil) =====
    # Eq.57: ψ(z) = ψ_s * (w_s(z)/w_sat)^(-B)
    # ψ_z is the water potential at deth Z
    # here the arithmetic operations code with '.', such as './'
    # is the vectorized version of the formula
    W = layer_water ./ porosity
    ψ = ψ_sat .* W.^(.-B_frac)

    # Rsoil: soil surface resistance r_s
    # 800 is the reference value
    # when top soil layer is fully dry, W[1] ≈ 0, Rsoil ≈ 800
    # when top soil layer is fully wet, W[1] ≈ 1, Rsoil ≈ 0, no evaporation
    # W[2] is the second top soil,  for some small correction maybe
    Rsoil = 800.0 * (1.0 - W[1]) / (0.2 + W[1])              # [s/m]

    # Eq.60: topsoil humidity h_ms (Kelvin equation)
    soil_humidity = exp(ψ[1] * g0 / R_air / temperature)                      # [-]

    # Air density
    ρ_air = pressure / R_air / temperature                                   # [kg/m³]

    # Aerodynamic resistance (FAO56 P20)
    if crop_height > 0.0
        Zd  = 2.0 / 3.0 * crop_height
        Zom = 0.123 * crop_height
        Zoh = 0.1 * Zom
        Rair  = log((ZA - Zd) / Zom) * log((ZA - Zd) / Zoh) / (karman_constant^2 * wind_speed)  # [s/m]
    else
        # Bare soil: Liu et al. (2007) HESS
        Rair = 94.909 * wind_speed^(-0.9036)                             # [s/m]
    end

    EVS_max = layer_water[1] * layer_width[1] * ρ_water / Float64(Δt) * L_vaporization       # [W/m²]

    EVS = 0.0
    if soil_humidity * saturation_vapor_pressure(temperature, pressure) > specific_humidity
        EVS = (soil_humidity * saturation_vapor_pressure(temperature, pressure) - specific_humidity) * L_vaporization * ρ_air / (Rsoil + Rair)
        EVS = min(EVS, EVS_max)
    end

    # Subtract evaporation from top layer (not for flooded paddy with plant)
    if is_irrigated == 0 || is_planted == 0
        layer_water[1] -= EVS / L_vaporization * Float64(Δt) / ρ_water / layer_width[1]
    end

    for i in 1:n_layer
        layer_water[i] = clamp(layer_water[i], ADSW, porosity)
    end

    # ===== 6. Eq.56-57  =====
    # W is an intermediate variable for next two steps
    W = layer_water ./ porosity
    # K(z) = K_s * (w_s/w_sat)^(2B+3)
    K = K_sat .* W.^(2.0 * B_frac + 3.0)
    # ψ(z) = ψ_s * (w_s/w_sat)^(-B)
    ψ = ψ_sat .* W.^(.-B_frac)

    # Inter-layer conductivity (thickness-weighted average)
    KB = [(K[i] * layer_width[i] + K[i+1] * layer_width[i+1]) / (layer_width[i] + layer_width[i+1]) for i in 1:n_layer-1]

    # Darcy flux between layers
    # Q = K * A * Δh / ΔL
    flux_water = zeros(Float64, n_layer + 1)
    # Surface infiltration
    flux_water[1] = min(W2SF, (porosity - layer_water[1]) * ρ_water * layer_width[1] / Float64(Δt))

    for i in 1:n_layer-1
        flux_water[i+1] = -KB[i] * ((ψ[i+1] - ψ[i]) / layer_width_ave[i] - g0)
    end

    # Bottom: no flux (base runoff disabled in this version)
    flux_water[n_layer+1] = 0.0

    # Build tridiagonal system: A * Δw = b
    # Taylor Expand
    # flux(w+Δw) ≈ flux(w) + ∂flux/∂w_i * Δw_i + ∂flux/∂w_(i+1) *Δw_(i+1)
    # ∂flux/∂w_i * Δw_i + ∂flux/∂w_(i+1) *Δw_(i+1) ≈ flux(w+Δw) - flux(w)
    # Δw_i & Δw_(i+1) is the Δw
    # ∂flux/∂w_i & ∂flux/∂w_(i+1) is the A
    # flux(w+Δw) - flux(w) is the b
    Δflux_water = zeros(Float64, n_layer, 2)
    for i in 1:n_layer-1
        Δflux_water[i, 1] = (-KB[i] / layer_width_ave[i]) * (B_frac * ψ_sat / porosity * (layer_water[i] / porosity)^(-B_frac - 1.0))
        Δflux_water[i, 2] = (-KB[i] / layer_width_ave[i]) * (B_frac * ψ_sat / porosity * (layer_water[i+1] / porosity)^(-B_frac - 1.0))
    end

    A = zeros(Float64, n_layer, n_layer)
    b = zeros(Float64, n_layer)

    for i in 1:n_layer
        if i == 1
            A[1, 1] = -ρ_water * layer_width[1] / Float64(Δt) - Δflux_water[1, 1]
            A[1, 2] =                                         - Δflux_water[1, 2]
            b[1]    = -flux_water[1] + flux_water[2]
        elseif i < n_layer
            A[i, i-1] =   Δflux_water[i-1, 1]
            A[i, i]   = - ρ_water * layer_width[i] / Float64(Δt) + Δflux_water[i-1, 2] - Δflux_water[i, 1]
            A[i, i+1] = - Δflux_water[i, 2]
            b[i]      = - flux_water[i] + flux_water[i+1]
        else  # i == n_layer
            A[n_layer, n_layer-1] =   Δflux_water[n_layer-1, 1]
            A[n_layer, n_layer]   = - ρ_water * layer_width[n_layer] / Float64(Δt) + Δflux_water[n_layer-1, 2]
            b[n_layer]            = - flux_water[n_layer] + flux_water[n_layer+1]
        end
    end

    # Solve tridiagonal system
    Δw = A \ b

    # Update soil water
    for i in 1:n_layer
        layer_water[i] += Δw[i]
    end

    # ===== 7. Irrigation / saturation (Eq.54) =====
    # Flooded paddy: w_s(z) = w_sat for all layers
    if is_irrigated == 1 && is_planted > 0
        for i in 1:n_layer
            layer_water[i] = porosity
        end
    end

    for i in 1:n_layer
        if i < n_layer
            layer_water[i] = min(layer_water[i], porosity)
        else
            layer_water[i] = min(layer_water[i], field_capacity)     # bottom layer drains to field capacity
        end
        layer_water[i] = max(layer_water[i], ADSW)
    end

    # ===== 8. Water stress =====
    # maize paper eq.15
    FAW     = [min(max(layer_water[i] - wilting_point, 0.0) / (field_capacity - wilting_point), 1.0) for i in 1:n_layer]
    FAWRICE = [min(max(layer_water[i] - wilting_point, 0.0) / (porosity - wilting_point), 1.0) for i in 1:n_layer]

    WSTRS = 0.0
    for i in 1:n_layer
        if crop_name == "Rice"
            # Paddy: FAW based on porosity (flooded = no stress), threshold=0.8
            FSTRS = FAWRICE[i] > 0.8 ? 1.0 : FAWRICE[i] / 0.8
        elseif crop_name == "Soybeans"
            # FAO56-style, threshold=0.5
            FSTRS = FAW[i] > 0.5 ? 1.0 : FAW[i] / 0.5
        elseif crop_name == "Wheat" || crop_name == "Maize"
            # FAO56 p=0.55 → threshold=0.45
            FSTRS = FAW[i] > 0.45 ? 1.0 : FAW[i] / 0.45
        else
            error("Unknown crop: $crop_name. Use Rice, Wheat, Soybeans, or Maize")
        end
        # Eq.78: weight by root distribution
        WSTRS += FSTRS * f_r[i]
    end

    return (water_stress=WSTRS, evaporation=EVS)
end