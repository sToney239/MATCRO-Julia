# CROP - Crop growth module
# Paper: Masutomi et al. (2016) https://doi.org/10.5194/gmd-9-4133-2016
# Merged from maize/SUB_CROP.f90 and soy/SUB_CROP.f90
# Crop-specific branches: Maize PRTSHT, CALSLN Soybean/Maize, CALLAI CO2 down-regulation


# CFSR2GL: starch to glucose conversion factor
const fraction_starch_convert = 1.11   

# ============ mutable crop state ============
@kwdef mutable struct CropState
    is_planted::Int = 0
    has_emerged::Bool = false
    is_grain_filling::Int = 0

    development_stage::Float64 = 0.0
    progress_at_grain_start::Float64 = 0.0

    accumulated_thermal_time::Float64 = 0.0
    accumulated_vern_days::Float64 = 0.0 # vernalization days accumulated

    leaf_biomass::Float64 = 0.0
    stem_biomass::Float64 = 0.0
    storage_organ_biomass::Float64 = 0.0
    root_biomass::Float64 = 0.0
    reserved_starch_pool::Float64 = 0.0       # WSR: Biomass shielded reserve, could not be used
    available_glucose_pool::Float64 = 0.0     # WAR: immediately available glucose pool
    dead_leaf_biomass::Float64 = 0.0
    grain_biomass::Float64 = 0.0
    potential_grain_biomass::Float64 = 0.0
    shoot_biomass::Float64 = 0.0

    LAI::Float64 = 0.0
    LAI_max::Float64 = 0.0
    LAI_max_season::Float64 = 0.0
    shoot_biomass_at_harvest::Float64 = 0.0
    crop_height::Float64 = 0.0
    root_length::Float64 = 0.0

    cold_damage_index::Float64 = 0.0
    heat_damage_index::Float64 = 0.0
    n_heat_events::Float64 = 0.0
    daily_avg_temp::Float64 = 0.0
    daily_max_temp::Float64 = 0.0

    leaf_nitrogen::Float64 = 0.0

    harvest_doy::Int = 0
    yield::Float64 = 0.0
end

# ============ piecewise linear interpolation ============
linear_interpolate(x, x1, y1, x2, y2) = (y2 - y1) / (x2 - x1) * (x - x2) + y2

# ============ 1. Judge planting (Eq.110 context) ============
function judge_planting!(crop::CropState, doy::Int, hour::Float64,
                         planting_doy_input::Int, is_irrigated::Int)
    # Normalize planting DOY to 1-365
    planting_doy = planting_doy_input
    if planting_doy < 1
        planting_doy += 365
    end
    if planting_doy > 365
        planting_doy -= 365
    end

    if doy == planting_doy && crop.is_planted == 0 && hour >= 12.0
        crop.is_planted = 1
        # need to refresh to 0 if multiple year of cropping
        crop.has_emerged = false
        crop.is_grain_filling = 0

        crop.leaf_biomass = 0.0
        crop.stem_biomass = 0.0
        crop.root_biomass = 0.0
        crop.storage_organ_biomass = 0.0
        crop.reserved_starch_pool = 0.0
        crop.available_glucose_pool = 0.0
        crop.dead_leaf_biomass = 0.0
        crop.grain_biomass = 0.0
        crop.potential_grain_biomass = 0.0

        crop.shoot_biomass = 0.0

        crop.accumulated_thermal_time = 0.0
        crop.accumulated_vern_days = 0.0
        crop.development_stage = 0.0
        crop.progress_at_grain_start = 0.0
        crop.root_length = 0.0

        crop.cold_damage_index = 0.0
        crop.heat_damage_index = 0.0
        crop.n_heat_events = 0.0

        crop.LAI = 0.0

        crop.LAI_max_season = 0.0
        crop.shoot_biomass_at_harvest = 0.0
    end
end

# ============ 2. Development stage (Eq.110-112) ============
# Eq.110: DVS = accumulated_thermal_time / thermal_time_requirement
# Eq.111: accumulated_thermal_time = integral of DVR dt
# Eq.112: DVR = piecewise linear function of temperature
function calc_development_stage!(crop::CropState, temperature::Float64, Δt::Int,
                                thermal_time_requirement::Float64, needs_vernalization::Int,
                                base_temp::Float64, optimal_temp::Float64,
                                ceiling_temp::Float64, vernalization_saturation::Float64)
    # Eq.112: Development rate (DVR)
    T_celsius = temperature - T_ice  # K → °C

    if T_celsius < base_temp
        dvr = 0.0
    elseif T_celsius < optimal_temp
        dvr = T_celsius - base_temp
    elseif T_celsius < ceiling_temp
        dvr = (base_temp - optimal_temp) / (ceiling_temp - optimal_temp) * (T_celsius - ceiling_temp)
    else
        dvr = 0.0
    end

    # Vernalization temperature parameters (Eq.112)
    TV1 = -4.0   # °C, lower bound
    TV2 = 3.0    # °C, range start
    TV3 = 10.0   # °C, range end
    TV4 = 17.0   # °C, upper bound
    # Vernalization factor (VF)
    vern_factor = 1.0
    if needs_vernalization > 0
        # Vernalization efficiency
        if T_celsius < TV1
            vern_efficiency = 0.0
        elseif T_celsius < TV2
            vern_efficiency = (T_celsius - TV1) / (TV2 - TV1)
        elseif T_celsius < TV3
            vern_efficiency = 1.0
        else
            vern_efficiency = max(0.0, (TV4 - T_celsius) / (TV4 - TV3))
        end

        crop.accumulated_vern_days += vern_efficiency * (Float64(Δt) / 3600.0) / 24.0

        vb = vernalization_saturation / 5.0
        if crop.accumulated_vern_days < vb
            vern_factor = 0.0
        elseif crop.accumulated_vern_days < vernalization_saturation
            vern_factor = (crop.accumulated_vern_days - vb) / (vernalization_saturation - vb)
        else
            vern_factor = 1.0
        end
    end

    # Eq.111: Accumulate thermal time
    crop.accumulated_thermal_time += dvr * (Float64(Δt) / 3600.0) / 24.0 * vern_factor

    # Eq.110: DVS
    crop.development_stage = crop.accumulated_thermal_time / thermal_time_requirement
end

# ============ 3. Judge emergence ============
function judge_emergence!(crop::CropState)
    if crop.development_stage > 0.012 && !crop.has_emerged
        crop.leaf_biomass = 1.0
        crop.stem_biomass = 1.0
        crop.root_biomass = 1.0
        crop.available_glucose_pool = 0.5
        crop.has_emerged = true
    end
end

# ============ 4. Loss rates (Eq.133-135) ============
# Eq.133: Leaf loss rate
# Eq.134: Starch remobilization rate
# Eq.135: Leaf death ratio (DVS-dependent)
function calc_loss_rates(crop::CropState, Δt::Int, half_progress::Float64,
                         dead_prgress_1::Float64, dead_ratio_1::Float64,
                         dead_prgress_2::Float64, dead_ratio_2::Float64,
                         dead_prgress_3::Float64, dead_ratio_3::Float64)
    growth_progress = crop.development_stage

    # Eq.135: Leaf death ratio (DVS-dependent piecewise linear)
    if growth_progress < dead_prgress_1
        leaf_loss_ratio = dead_ratio_1
    elseif growth_progress < dead_prgress_2
        leaf_loss_ratio = linear_interpolate(growth_progress, dead_prgress_1, dead_ratio_1, dead_prgress_2, dead_ratio_2)
    elseif growth_progress < dead_prgress_3
        leaf_loss_ratio = linear_interpolate(growth_progress, dead_prgress_2, dead_ratio_2, dead_prgress_3, dead_ratio_3)
    else
        leaf_loss_ratio = dead_ratio_3
    end

    # Eq.133: Leaf & root Loss
    leaf_loss = (crop.leaf_biomass + crop.available_glucose_pool) * leaf_loss_ratio * Float64(Δt)
    root_loss_ratio = 0.0
    root_loss = crop.root_biomass * root_loss_ratio * Float64(Δt)
    
    # Time for Starch remobilization (10 days after heading)
    # Preparing ration of leaf remobilization for Eq.134
    if growth_progress < half_progress
        starch_mobilization_time = 0.0
    else
        starch_mobilization_time = 10.0   # Bouman (2001) [days]
    end
    # Eq.134: Starch remobilization
    if starch_mobilization_time > 0.0
        starch_loss = crop.reserved_starch_pool * (1.0 / (starch_mobilization_time * seconds_per_day)) * Float64(Δt)
    else
        starch_loss = 0.0
    end

    return (leaf_loss=leaf_loss, root_loss=root_loss, starch_loss=starch_loss)
end

# ============ 5. Growth and partitioning (Eq.119-132) ============
# Eq.119: Glucose supply
# Eq.120-121: Critical threshold and available carbohydrate
# Eq.130: Shoot partition ratio (Pr,sh)
# Eq.131: Leaf partition ratio (Pr,lef)
# Eq.132: Panicle partition ratio (Pr,pnc)
function calc_growth_partitioning!(crop::CropState, net_assimilation::Float64,
                                  starch_loss::Float64, Δt::Int,
                                  half_progress::Float64,
                                  k_leaf_convert::Float64, k_stem_convert::Float64,
                                  k_root_convert::Float64, k_grain_convert::Float64,
                                  fraction_starch_reserve::Float64,
                                  shoot_progress_1::Float64,
                                  shoot_alloc_ratio_1::Float64, shoot_progress_2::Float64,
                                  leaf_alloc_ratio_0::Float64, 
                                  leaf_progress_1::Float64, leaf_alloc_ratio_1::Float64,
                                  leaf_progress_2::Float64, leaf_alloc_ratio_2::Float64,
                                  panicle_progress_1::Float64, panicle_alloc_ratio_1::Float64,
                                  panicle_progress_2::Float64, panicle_alloc_ratio_2::Float64,
                                  panicle_progress_3::Float64, panicle_alloc_ratio_3::Float64,
                                  LAI_threshold_grain::Float64,
                                  k_leaf_loss::Float64,
                                  crop_name::String)
    growth_progress = crop.development_stage
    glucose_pool = crop.available_glucose_pool
    glucose_pool_initial = glucose_pool

    # Eq.119: glucose_pool(Glucose supply) = Net assimilation + remobilization from starch
    remobilized_glucose = starch_loss * fraction_starch_convert

    # Net assimilation → glucose: [mol(CO2)/m²/s] → [kg(CH2O)/ha/Δt]
    assimilation_glucose = net_assimilation * Float64(Δt) * 300.0

    glucose_pool = glucose_pool + assimilation_glucose + remobilized_glucose
    glucose_pool = max(glucose_pool, 0.0)

    # Eq.120-121: Critical threshold for growth
    if glucose_pool > 0.1 * crop.leaf_biomass
        carbon_total_growth = (glucose_pool - 0.1 * crop.leaf_biomass) / (seconds_per_day * 0.5 / Float64(Δt))
    else
        carbon_total_growth = 0.0
    end
    carbon_total_growth = max(0.0, carbon_total_growth)
    glucose_pool -= carbon_total_growth

    glucose_consumed = glucose_pool_initial - glucose_pool

    #############################################################
    # Eq.130: Shoot partition ratio (Pr,sh)
    # The proportion between total and shoot carbonhydrate
    if crop_name == "Maize"
        if growth_progress < shoot_progress_1
            ratio_shoot_alloc = linear_interpolate(growth_progress, 0.0, 0.5, shoot_progress_1, 1.0 - shoot_alloc_ratio_1)
        elseif growth_progress < shoot_progress_2
            ratio_shoot_alloc = linear_interpolate(growth_progress, shoot_progress_1, 1.0 - shoot_alloc_ratio_1, shoot_progress_2, 1.0)
        else
            ratio_shoot_alloc = 1.0
        end
    elseif crop_name == "Rice"
        progress_transplant_start = 0.06     # DVS at transplanting and  (Dvs,tr)
        progress_transplant_end = 0.08       # DVS at transplanting shock ends (Dvs,te)
        if growth_progress < progress_transplant_start
            ratio_shoot_alloc = 1.0 - shoot_alloc_ratio_1
        elseif growth_progress < progress_transplant_end
            ratio_shoot_alloc = 0
        elseif growth_progress < shoot_progress_2
            ratio_shoot_alloc = 1.0 - shoot_alloc_ratio_1
        elseif growth_progress < shoot_progress_2
            ratio_shoot_alloc = linear_interpolate(growth_progress, shoot_progress_1, 1.0 - shoot_alloc_ratio_1, shoot_progress_2, 1.0)
        else
            ratio_shoot_alloc = 1.0
        end
    else
        if growth_progress < shoot_progress_1
            ratio_shoot_alloc = 1.0 - shoot_alloc_ratio_1
        elseif growth_progress < shoot_progress_2
            ratio_shoot_alloc = linear_interpolate(growth_progress, shoot_progress_1, 1.0 - shoot_alloc_ratio_1, shoot_progress_2, 1.0)
        else
            ratio_shoot_alloc = 1.0
        end
    end
 
    # Grain filling initiation
    if crop.LAI > LAI_threshold_grain && leaf_progress_1 < growth_progress && crop.is_grain_filling == 0
        crop.progress_at_grain_start = growth_progress
        crop.is_grain_filling = 1
    end

    #############################################################
    # Eq.131: Leaf partition ratio (Pr,lef)
    # Between shoot and leaf, Fig 2a in Maize and Soy paper
    if growth_progress < leaf_progress_1
        ratio_leaf_alloc = linear_interpolate(growth_progress, 0.0, leaf_alloc_ratio_0, leaf_progress_1, leaf_alloc_ratio_1)
    elseif growth_progress < leaf_progress_2
        if crop.LAI < LAI_threshold_grain
            ratio_leaf_alloc = linear_interpolate(growth_progress, leaf_progress_1, leaf_alloc_ratio_1, leaf_progress_2, leaf_alloc_ratio_2)
        else
            # Adjusted leaf partition breakpoints after grain start
            adjusted_leaf_progress_1 = crop.progress_at_grain_start
            adjusted_leaf_progress_2 = (leaf_progress_2 - adjusted_leaf_progress_1) * k_leaf_loss + adjusted_leaf_progress_1
            adjusted_leaf_alloc_ratio_1 = linear_interpolate(adjusted_leaf_progress_1, leaf_progress_1, leaf_alloc_ratio_1, leaf_progress_2, leaf_alloc_ratio_2)

            ratio_leaf_alloc = linear_interpolate(growth_progress, adjusted_leaf_progress_1, adjusted_leaf_alloc_ratio_1, adjusted_leaf_progress_2, leaf_alloc_ratio_2)
        end
    else
        ratio_leaf_alloc = 0.0
    end
    ratio_leaf_alloc = max(0.0, ratio_leaf_alloc)

    #############################################################
    # Eq.132: Panicle partition ratio (Pr,pnc)
    # Between shoot and Panicle, Fig 2b in Maize and Soy paper
    if growth_progress < panicle_progress_1
        ratio_panicle_alloc = panicle_alloc_ratio_1
    elseif growth_progress < panicle_progress_2
        ratio_panicle_alloc = linear_interpolate(growth_progress, panicle_progress_1, panicle_alloc_ratio_1, panicle_progress_2, panicle_alloc_ratio_2)
    elseif growth_progress < panicle_progress_3
        ratio_panicle_alloc = linear_interpolate(growth_progress, panicle_progress_2, panicle_alloc_ratio_2, panicle_progress_3, panicle_alloc_ratio_3)
    else
        ratio_panicle_alloc = 0.0
    end

    #############################################################
    
    # Normalize if leaf + panicle > 1
    if (ratio_leaf_alloc + ratio_panicle_alloc) > 1.0
        ratio_leaf_alloc    = ratio_leaf_alloc    / (ratio_leaf_alloc + ratio_panicle_alloc)
        ratio_panicle_alloc = ratio_panicle_alloc / (ratio_leaf_alloc + ratio_panicle_alloc)
    end

    # Eq.122-127: Allocate carbohydrate
    carbon_shoot_growh    = carbon_total_growth * ratio_shoot_alloc         
    carbon_root_growth    = carbon_total_growth - carbon_shoot_growh        
    carbon_leaf_growth    = carbon_shoot_growh * ratio_leaf_alloc        
    carbon_panicle_growth = carbon_shoot_growh * ratio_panicle_alloc    
    carbon_stem_growth    = carbon_shoot_growh - carbon_leaf_growth - carbon_panicle_growth 

    leaf_growth           = carbon_leaf_growth    * k_leaf_convert
    stem_growth           = carbon_stem_growth    * (1.0 - fraction_starch_reserve) * k_stem_convert
    starch_growth         = carbon_stem_growth    * fraction_starch_reserve / fraction_starch_convert
    storage_organ_growth  = carbon_panicle_growth * k_grain_convert
    root_growth           = carbon_root_growth    * k_root_convert

    # Grain growth near maturity
    if growth_progress > half_progress * 0.95
        grain_growth = storage_organ_growth
    else
        grain_growth = 0.0
    end

    crop.available_glucose_pool = glucose_pool

    return (leaf_growth          = leaf_growth, 
            stem_growth          = stem_growth,
            starch_growth        = starch_growth, 
            storage_organ_growth = storage_organ_growth, 
            root_growth          = root_growth, 
            grain_growth         = grain_growth,
            reserve_change       = glucose_consumed)
end

# ============ 6. Update biomass (Eq.113-118) ============
# Eq.113-118: Update dry weight of each organ
function update_biomass!(crop::CropState, leaf_growth::Float64, stem_growth::Float64,
                         starch_growth::Float64, storage_organ_growth::Float64, root_growth::Float64,
                         grain_growth::Float64, leaf_loss::Float64, root_loss::Float64,
                         starch_loss::Float64)
    # Eq.113: Leaf
    if (crop.leaf_biomass + leaf_growth - leaf_loss) > 0.0
        crop.leaf_biomass += leaf_growth - leaf_loss
    else
        crop.leaf_biomass += leaf_growth
    end
    
    crop.stem_biomass          += stem_growth                 # Eq.114: Stem
    crop.storage_organ_biomass += storage_organ_growth        # Eq.115: Storage organ (panicle)
    crop.root_biomass          += root_growth - root_loss     # Eq.116: Root
    crop.reserved_starch_pool  += starch_growth - starch_loss # Eq.117: Reserve
    crop.grain_biomass         += grain_growth                # Grain
    crop.dead_leaf_biomass     += leaf_loss                   # Dead leaf accumulation

    # Eq.136: Shoot biomass
    crop.shoot_biomass = crop.leaf_biomass + 
                         crop.stem_biomass +
                         crop.storage_organ_biomass + 
                         crop.reserved_starch_pool + 
                         crop.available_glucose_pool
end

# ============ 7. LAI calculation (Eq.137-138) ============
# Eq.137: L = (Wlef + Wglu) / Slw
# S_lw: specific leaf weight [kg / m^2]
# Eq.138: Slw = Slw,mx + (Slw,mn - Slw,mx) * exp(-kSlw * DVS)
function calc_LAI!(crop::CropState, leaf_weight_min::Float64, leaf_weight_max::Float64,
                   leaf_weight_decay_rate::Float64, co2_ppm::Float64, crop_name::String)
    growth_progress = crop.development_stage
    # Eq.138: Specific leaf weight [kg/m²] (SLW)
    s_lw = leaf_weight_max + (leaf_weight_min - leaf_weight_max) * exp(-leaf_weight_decay_rate * growth_progress)

    # CO2 down-regulation for C3 crops
    if crop_name != "Maize"
        s_lw = s_lw / ((0.856 * (1.0 + 1.035 * exp(-4.35e-3 * co2_ppm))) / (0.856 * (1.0 + 1.035 * exp(-4.35e-3 * 368.87))))
    end

    # Eq.137: LAI
    crop.LAI = (crop.leaf_biomass + crop.available_glucose_pool) / s_lw

    if crop.LAI > crop.LAI_max
        crop.LAI_max = crop.LAI
    end
end

# ============ 8. Crop height (Eq.139) ============
# Eq.139: hgt = haa * (DVS/hDVS) before heading; hgt = haa after heading
function calc_height!(crop::CropState, half_progress::Float64, max_crop_height::Float64, crop_name::String)
    growth_progress = crop.development_stage
    
    if crop_name == "Rice"
        if growth_progress < half_progress
            height_coeff_a1 = 0.439  
            height_coeff_b1 = 0.675            
            crop.crop_height = height_coeff_a1 * crop.LAI^height_coeff_b1 
        else
            height_coeff_a2 = 0.366      
            height_coeff_b2 = 0.318            
            crop.crop_height = height_coeff_a2 * crop.LAI^height_coeff_b2
        end
    else
        if growth_progress < half_progress
            crop.crop_height = max_crop_height * (growth_progress / half_progress)
        else
            crop.crop_height = max_crop_height
        end
    end
     
end

# ============ 9. Root length (Eq.140) ============
# Eq.140: zrt = min(zrt,mx, rrt * dt)
function calc_root_length!(crop::CropState, root_growth_rate::Float64,
                          max_root_length::Float64, Δt::Int)
    crop.root_length += root_growth_rate / seconds_per_day * Float64(Δt)
    crop.root_length = min(crop.root_length, max_root_length)
end

# ============ 10. Judge harvest ============
function judge_harvest!(crop::CropState, doy::Int, thermal_time_requirement::Float64,
                        five_day_avg_temp::Float64, harvest_temp_threshold::Float64,
                        heat_damage_threshold::Float64, harvest_index::Float64)
    # Harvest condition: maturity or cold-induced
    if (crop.is_planted == 1 && crop.accumulated_thermal_time > thermal_time_requirement) ||
       (crop.is_planted == 1 && five_day_avg_temp < harvest_temp_threshold && crop.development_stage >= 0.5)

        # Cold-adjusted harvest index
        harvest_index_cold_adjusted = harvest_index * (1.0 - (0.054 * crop.cold_damage_index^1.56) / 100.0)
        harvest_index_cold_adjusted = max(harvest_index_cold_adjusted, 0.0)

        # Heat-adjusted harvest index
        if crop.heat_damage_index > 0.0
            heat_damage_avg = crop.heat_damage_index / crop.n_heat_events - T_ice
        else
            heat_damage_avg = 0.0
        end
        harvest_index_heat_adjusted = harvest_index * 1.0 / (1.0 + exp(0.853 * (heat_damage_avg - heat_damage_threshold)))

        # Final harvest index
        harvest_index_final = min(harvest_index_cold_adjusted, harvest_index_heat_adjusted)

        # Record harvest DOY
        if doy > 0
            crop.harvest_doy = doy
        end

        # Eq.141: Yield
        crop.yield = crop.storage_organ_biomass * harvest_index_final

        # Save season peak values before reset
        crop.LAI_max_season = crop.LAI_max
        crop.shoot_biomass_at_harvest = crop.shoot_biomass

        # Full reset after harvest (matching Fortran JUDHVT subroutine)
        crop.development_stage = 0.0
        crop.is_planted = 0
        crop.has_emerged = false
        crop.is_grain_filling = 0
        # Reset all biomass variables
        crop.leaf_biomass = 0.0
        crop.stem_biomass = 0.0
        crop.storage_organ_biomass = 0.0
        crop.root_biomass = 0.0
        crop.reserved_starch_pool = 0.0
        crop.available_glucose_pool = 0.0
        crop.dead_leaf_biomass = 0.0
        crop.grain_biomass = 0.0
        crop.potential_grain_biomass = 0.0
        # Reset other state variables
        crop.LAI = 0.0
        crop.LAI_max = 0.0
        crop.crop_height = 0.0
        crop.root_length = 0.0
        crop.cold_damage_index = 0.0
        crop.heat_damage_index = 0.0
        crop.n_heat_events = 0
    end
end

# ============ 11. Cold damage index (CDI) ============
function calc_cold_damage_index!(crop::CropState, temperature::Float64,
                                 half_progress::Float64, panicle_progress_1::Float64,
                                 hour::Float64, Δt::Int,
                                 cold_damage_threshold::Float64)
    # Reset daily average at start of day
    if hour < Float64(Δt) / seconds_per_day * 24.0
        crop.daily_avg_temp = temperature
    else
        crop.daily_avg_temp += temperature
    end

    growth_progress = crop.development_stage

    # CDI active between flowering and grain filling
    if (panicle_progress_1 + 0.05) < growth_progress && growth_progress < ((1.0 - half_progress) * 0.2 + half_progress)
        if (24.0 - Float64(Δt) / seconds_per_day * 24.0) <= hour
            daily_avg_temp = crop.daily_avg_temp * Float64(Δt) / seconds_per_day

            if daily_avg_temp < (cold_damage_threshold + T_ice)
                crop.cold_damage_index += (cold_damage_threshold + T_ice) - daily_avg_temp
            end
        end
    end
end

# ============ 12. Heat damage index (HDI) ============
function calc_heat_damage_index!(crop::CropState, temperature::Float64,
                                 half_progress::Float64, hour::Float64, Δt::Int)
    # Track daily max temperature
    if hour < Float64(Δt) / seconds_per_day * 24.0
        crop.daily_max_temp = temperature
    elseif crop.daily_max_temp < temperature
        crop.daily_max_temp = temperature
    end

    growth_progress = crop.development_stage

    # HDI active near anthesis
    if half_progress * 0.96 < growth_progress && growth_progress < ((1.0 - half_progress) * 0.2 + half_progress)
        if (24.0 - Float64(Δt) / seconds_per_day * 24.0) <= hour
            crop.heat_damage_index += crop.daily_max_temp
            crop.n_heat_events += 1.0
        end
    end
end

# ============ 13. Specific leaf nitrogen [g N/m² leaf] (SLN, Eq.138 context) ============
function calc_specific_leaf_nitrogen!(crop::CropState, n_fertilizer::Float64,
                                      leaf_nitrogen_x1::Float64, leaf_nitrogen_x2::Float64, leaf_nitrogen_x3::Float64,
                                      leaf_nitrogen_max::Float64, leaf_nitrogen_min::Float64,
                                      leaf_nitrogen_sensitivity::Float64,
                                      co2_ppm::Float64, crop_name::String)
    growth_progress = crop.development_stage
    x1 = leaf_nitrogen_x1
    x2 = leaf_nitrogen_x2
    x3 = leaf_nitrogen_x3

    if crop_name == "Rice" || crop_name == "Wheat"
        y1 = 0.7742242627
        y2 = leaf_nitrogen_max - (leaf_nitrogen_max - leaf_nitrogen_min) * exp(-leaf_nitrogen_sensitivity * n_fertilizer)
        y3 = 0.5

        if growth_progress < 0.0
            sln = 0.0
        elseif growth_progress < x1
            sln = y1
        elseif growth_progress < x2
            sln = linear_interpolate(growth_progress, x1, y1, x2, y2)
        elseif growth_progress < x3
            sln = linear_interpolate(growth_progress, x2, y2, x3, y3)
        else
            sln = 0.0
        end

        # CO2 down-regulation (C3 crops)
        sln = sln * (1.037 - 8.33e-5 * co2_ppm) / (1.037 - 8.33e-5 * 368.87)

    elseif crop_name == "Soybean"
        # Soybean: 4-segment piecewise (matching Fortran CALSLN)
        # Breakpoints: 0, x1=SLNX1, x2=SLNX2, x3=SLNX3, 1.0
        y0 = 0.75
        y1 = 2.25
        y2 = 1.7
        y3 = (2.25 - 1.8) / 300.0 * n_fertilizer + 1.8
        y4 = 0.75

        if growth_progress < 0.0
            sln = 0.0
        elseif growth_progress < x1
            sln = linear_interpolate(growth_progress, 0.0, y0, x1, y1)
        elseif growth_progress < x2
            sln = linear_interpolate(growth_progress, x1, y1, x2, y2)
        elseif growth_progress < x3
            sln = linear_interpolate(growth_progress, x2, y2, x3, y3)
        else
            sln = linear_interpolate(growth_progress, x3, y3, 1.0, y4)
        end

        # CO2 down-regulation (C3)
        sln = sln * (1.037 - 8.33e-5 * co2_ppm) / (1.037 - 8.33e-5 * 368.87)

    elseif crop_name == "Maize"
        y1 = 0.5
        y2 = 2.1 - 1.4 * exp(-leaf_nitrogen_sensitivity * n_fertilizer)
        y3 = 0.0013 * n_fertilizer + 0.5295

        if growth_progress < 0.0
            sln = 0.0
        elseif growth_progress < x1
            sln = y1
        elseif growth_progress < x2
            sln = linear_interpolate(growth_progress, x1, y1, x2, y2)
        elseif growth_progress < x3
            sln = linear_interpolate(growth_progress, x2, y2, x3, y3)
        else
            sln = 0.0
        end
        # No CO2 down-regulation for Maize (C4 crop)

    else
        error("Unknown crop: $crop_name. Use Rice, Wheat, Soybean, or Maize")
    end

    crop.leaf_nitrogen = sln
end

# ============ Top-level crop step function ============
function crop_step!(crop::CropState;
                    doy::Int, hour::Float64, Δt::Int,
                    temperature::Float64, gpp::Float64, rsp::Float64,
                    planting_doy::Int, thermal_time_requirement::Float64, half_progress::Float64,
                    needs_vernalization::Int, base_temp::Float64,
                    optimal_temp::Float64, ceiling_temp::Float64,
                    vernalization_saturation::Float64,
                    k_leaf_convert::Float64, k_stem_convert::Float64,
                    k_root_convert::Float64, k_grain_convert::Float64,
                    fraction_starch_reserve::Float64,
                    
                    # progress key breaks for different organs
                    shoot_progress_1::Float64, 
                    shoot_alloc_ratio_1::Float64, shoot_progress_2::Float64,

                    leaf_alloc_ratio_0::Float64, 
                    leaf_progress_1::Float64, leaf_alloc_ratio_1::Float64,
                    leaf_progress_2::Float64, leaf_alloc_ratio_2::Float64,

                    panicle_progress_1::Float64, panicle_alloc_ratio_1::Float64,
                    panicle_progress_2::Float64, panicle_alloc_ratio_2::Float64,
                    panicle_progress_3::Float64, panicle_alloc_ratio_3::Float64,

                    dead_prgress_1::Float64, dead_ratio_1::Float64,
                    dead_prgress_2::Float64, dead_ratio_2::Float64,
                    dead_prgress_3::Float64, dead_ratio_3::Float64,
                    
                    leaf_nitrogen_x1::Float64, leaf_nitrogen_x2::Float64, leaf_nitrogen_x3::Float64,
                    leaf_nitrogen_max::Float64, leaf_nitrogen_min::Float64,
                    leaf_nitrogen_sensitivity::Float64,

                    LAI_threshold_grain::Float64, k_leaf_loss::Float64,
                    leaf_weight_min::Float64, leaf_weight_max::Float64,
                    leaf_weight_decay_rate::Float64,
                    max_crop_height::Float64,
                    root_growth_rate::Float64, max_root_length::Float64,
                    
                    n_fertilizer::Float64, co2_ppm::Float64,
                    is_irrigated::Int,
                    cold_damage_threshold::Float64,
                    heat_damage_threshold::Float64,
                    harvest_index::Float64,
                    harvest_temp_threshold::Float64,
                    five_day_temp_buffer::Vector{Float64},
                    five_day_temp_count::Int,
                    crop_name::String)

    # 1. Judge planting
    judge_planting!(crop, doy, hour, planting_doy, is_irrigated)

    if crop.is_planted > 0
        # 2. Development stage (CRODVS)
        calc_development_stage!(crop, temperature, Δt, thermal_time_requirement,
                                needs_vernalization, base_temp,
                                optimal_temp, ceiling_temp,
                                vernalization_saturation)

        if crop.has_emerged
            # 4. Loss rates
            losses = calc_loss_rates(crop, Δt, half_progress,
                                     dead_prgress_1, dead_ratio_1, dead_prgress_2, dead_ratio_2,
                                     dead_prgress_3, dead_ratio_3)

            # 5. Growth and partitioning
            growth = calc_growth_partitioning!(crop, gpp - rsp, losses.starch_loss, Δt,
                                              half_progress, k_leaf_convert, k_stem_convert, k_root_convert,
                                              k_grain_convert, fraction_starch_reserve,
                                              shoot_progress_1, shoot_alloc_ratio_1, shoot_progress_2,
                                              leaf_alloc_ratio_0, leaf_progress_1, leaf_alloc_ratio_1,
                                              leaf_progress_2, leaf_alloc_ratio_2,
                                              panicle_progress_1, panicle_alloc_ratio_1,
                                              panicle_progress_2, panicle_alloc_ratio_2,
                                              panicle_progress_3, panicle_alloc_ratio_3,
                                              LAI_threshold_grain, k_leaf_loss, crop_name)

            # 6. Update biomass
            update_biomass!(crop, growth.leaf_growth, growth.stem_growth,
                           growth.starch_growth, growth.storage_organ_growth,
                           growth.root_growth, growth.grain_growth,
                           losses.leaf_loss, losses.root_loss, losses.starch_loss)

            # 7. LAI
            calc_LAI!(crop, leaf_weight_min, leaf_weight_max, leaf_weight_decay_rate, co2_ppm, crop_name)

            # 8. Height
            calc_height!(crop, half_progress, max_crop_height, crop_name)

            # 9. Root length
            calc_root_length!(crop, root_growth_rate, max_root_length, Δt)

            # 11. Cold damage
            calc_cold_damage_index!(crop, temperature, half_progress, panicle_progress_1,
                                   hour, Δt, cold_damage_threshold)

            # 12. Heat damage
            calc_heat_damage_index!(crop, temperature, half_progress, hour, Δt)
        end

        # 10. Judge harvest
        five_day_avg = sum(five_day_temp_buffer[1:five_day_temp_count]) / Float64(five_day_temp_count)
        judge_harvest!(crop, doy, thermal_time_requirement, five_day_avg,
                       harvest_temp_threshold, heat_damage_threshold, harvest_index)

        # 3. Judge emergence (JUDEMR: after JUDHVT, same order as Fortran)
        # Fortran calls JUDHVT first, then JUDEMR. JUDEMR checks DVS>0.012
        # and sets WLF=1,WST=1,WRT=1,WAR=0.5. This must happen AFTER CRODVS
        # has updated DVS, so emergence can be detected in the same timestep
        # that DVS crosses the 0.012 threshold.
        judge_emergence!(crop)
    end

    # 13. Specific leaf nitrogen (always calculated)
    calc_specific_leaf_nitrogen!(crop, n_fertilizer, leaf_nitrogen_x1, leaf_nitrogen_x2, leaf_nitrogen_x3,
                                 leaf_nitrogen_max, leaf_nitrogen_min, leaf_nitrogen_sensitivity,
                                 co2_ppm, crop_name)

    return crop
end
