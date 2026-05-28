# MATCRO - Unified Julia Crop Model
# Combines: constants, tinterp, radiation, photosynthesis (C3/C4), crop, soil, io
# Supports both CSV and NetCDF input formats
# Usage: julia matcro.jl config.toml

# ============================================================
# Module includes (in dependency order)
# ============================================================
include("01_constants.jl")
include("02_io.jl")
include("02_1_point_io.jl")
include("02_2_spatial_io.jl")
include("03_hour_interpolation.jl")
include("04_radiation.jl")
include("05_1_photosynthesis_C3.jl")
include("05_2_photosynthesis_C4.jl")
include("05_photosynthesis.jl")
include("06_crop.jl")
include("07_soil.jl")

using Printf

# ============================================================
# run_simulation — main entry point
# ============================================================
function run_simulation(config_path::String)
    println("=" ^ 60)
    println("                   MATCRO (Julia Version)")
    println("=" ^ 60)
    flush(stdout)

    # 1. Read configuration
    config = read_config(config_path)
    println("  Config:")
    println("    Crop: ", config.crop_name)
    println("    Period: $(config.start_year)-$(config.end_year)")
    println("    Input format: ", config.input_format)

    # 2. Read crop parameters
    params = read_crop_params(config.crop_param_file)
    println("    Parameters loaded from: ", config.crop_param_file)

    # 3. Load forcing data
    if config.input_format == "point"
        println("  Loading CSV forcing: ", config.csv_path)
        forcing_data = read_forcing_csv(config.csv_path)
    elseif config.input_format == "raster"
        println("  Will load NetCDF forcing per year")
        forcing_data = nothing
    else
        error("Unknown input format: $(config.input_format)")
    end

    # 3b. Spatial parallel mode (NetCDF)
    if config.input_format == "raster"
        return run_spatial_simulation(config)
    end

    # 4. Initialize crop state
    crop = CropState()

    # 5. Initialize soil water (5 layers, saturated → clamped to porosity on first SOIL call)
    # Fortran initializes WSL(1:NSL) = 1.D0 (saturated); SOIL clamps to porosity
    layer_water = [1.0, 1.0, 1.0, 1.0, 1.0]
    water_stress = 1.0

    # 6. Five-day temperature buffer (pre-fill with 295.15 K = 22°C, matching Fortran)
    buf_len = 86400 ÷ config.time_step * 5
    five_day_buffer = fill(0.0, buf_len)
    five_day_count = 0

    # 7. Output storage
    daily_records = NamedTuple[]
    yearly_results = NamedTuple[]

    # 8. Year loop
    for year in config.start_year:config.end_year
        println("\n  Year $year ...")

        # Read CO2
        co2_ppm = read_co2(config, year)

        # Load this year's forcing
        if config.input_format == "point"
            year_forcing = forcing_data
        elseif config.input_format == "raster"
            year_forcing = Dict{Int,Dict{Int,DailyForcing}}()
            year_forcing[year] = read_forcing_netcdf(config, year)
        else
            error!("Invalid input format: $(config.input_format)")
        end

        # Determine DOY range
        if config.start_year == config.end_year
            stdd = config.start_doy; endd = config.end_doy
        elseif year == config.start_year
            stdd = config.start_doy; endd = 365
        elseif year == config.end_year
            stdd = 1; endd = config.end_doy
        else
            stdd = 1; endd = 365
        end

        # DOY loop
        for doy in stdd:endd
            # Adjacent-day forcing for temperature interpolation
            # At year boundaries (DOY 1/365), try previous/next year's data if available
            f_today = get_forcing(year_forcing, year, doy)
            if doy == 1
                f_prev = get_forcing(year_forcing, year - 1, 365)
                f_prev === nothing && (f_prev = f_today)
            else
                f_prev = get_forcing(year_forcing, year, doy - 1)
            end
            if doy == 365
                f_next = get_forcing(year_forcing, year + 1, 1)
                f_next === nothing && (f_next = f_today)
            else
                f_next = get_forcing(year_forcing, year, doy + 1)
            end

            if f_today === nothing
                continue
            end

            # Integrated sin(solar elevation)
            int_sinb = calc_int_sinb(doy, config.latitude, config.time_step)

            # Hourly loop
            n_steps = 86400 ÷ config.time_step

            # Debug output file for module-level comparison
            debug_path = joinpath(config.output_dir, "debug_module_outputs.csv")
            mkpath(config.output_dir)
            if !isfile(debug_path)
                debug_file = open(debug_path, "w")
                println(debug_file, "year,doy,hour,TMP,RSD,SHM,WND,PRS,WSTRS_in,LAI,DVS,GPP,RSP,TSP,ROT,QPARSNLF,QPARSHLF,VMXSNLF,VMXSHLF,LAISN,LAISH,WLF,WST,WSO,WRT,WSR,WAR,WDL")
            else
                debug_file = open(debug_path, "a")
            end

            for ihour in 1:n_steps
                hour = (Float64(ihour) - 0.5) * Float64(config.time_step) / 3600.0

                # ----- Time interpolation (daily → hourly) -----
                hourly = interpolate_time(;
                    doy=doy, prev_doy=(doy == 1 ? 365 : doy - 1), next_doy=(doy == 365 ? 1 : doy + 1), hour=hour,
                    lat=config.latitude, Δt=config.time_step,
                    tmax_prev=f_prev.tmax, tmax=f_today.tmax, tmax_next=f_next.tmax,
                    tmin_prev=f_prev.tmin, tmin=f_today.tmin, tmin_next=f_next.tmin,
                    radiation=f_today.radiation, precip=f_today.precip,
                    humidity=f_today.humidity, wind=f_today.wind,
                    pressure=f_today.pressure, ozone=f_today.ozone,
                    int_sinb=int_sinb,
                    wind_height=get(config.raster_vars, "wnd", Dict("height"=>10.0))["height"]
                )

                tmp_K = hourly.temperature  # already in K from tinterp
                wnd = max(hourly.wind, 0.001)
                prc = hourly.precipitation
                rsd = hourly.radiation
                shm = hourly.humidity
                prs = hourly.pressure

                # ----- 5-day temperature buffer -----
                if five_day_count < buf_len
                    five_day_count += 1
                    five_day_buffer[five_day_count] = tmp_K - T_ice
                else
                    five_day_buffer[1:(five_day_count-1)] = five_day_buffer[2:five_day_count]
                    five_day_buffer[five_day_count] = tmp_K - T_ice
                end

                # ----- RAD (radiation transfer) -----
                rad_result = calc_radiation(;
                    leaf_nitrogen=crop.leaf_nitrogen,
                    kn=params.k_nitrogen,
                    shortwave_radiation=rsd,
                    LAI=crop.LAI,
                    RLFv=params.leaf_PAR_reflectance,
                    TLFv=params.leaf_PAR_transmittance,
                    RLFn=params.leaf_NIR_reflectance,
                    TLFn=params.leaf_NIR_transmittance,
                    lat=config.latitude,
                    doy=doy,
                    hour=hour,
                    crop_name=config.crop_name,
                    development_stage=crop.development_stage
                )

                # ----- PHSYN (photosynthesis → GPP, RSP, TSP) -----
                phsyn_result = calc_photosynthesis(;
                    Qp_sunlit=rad_result.PAR_abs_sunlit_leaf,
                    Qp_shade=rad_result.PAR_abs_shade_leaf,
                    Vmax25_sunlit=rad_result.Vmax_sunlit_leaf,
                    Vmax25_shade=rad_result.Vmax_shade_leaf,
                    LAI_sunlit=rad_result.LAI_sunlit,
                    LAI_shade=rad_result.LAI_shade,
                    leaf_temperature=tmp_K,
                    wind_speed=wnd,
                    specific_humidity=shm,
                    pressure=prs,
                    co2_ppm=co2_ppm,
                    water_stress=water_stress,
                    crop_height=crop.crop_height,
                    EFFCON=params.quantum_efficiency,
                    atheta=params.a_theta,
                    btheta=params.b_theta,
                    m_H2O=params.m_H2O,
                    b_H2O=params.b_H2O,
                    crop_name=config.crop_name
                )

                gpp = phsyn_result.gpp
                rsp = phsyn_result.rsp
                tsp = phsyn_result.tsp

                # ----- CROP (growth simulation) -----
                crop_step!(crop;
                    doy=doy, hour=hour, Δt=config.time_step,
                    temperature=tmp_K, gpp=gpp, rsp=rsp,
                    planting_doy=config.planting_doy + Int(params.planting_offset),
                    thermal_time_requirement=config.thermal_time_requirement,
                    half_progress=params.half_progress,
                    needs_vernalization=params.needs_vernalization,
                    base_temp=params.base_temp,
                    optimal_temp=params.optimal_temp,
                    ceiling_temp=params.ceiling_temp,
                    vernalization_saturation=params.vernalization_saturation,
                    k_leaf_convert=params.k_leaf_convert,
                    k_stem_convert=params.k_stem_convert,
                    k_root_convert=params.k_root_convert,
                    k_grain_convert=params.k_grain_convert,
                    fraction_starch_reserve=params.fraction_starch_reserve,
                    shoot_progress_1=params.shoot_progress_1,
                    shoot_alloc_ratio_1=params.shoot_alloc_ratio_1,
                    shoot_progress_2=params.shoot_progress_2,
                    leaf_alloc_ratio_0=params.leaf_alloc_ratio_0,
                    leaf_progress_1=params.leaf_progress_1,
                    leaf_alloc_ratio_1=params.leaf_alloc_ratio_1,
                    leaf_progress_2=params.leaf_progress_2,
                    leaf_alloc_ratio_2=params.leaf_alloc_ratio_2,
                    panicle_progress_1=params.panicle_progress_1,
                    panicle_alloc_ratio_1=params.panicle_alloc_ratio_1,
                    panicle_progress_2=params.panicle_progress_2,
                    panicle_alloc_ratio_2=params.panicle_alloc_ratio_2,
                    panicle_progress_3=params.panicle_progress_3,
                    panicle_alloc_ratio_3=params.panicle_alloc_ratio_3,
                    dead_prgress_1=params.dead_progress_1,
                    dead_ratio_1=params.dead_ratio_1,
                    dead_prgress_2=params.dead_progress_2,
                    dead_ratio_2=params.dead_ratio_2,
                    dead_prgress_3=params.dead_progress_3,
                    dead_ratio_3=params.dead_ratio_3,
                    leaf_nitrogen_x1=params.leaf_nitrogen_x1,
                    leaf_nitrogen_x2=params.leaf_nitrogen_x2,
                    leaf_nitrogen_x3=params.leaf_nitrogen_x3,
                    leaf_nitrogen_max=params.leaf_nitrogen_max,
                    leaf_nitrogen_min=params.leaf_nitrogen_min,
                    leaf_nitrogen_sensitivity=params.leaf_nitrogen_sensitivity,
                    LAI_threshold_grain=params.LAI_threshold_grain,
                    k_leaf_loss=params.k_leaf_loss,
                    leaf_weight_asymptote=params.leaf_weight_asymptote,
                    leaf_weight_intercept=params.leaf_weight_intercept,
                    leaf_weight_decay_rate=params.leaf_weight_decay_rate,
                    max_crop_height=params.height_coeff_a1,
                    root_growth_rate=params.root_growth_rate,
                    max_root_length=params.max_root_length,
                    n_fertilizer=config.n_fertilizer,
                    co2_ppm=co2_ppm,
                    is_irrigated=config.is_irrigated,
                    cold_damage_threshold=params.cold_damage_threshold,
                    heat_damage_threshold=params.heat_damage_threshold,
                    harvest_index=params.harvest_index,
                    harvest_temp_threshold=params.harvest_temp_threshold,
                    five_day_temp_buffer=five_day_buffer,
                    five_day_temp_count=five_day_count,
                    crop_name=config.crop_name
                )

                # ----- SOIL (water balance) -----
                # Convert TSP from kg/m2/s to W/m2 for soil module
                tsp_W = tsp * L_vaporization

                soil_result = calc_soil_water(;
                    layer_water=layer_water,
                    transpiration=tsp_W,
                    W2SF=prc,
                    z_rt=crop.root_length,
                    is_irrigated=config.is_irrigated,
                    Δt=config.time_step,
                    soil_type_i=config.soil_type,
                    temperature=tmp_K,
                    pressure=prs,
                    wind_speed=wnd,
                    specific_humidity=shm,
                    crop_height=crop.crop_height,
                    is_planted=crop.is_planted,
                    crop_name=config.crop_name
                )
                water_stress = soil_result.water_stress

                # Write hourly debug output (all modules, matching Fortran debug_module_outputs.txt)
                #= Printf.@printf(debug_file, "%d,%d,%.1f,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e\n",
                    year, doy, hour,
                    tmp_K, rsd, shm, wnd, prs,
                    water_stress_before_phsyn, crop.LAI, crop.development_stage,
                    gpp, rsp, tsp,
                    crop.root_length,
                    rad_result.PAR_abs_sunlit_leaf, rad_result.PAR_abs_shade_leaf,
                    rad_result.Vmax_sunlit_leaf, rad_result.Vmax_shade_leaf,
                    rad_result.LAI_sunlit, rad_result.LAI_shade,
                    crop.leaf_biomass, crop.stem_biomass, crop.storage_organ_biomass,
                    crop.root_biomass, crop.reserved_starch_pool,
                    crop.available_glucose_pool, crop.dead_leaf_biomass) =#

            end # hourly loop

            close(debug_file)

            # Record daily state
            push!(daily_records, (
                year=year, doy=doy,
                yield=crop.yield,
                LAI=crop.LAI,
                development_stage=crop.development_stage,
                biomass_aboveground=crop.shoot_biomass,
                leaf_biomass=crop.leaf_biomass,
                stem_biomass=crop.stem_biomass,
                root_biomass=crop.root_biomass,
                storage_organ_biomass=crop.storage_organ_biomass,
                crop_height=crop.crop_height,
                root_length=crop.root_length,
                water_stress=water_stress,
                harvest_doy=crop.harvest_doy,
            ))

        end # DOY loop

        # Record yearly result
        push!(yearly_results, (
            year=year,
            yield=crop.yield,
            LAI_max=crop.LAI_max_season,
            biomass_aboveground=crop.shoot_biomass_at_harvest,
            harvest_doy=crop.harvest_doy,
        ))

        @printf("    Yield: %.2f kg/ha, Harvest DOY: %d, LAI_max: %.2f, Biomass_aboveground: %.2f kg/ha\n",
                crop.yield, crop.harvest_doy, crop.LAI_max_season, crop.shoot_biomass_at_harvest)

    end # Year loop

    # 9. Write output
    mkpath(config.output_dir)
    write_output_csv(yearly_results, joinpath(config.output_dir, "yield_summary.csv"))
    write_daily_csv(daily_records, joinpath(config.output_dir, "daily_output.csv"))

    println("\n" * "=" ^ 60)
    println("  Simulation Complete")
    println("  Output saved to: ", abspath(config.output_dir))
    println("=" ^ 60)
    return yearly_results, daily_records
end

# ============================================================
# run_spatial_simulation — NetCDF spatial parallel mode
# Reads entire spatial fields, runs pixels in parallel with Threads
# ============================================================

function run_spatial_simulation(config::Config)
    n_actual_threads = Threads.nthreads()
    println("\n  Spatial parallel mode: $n_actual_threads thread(s)")

    # 1. Read crop parameters
    params = read_crop_params(config.crop_param_file)

    # 2. Get grid info
    lats, lons = get_grid_info(config)
    n_lon = length(lons)
    n_lat = length(lats)
    n_pixels = n_lon * n_lat
    println("  Grid: $n_lon x $n_lat = $n_pixels pixels")

    # 3. Create boundary mask (if boundary file is specified)
    boundary_mask = nothing
    if !isempty(config.boundary_file)
        if isfile(config.boundary_file)
            println("  Loading boundary: ", config.boundary_file)
            boundary_mask = create_boundary_mask(lons, lats, config.boundary_file; buffer_deg=config.boundary_buffer)
            n_in_boundary = sum(boundary_mask)
            println("  Boundary mask: $n_in_boundary / $n_pixels pixels within boundary")
        else
            println("  [WARN] Boundary file not found: $(config.boundary_file), running all pixels")
        end
    end

    # 4. Year loop
    years = collect(config.start_year:config.end_year)
    yield_3d = Array{Float64,3}(undef, n_lon, n_lat, length(years))
    harvest_doy_3d = Array{Float64,3}(undef, n_lon, n_lat, length(years))
    LAI_max_3d = Array{Float64,3}(undef, n_lon, n_lat, length(years))
    biomass_aboveground_3d = Array{Float64,3}(undef, n_lon, n_lat, length(years))

    for (i_year, year) in enumerate(years)
        println("\n  Year $year ...")
        co2_ppm = read_co2(config, year)

        # Read spatial forcing
        spatial_forcing = read_forcing_netcdf_spatial(config, year)
        n_days = length(spatial_forcing)

        # Read management params for this year (from NC files or defaults)
        planting_doy = load_management_param(config, "planting_doy", year, n_lon, n_lat; lats=lats, lons=lons)
        is_irrigated  = load_management_param(config, "is_irrigated", year, n_lon, n_lat; lats=lats, lons=lons)
        soil_type  = load_management_param(config, "soil_type", year, n_lon, n_lat; lats=lats, lons=lons)
        n_fertilizer  = load_management_param(config, "n_fertilizer", year, n_lon, n_lat; lats=lats, lons=lons)
        thermal_time_requirement  = load_management_param(config, "thermal_time_requirement", year, n_lon, n_lat; lats=lats, lons=lons)

        # Parallel pixel loop
        indices = [(i_lon, i_lat) for i_lon in 1:n_lon for i_lat in 1:n_lat]
        results = Vector{NamedTuple}(undef, length(indices))

        Threads.@threads for index in eachindex(indices)
            i_lon, i_lat = indices[index]
            # Skip pixels outside boundary
            if boundary_mask !== nothing && !boundary_mask[i_lon, i_lat]
                results[index] = (yield=NaN, harvest_doy=NaN, LAI_max=NaN, biomass_aboveground=NaN)
                continue
            end
            lat = lats[i_lat]
            lon = lons[i_lon]
            try
                results[index] = run_pixel_spatial(
                    config, params, year, lat, lon,
                    spatial_forcing, n_days, co2_ppm,
                    Float64(planting_doy[i_lon, i_lat]), Float64(is_irrigated[i_lon, i_lat]),
                    Float64(soil_type[i_lon, i_lat]), Float64(n_fertilizer[i_lon, i_lat]),
                    Float64(thermal_time_requirement[i_lon, i_lat]);
                    i_lon=i_lon, i_lat=i_lat)
            catch
                results[index] = (yield=NaN, harvest_doy=NaN, LAI_max=NaN, biomass_aboveground=NaN)
            end
        end

        for (index, (i_lon, i_lat)) in enumerate(indices)
            yield_3d[i_lon, i_lat, i_year] = results[index].yield
            harvest_doy_3d[i_lon, i_lat, i_year] = results[index].harvest_doy
            LAI_max_3d[i_lon, i_lat, i_year] = results[index].LAI_max
            biomass_aboveground_3d[i_lon, i_lat, i_year] = results[index].biomass_aboveground
        end

        valid_mask = .!isnan.(yield_3d[:, :, i_year])
        if any(valid_mask)
            v_yield = yield_3d[:, :, i_year][valid_mask]
            v_harvest = harvest_doy_3d[:, :, i_year][valid_mask]
            v_lai = LAI_max_3d[:, :, i_year][valid_mask]
            v_biomass = biomass_aboveground_3d[:, :, i_year][valid_mask]
            @printf("    Mean: Yield=%.2f kg/ha, Harvest DOY=%.1f, LAI_max=%.2f, Biomass_aboveground=%.2f kg/ha\n",
                    sum(v_yield)/length(v_yield), sum(v_harvest)/length(v_harvest),
                    sum(v_lai)/length(v_lai), sum(v_biomass)/length(v_biomass))
        end
    end

    # 5. Write output
    mkpath(config.output_dir)
    # Write per-year TIF
    for (i_year, year) in enumerate(years)
        yield_path = joinpath(config.output_dir, "yield_$(year).tif")
        harvest_path = joinpath(config.output_dir, "harvest_doy_$(year).tif")
        LAI_max_path = joinpath(config.output_dir, "LAI_max_$(year).tif")
        biomass_path = joinpath(config.output_dir, "biomass_aboveground_$(year).tif")
        write_float64_tif(yield_3d[:, :, i_year], lats, lons, year, yield_path; unit="kg/ha")
        write_harvest_doy_tif(harvest_doy_3d[:, :, i_year], lats, lons, year, harvest_path)
        write_float64_tif(LAI_max_3d[:, :, i_year], lats, lons, year, LAI_max_path; unit="m2/m2")
        write_float64_tif(biomass_aboveground_3d[:, :, i_year], lats, lons, year, biomass_path; unit="kg/ha")
    end

    println("\n" * "=" ^ 60)
    println("  Spatial Simulation Complete")
    println("  Output saved to: ", abspath(config.output_dir))
    println("=" ^ 60)
    return yield_3d
end

# ============================================================
# run_pixel_spatial — run MATCRO for a single pixel in spatial mode
# ============================================================
function run_pixel_spatial(config::Config, params::CropParameters,
                           year::Int, lat::Float64, lon::Float64,
                           spatial_forcing::Dict, n_days::Int, co2_ppm::Float64,
                           p_planting_doy::Float64, p_is_irrigated::Float64,
                           p_soil_type::Float64, p_n_fertilizer::Float64,
                           p_thermal_time::Float64;
                           i_lon::Int=1, i_lat::Int=1)
    # Initialize state
    crop = CropState()
    layer_water = [1.0, 1.0, 1.0, 1.0, 1.0]
    water_stress = 1.0
    Δt = config.time_step
    buf_len = 86400 ÷ Δt * 5
    five_day_buffer = fill(0.0, buf_len)
    five_day_count = 0

    for doy in 1:n_days
        weather_spatial = spatial_forcing[doy]
        weather_today = DailyForcing(;
            doy=doy,
            tmax=weather_spatial.tmax[i_lon, i_lat], tmin=weather_spatial.tmin[i_lon, i_lat],
            radiation=weather_spatial.radiation[i_lon, i_lat], precip=weather_spatial.precip[i_lon, i_lat],
            humidity=weather_spatial.humidity[i_lon, i_lat], wind=weather_spatial.wind[i_lon, i_lat],
            pressure=weather_spatial.pressure[i_lon, i_lat],
        )
        int_sinb = calc_int_sinb(doy, lat, Δt)
        n_steps = 86400 ÷ Δt

        for ihour in 1:n_steps
            hour = (Float64(ihour) - 0.5) * Float64(Δt) / 3600.0

            hourly = interpolate_time(;
                doy=doy, prev_doy=(doy == 1 ? n_days : doy - 1),
                next_doy=(doy == n_days ? 1 : doy + 1),
                hour=hour, lat=lat, Δt=Δt,
                tmax_prev=weather_today.tmax, tmax=weather_today.tmax, tmax_next=weather_today.tmax,
                tmin_prev=weather_today.tmin, tmin=weather_today.tmin, tmin_next=weather_today.tmin,
                radiation=weather_today.radiation, precip=weather_today.precip,
                humidity=weather_today.humidity, wind=weather_today.wind,
                pressure=weather_today.pressure, ozone=weather_today.ozone,
                int_sinb=int_sinb, wind_height=10.0,
            )

            temperature_K = hourly.temperature
            wind_hourly = max(hourly.wind, 0.001)
            precip_hourly = hourly.precipitation
            radiation_hourly = hourly.radiation
            humidity_hourly = hourly.humidity
            pressure_hourly = hourly.pressure

            # 5-day temp buffer
            if five_day_count < buf_len
                five_day_count += 1
                five_day_buffer[five_day_count] = temperature_K - T_ice
            else
                five_day_buffer[1:(five_day_count-1)] = five_day_buffer[2:five_day_count]
                five_day_buffer[five_day_count] = temperature_K - T_ice
            end

            # Radiation
            rad_result = calc_radiation(;
                leaf_nitrogen=crop.leaf_nitrogen, kn=params.k_nitrogen,
                shortwave_radiation=radiation_hourly, LAI=crop.LAI,
                RLFv=params.leaf_PAR_reflectance, TLFv=params.leaf_PAR_transmittance,
                RLFn=params.leaf_NIR_reflectance, TLFn=params.leaf_NIR_transmittance,
                lat=lat, doy=doy, hour=hour,
                crop_name=config.crop_name, development_stage=crop.development_stage,
            )

            # Photosynthesis
            phsyn_result = calc_photosynthesis(;
                Qp_sunlit=rad_result.PAR_abs_sunlit_leaf,
                Qp_shade=rad_result.PAR_abs_shade_leaf,
                Vmax25_sunlit=rad_result.Vmax_sunlit_leaf,
                Vmax25_shade=rad_result.Vmax_shade_leaf,
                LAI_sunlit=rad_result.LAI_sunlit, LAI_shade=rad_result.LAI_shade,
                leaf_temperature=temperature_K, wind_speed=wind_hourly,
                specific_humidity=humidity_hourly, pressure=pressure_hourly,
                co2_ppm=co2_ppm, water_stress=water_stress,
                crop_height=crop.crop_height,
                EFFCON=params.quantum_efficiency,
                atheta=params.a_theta, btheta=params.b_theta,
                m_H2O=params.m_H2O, b_H2O=params.b_H2O,
                crop_name=config.crop_name,
            )

            gpp = phsyn_result.gpp
            rsp = phsyn_result.rsp
            tsp = phsyn_result.tsp

            # Crop
            crop_step!(crop;
                doy=doy, hour=hour, Δt=Δt,
                temperature=temperature_K, gpp=gpp, rsp=rsp,
                planting_doy=Int(p_planting_doy) + Int(params.planting_offset),
                thermal_time_requirement=p_thermal_time,
                half_progress=params.half_progress,
                needs_vernalization=params.needs_vernalization,
                base_temp=params.base_temp, optimal_temp=params.optimal_temp,
                ceiling_temp=params.ceiling_temp,
                vernalization_saturation=params.vernalization_saturation,
                k_leaf_convert=params.k_leaf_convert, k_stem_convert=params.k_stem_convert,
                k_root_convert=params.k_root_convert, k_grain_convert=params.k_grain_convert,
                fraction_starch_reserve=params.fraction_starch_reserve,
                shoot_progress_1=params.shoot_progress_1, shoot_alloc_ratio_1=params.shoot_alloc_ratio_1,
                shoot_progress_2=params.shoot_progress_2,
                leaf_alloc_ratio_0=params.leaf_alloc_ratio_0,
                leaf_progress_1=params.leaf_progress_1, leaf_alloc_ratio_1=params.leaf_alloc_ratio_1,
                leaf_progress_2=params.leaf_progress_2, leaf_alloc_ratio_2=params.leaf_alloc_ratio_2,
                panicle_progress_1=params.panicle_progress_1, panicle_alloc_ratio_1=params.panicle_alloc_ratio_1,
                panicle_progress_2=params.panicle_progress_2, panicle_alloc_ratio_2=params.panicle_alloc_ratio_2,
                panicle_progress_3=params.panicle_progress_3, panicle_alloc_ratio_3=params.panicle_alloc_ratio_3,
                dead_prgress_1=params.dead_progress_1, dead_ratio_1=params.dead_ratio_1,
                dead_prgress_2=params.dead_progress_2, dead_ratio_2=params.dead_ratio_2,
                dead_prgress_3=params.dead_progress_3, dead_ratio_3=params.dead_ratio_3,
                leaf_nitrogen_x1=params.leaf_nitrogen_x1, leaf_nitrogen_x2=params.leaf_nitrogen_x2,
                leaf_nitrogen_x3=params.leaf_nitrogen_x3,
                leaf_nitrogen_max=params.leaf_nitrogen_max, leaf_nitrogen_min=params.leaf_nitrogen_min,
                leaf_nitrogen_sensitivity=params.leaf_nitrogen_sensitivity,
                LAI_threshold_grain=params.LAI_threshold_grain, k_leaf_loss=params.k_leaf_loss,
                leaf_weight_asymptote=params.leaf_weight_asymptote,
                leaf_weight_intercept=params.leaf_weight_intercept,
                leaf_weight_decay_rate=params.leaf_weight_decay_rate,
                max_crop_height=2.5,
                root_growth_rate=params.root_growth_rate, max_root_length=params.max_root_length,
                n_fertilizer=p_n_fertilizer, co2_ppm=co2_ppm,
                is_irrigated=Int(p_is_irrigated),
                cold_damage_threshold=params.cold_damage_threshold,
                heat_damage_threshold=params.heat_damage_threshold,
                harvest_index=params.harvest_index, harvest_temp_threshold=params.harvest_temp_threshold,
                five_day_temp_buffer=five_day_buffer, five_day_temp_count=five_day_count,
                crop_name=config.crop_name,
            )

            # Soil
            tsp_W = tsp * L_vaporization
            soil_result = calc_soil_water(;
                layer_water=layer_water, transpiration=tsp_W,
                W2SF=precip_hourly, depth_root=crop.root_length,
                is_irrigated=Int(p_is_irrigated), Δt=Δt,
                soil_type_i=Int(p_soil_type), temperature=temperature_K,
                pressure=pressure_hourly, wind_speed=wind_hourly,
                specific_humidity=humidity_hourly, crop_height=crop.crop_height,
                is_planted=crop.is_planted, crop_name=config.crop_name,
            )
            water_stress = soil_result.water_stress
        end
    end

    return (yield=crop.yield, harvest_doy=crop.harvest_doy,
            LAI_max=crop.LAI_max_season, biomass_aboveground=crop.shoot_biomass_at_harvest)
end

# Command-line entry point (skip when included from other scripts)
# Check if this file is being run as the main program
const _RUN_AS_MAIN = length(ARGS) > 0 && endswith(lowercase(ARGS[1]), ".toml")
if _RUN_AS_MAIN
    run_simulation(ARGS[1])
end