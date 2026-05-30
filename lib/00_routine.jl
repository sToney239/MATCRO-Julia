# 00_routine.jl — Core simulation routines
# Contains: ManagementParams struct, init_simulation_state, run_hourly_step!,
#            run_point_simulation, run_spatial_simulation, run_pixel_spatial

using Printf

# ============================================================
# ManagementParams — unified struct for per-site management
# ============================================================
struct ManagementParams
    planting_doy::Int
    is_irrigated::Int
    soil_type::Int
    n_fertilizer::Float64
    thermal_time_requirement::Float64
    max_crop_height::Float64
    wind_height::Float64
end

# ============================================================
# init_simulation_state — shared initialization
# ============================================================
function init_simulation_state(Δt::Int)
    crop = CropState()
    layer_water = [1.0, 1.0, 1.0, 1.0, 1.0]
    water_stress = 1.0
    buffer_len = 86400 ÷ Δt * 5   # 86400 sec/day = 24 h/day * 3600 sec/h, and default 5 days as buffer
    five_day_buffer = fill(0.0, buffer_len)
    five_day_count = 0
    return (crop, layer_water, water_stress, five_day_buffer, five_day_count)
end

# ============================================================
# run_hourly_step! — single hourly time step
# Runs: interpolate_time → 5-day buffer → radiation → photosynthesis
#       → crop_step! → soil_water → water_stress update
# ============================================================
function run_hourly_step!(crop, params, doy, hour, Δt, lat, crop_name,
                          # Forcing (from interpolate_time output)
                          tmp_K, wnd, prc, rsd, shm, prs,
                          co2_ppm, water_stress,
                          layer_water, five_day_buffer, five_day_count, buffer_len,
                          mgmt::ManagementParams)
    # ----- 5-day temperature buffer -----
    if five_day_count < buffer_len
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
        lat=lat,
        doy=doy,
        hour=hour,
        crop_name=crop_name,
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
        crop_name=crop_name
    )

    gpp = phsyn_result.gpp
    rsp = phsyn_result.rsp
    tsp = phsyn_result.tsp

    # ----- CROP (growth simulation) -----
    crop_step!(crop;
        doy=doy, hour=hour, Δt=Δt,
        temperature=tmp_K, gpp=gpp, rsp=rsp,
        planting_doy=mgmt.planting_doy + Int(params.planting_offset),
        thermal_time_requirement=mgmt.thermal_time_requirement,
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
        leaf_weight_min=params.leaf_weight_min,
        leaf_weight_max=params.leaf_weight_max,
        leaf_weight_decay_rate=params.leaf_weight_decay_rate,
        max_crop_height=mgmt.max_crop_height,
        root_growth_rate=params.root_growth_rate,
        max_root_length=params.max_root_length,
        n_fertilizer=mgmt.n_fertilizer,
        co2_ppm=co2_ppm,
        is_irrigated=mgmt.is_irrigated,
        cold_damage_threshold=params.cold_damage_threshold,
        heat_damage_threshold=params.heat_damage_threshold,
        harvest_index=params.harvest_index,
        harvest_temp_threshold=params.harvest_temp_threshold,
        five_day_temp_buffer=five_day_buffer,
        five_day_temp_count=five_day_count,
        crop_name=crop_name
    )

    # ----- SOIL (water balance) -----
    tsp_W = tsp * L_vaporization
    soil_result = calc_soil_water(;
        layer_water=layer_water,
        transpiration=tsp_W,
        W2SF=prc,
        depth_root=crop.root_length,
        is_irrigated=mgmt.is_irrigated,
        Δt=Δt,
        soil_type_i=mgmt.soil_type,
        temperature=tmp_K,
        pressure=prs,
        wind_speed=wnd,
        specific_humidity=shm,
        crop_height=crop.crop_height,
        is_planted=crop.is_planted,
        crop_name=crop_name
    )

    return soil_result.water_stress, five_day_count
end

# ============================================================
# run_point_simulation — point (CSV) simulation
# ============================================================
function run_point_simulation(config::Config)
    # Read crop parameters
    params = read_crop_params(config.crop_param_file)

    # Load CSV forcing
    forcing_data = read_forcing_csv(config.csv_path)

    # Initialize simulation state
    crop, layer_water, water_stress, five_day_buffer, five_day_count =
        init_simulation_state(config.time_step)

    # Build management params from config
    mgmt = ManagementParams(
        config.planting_doy,
        config.is_irrigated,
        config.soil_type,
        config.n_fertilizer,
        config.thermal_time_requirement,
        params.height_coeff_a1,
        get(config.raster_vars, "wnd", Dict("height"=>10.0))["height"]
    )

    # Output storage
    daily_records = NamedTuple[]
    yearly_results = NamedTuple[]

    # Year loop
    for year in config.start_year:config.end_year
        println("\n  Year $year ...")
        co2_ppm = read_co2(config, year)
        year_forcing = forcing_data

        # Determine DOY range
        if config.start_year == config.end_year
            start_day = config.start_doy; end_day = config.end_doy
        elseif year == config.start_year
            start_day = config.start_doy; end_day = 365
        elseif year == config.end_year
            start_day = 1; end_day = config.end_doy
        else
            start_day = 1; end_day = 365
        end

        # DOY loop
        for doy in start_day:end_day
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
            Δt = config.time_step

            for ihour in 1:n_steps
                hour = (Float64(ihour) - 0.5) * Float64(Δt) / 3600.0

                # ----- Time interpolation (daily → hourly) -----
                hourly = interpolate_time(;
                    doy=doy, prev_doy=(doy == 1 ? 365 : doy - 1),
                    next_doy=(doy == 365 ? 1 : doy + 1), hour=hour,
                    lat=config.latitude, Δt=Δt,
                    tmax_prev=f_prev.tmax, tmax=f_today.tmax, tmax_next=f_next.tmax,
                    tmin_prev=f_prev.tmin, tmin=f_today.tmin, tmin_next=f_next.tmin,
                    radiation=f_today.radiation, precip=f_today.precip,
                    humidity=f_today.humidity, wind=f_today.wind,
                    pressure=f_today.pressure, ozone=f_today.ozone,
                    int_sinb=int_sinb,
                    wind_height=mgmt.wind_height
                )

                tmp_K = hourly.temperature
                wnd = max(hourly.wind, 0.001)
                prc = hourly.precipitation
                rsd = hourly.radiation
                shm = hourly.humidity
                prs = hourly.pressure

                buffer_len = 86400 ÷ Δt * 5
                water_stress, five_day_count = run_hourly_step!(
                    crop, params, doy, hour, Δt, config.latitude, config.crop_name,
                    tmp_K, wnd, prc, rsd, shm, prs,
                    co2_ppm, water_stress,
                    layer_water, five_day_buffer, five_day_count, buffer_len,
                    mgmt
                )
            end # hourly loop

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

    # Write output
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
    println("\n  Spatial mode: $n_actual_threads thread(s) running")

    # 1. Read crop parameters
    params = read_crop_params(config.crop_param_file)

    # 2. Get grid info
    lats, lons = get_grid_info(config)
    n_lon = length(lons)
    n_lat = length(lats)
    n_pixels = n_lon * n_lat
    println("    Grid: $n_lon x $n_lat = $n_pixels pixels")

    # 3. Create boundary mask (if boundary file is specified)
    boundary_mask = nothing
    if !isempty(config.boundary_file)
        if isfile(config.boundary_file)
            println("    Loading boundary: ", config.boundary_file)
            boundary_mask = create_boundary_mask(lons, lats, config.boundary_file; buffer_deg=config.boundary_buffer)
            n_in_boundary = sum(boundary_mask)
            println("    Boundary mask: $n_in_boundary / $n_pixels pixels within boundary")
        else
            println("    [WARN] Boundary file not found: $(config.boundary_file), running all pixels")
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
        planting_doy_field = load_management_param(config, "planting_doy", year, n_lon, n_lat; lats=lats, lons=lons)
        is_irrigated_field = load_management_param(config, "is_irrigated", year, n_lon, n_lat; lats=lats, lons=lons)
        soil_type_field = load_management_param(config, "soil_type", year, n_lon, n_lat; lats=lats, lons=lons)
        n_fertilizer_field = load_management_param(config, "n_fertilizer", year, n_lon, n_lat; lats=lats, lons=lons)
        thermal_time_field = load_management_param(config, "thermal_time_requirement", year, n_lon, n_lat; lats=lats, lons=lons)

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

            # Build per-pixel management params
            pixel_mgmt = ManagementParams(
                Int(planting_doy_field[i_lon, i_lat]),
                Int(is_irrigated_field[i_lon, i_lat]),
                Int(soil_type_field[i_lon, i_lat]),
                Float64(n_fertilizer_field[i_lon, i_lat]),
                Float64(thermal_time_field[i_lon, i_lat]),
                params.height_coeff_a1, 
                get(config.raster_vars, "wnd", Dict("height"=>10.0))["height"]  # wind_height from config
            )

            try
                results[index] = run_pixel_spatial(
                    config, params, year, lat, lon,
                    spatial_forcing, n_days, co2_ppm,
                    pixel_mgmt;
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
            @printf("    Average Yield=%.2f kg/ha, Average Harvest DOY=%.1f\n",
                    sum(v_yield)/length(v_yield), sum(v_harvest)/length(v_harvest))
            @printf("    Average Aboveground Biomass=%.2f kg/ha, Average Max LAI=%.2f\n",
                    sum(v_biomass)/length(v_biomass),sum(v_lai)/length(v_lai))
        end
    end

    # 5. Write output
    mkpath(config.output_dir)
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
                           mgmt::ManagementParams;
                           i_lon::Int=1, i_lat::Int=1)
    # Initialize simulation state
    Δt = config.time_step
    crop, layer_water, water_stress, five_day_buffer, five_day_count = init_simulation_state(Δt)
    buffer_len = 86400 ÷ Δt * 5

    for doy in config.start_doy:min(config.end_doy, n_days)
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

            # ----- Time interpolation -----
            # Spatial mode: use same-day tmax/tmin for prev/next (no adjacent-day data)
            hourly = interpolate_time(;
                doy=doy, prev_doy=(doy == 1 ? n_days : doy - 1),
                next_doy=(doy == n_days ? 1 : doy + 1),
                hour=hour, lat=lat, Δt=Δt,
                tmax_prev=weather_today.tmax, tmax=weather_today.tmax, tmax_next=weather_today.tmax,
                tmin_prev=weather_today.tmin, tmin=weather_today.tmin, tmin_next=weather_today.tmin,
                radiation=weather_today.radiation, precip=weather_today.precip,
                humidity=weather_today.humidity, wind=weather_today.wind,
                pressure=weather_today.pressure, ozone=weather_today.ozone,
                int_sinb=int_sinb, wind_height=mgmt.wind_height,
            )

            tmp_K = hourly.temperature
            wnd = max(hourly.wind, 0.001)
            prc = hourly.precipitation
            rsd = hourly.radiation
            shm = hourly.humidity
            prs = hourly.pressure

            water_stress, five_day_count = run_hourly_step!(
                crop, params, doy, hour, Δt, lat, config.crop_name,
                tmp_K, wnd, prc, rsd, shm, prs,
                co2_ppm, water_stress,
                layer_water, five_day_buffer, five_day_count, buffer_len,
                mgmt
            )
        end # hourly loop
    end # DOY loop

    return (yield=crop.yield, harvest_doy=crop.harvest_doy,
            LAI_max=crop.LAI_max_season, biomass_aboveground=crop.shoot_biomass_at_harvest)
end