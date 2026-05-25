# IO - Core data structures, TOML config/parameter reading, CO2 reading
# Companion files: io_csv.jl (CSV I/O), io_netcdf.jl (NetCDF I/O)

using TOML


# ============================================================
# Config — simulation settings (replaces Fortran SETTING.txt)
# ============================================================
@kwdef mutable struct Config
    # Time
    start_year::Int
    end_year::Int
    start_doy::Int
    end_doy::Int
    time_step::Int             # Δt [seconds] (Fortran: TRES)

    # Location (single point)
    latitude::Float64          # [degree]

    # Crop
    crop_name::String          # "Maize", "Rice", "Wheat", "Soybeans"
    crop_param_file::String    # path to TOML parameter file
    planting_doy::Int          # planting day of year
    is_irrigated::Int          # 0 = rainfed, 1 = irrigated

    # Soil & management
    soil_type::Int             # soil texture index (1-13)
    n_fertilizer::Float64      # nitrogen fertilizer [kg N/ha]
    thermal_time_requirement::Float64  # GDH at maturity

    # CO2
    co2_file::String           # path to CO2 file (empty → use fixed_ppm)
    co2_fixed_ppm::Float64     # fallback CO2 [ppm]

    # Input
    input_format::String       # "csv" or "netcdf" (auto-detected)
    csv_path::String           # CSV forcing file path (if format=csv)

    # NetCDF input metadata (populated if format=netcdf)
    netcdf_vars::Dict{String,Dict{String,Any}}  # variable → metadata dict
    nthreads::Int              # parallel threads for spatial NetCDF
    lon_dim::String            # longitude dimension name in NetCDF
    lat_dim::String            # latitude dimension name in NetCDF
    time_dim::String           # time dimension name in NetCDF

    # Output
    output_dir::String
    output_format::String      # "csv" or "netcdf"
end

# ============================================================
# CropParameters — all crop parameters (replaces Fortran RDPRM)
# Key names match Julia crop_step! / rad / leaf_photosynthesis_c4 args
# ============================================================
@kwdef mutable struct CropParameters
    # ----- Photosynthesis (→ leaf_photosynthesis_c4) -----
    respiration_coeff::Float64 = 0.015         # Fortran: RESPCP
    quantum_efficiency::Float64 = 0.05         # Fortran: EFFCON
    a_theta::Float64 = 0.8                     # Fortran: ATHETA
    b_theta::Float64 = 0.95                    # Fortran: BTHETA
    m_H2O::Float64 = 4.0                       # Fortran: MH2O
    b_H2O::Float64 = 0.04                      # Fortran: BH2O

    # ----- Radiation (→ calc_radiation) -----
    k_nitrogen::Float64 = 0.3                  # Fortran: KN
    leaf_PAR_reflectance::Float64 = 0.105      # Fortran: RLFV
    leaf_PAR_transmittance::Float64 = 0.07     # Fortran: TLFV
    leaf_NIR_reflectance::Float64 = 0.58       # Fortran: RLFN
    leaf_NIR_transmittance::Float64 = 0.25     # Fortran: TLFN

    # ----- Development (→ crop_step!) -----
    half_progress::Float64 = 0.52              # Fortran: hDVS
    needs_vernalization::Int = 0               # Fortran: VN
    base_temp::Float64 = 8.0                   # Fortran: TB
    optimal_temp::Float64 = 30.0               # Fortran: TO
    ceiling_temp::Float64 = 42.0               # Fortran: TH
    vernalization_saturation::Float64 = 60.0

    # ----- Biomass conversion (→ crop_step!) -----
    k_leaf_convert::Float64 = 0.871            # Fortran: CFLF
    k_stem_convert::Float64 = 0.810            # Fortran: CFST
    k_root_convert::Float64 = 0.857            # Fortran: CFRT
    k_grain_convert::Float64 = 0.815           # Fortran: CFSO

    # ----- Partitioning: Root (→ crop_step!) -----
    shoot_progress_1::Float64 = 0.35           # Fortran: RTX
    shoot_alloc_ratio_1::Float64 = 0.25        # Fortran: RTY
    shoot_progress_2::Float64 = 0.72           # Fortran: RTX2

    # ----- Partitioning: Leaf (→ crop_step!) -----
    leaf_alloc_ratio_0::Float64 = 0.49         # Fortran: LEVY0
    leaf_progress_1::Float64 = 0.25            # Fortran: LEVX1
    leaf_alloc_ratio_1::Float64 = 0.49         # Fortran: LEVY1
    leaf_progress_2::Float64 = 0.48            # Fortran: LEVX2
    leaf_alloc_ratio_2::Float64 = 0.0          # Fortran: LEVY2

    # ----- Partitioning: Panicle (→ crop_step!) -----
    panicle_progress_1::Float64 = 0.37         # Fortran: PNCLX1
    panicle_alloc_ratio_1::Float64 = 0.0       # Fortran: PNCLY1
    panicle_progress_2::Float64 = 0.6          # Fortran: PNCLX2
    panicle_alloc_ratio_2::Float64 = 1.0       # Fortran: PNCLY2
    panicle_progress_3::Float64 = 1.0          # Fortran: PNCLX3
    panicle_alloc_ratio_3::Float64 = 1.0       # Fortran: PNCLY3

    # ----- Partitioning: Dead leaf (→ crop_step!) -----
    dead_progress_1::Float64 = 0.0             # Fortran: DLVX1
    dead_ratio_1::Float64 = 0.0                # Fortran: DLVY1
    dead_progress_2::Float64 = 0.65            # Fortran: DLVX2
    dead_ratio_2::Float64 = 0.0                # Fortran: DLVY2
    dead_progress_3::Float64 = 1.0             # Fortran: DLVX3
    dead_ratio_3::Float64 = 0.0000003          # Fortran: DLVY3

    # ----- Stem reserve (→ crop_step!) -----
    fraction_starch_reserve::Float64 = 0.35    # Fortran: FSTR

    # ----- Specific leaf weight (→ crop_step!) -----
    leaf_weight_asymptote::Float64 = 400.0     # Fortran: SLWYA
    leaf_weight_intercept::Float64 = 700.0     # Fortran: SLWYB
    leaf_weight_decay_rate::Float64 = 3.0      # Fortran: SLWX

    # ----- Plant height -----
    height_coeff_a1::Float64 = 2.0              # Fortran: HGTAA
    height_coeff_a2::Float64 = 0.6753014        # Fortran: HGTAB
    height_coeff_b1::Float64 = 0.3664009        # Fortran: HGTBA
    height_coeff_b2::Float64 = 0.3175134        # Fortran: HGTBB

    # ----- Root growth (→ crop_step!) -----
    root_growth_rate::Float64 = 0.06           # Fortran: GZRT [m/day]
    max_root_length::Float64 = 1.5             # Fortran: MXRT [m]

    # ----- Soil water stress -----
    soil_stress_param::Float64 = 0.003297       # Fortran: GMMSL

    # ----- Specific leaf nitrogen (→ crop_step!) -----
    leaf_nitrogen_x1::Float64 = 0.0            # Fortran: SLNX1
    leaf_nitrogen_x2::Float64 = 0.52           # Fortran: SLNX2
    leaf_nitrogen_x3::Float64 = 1.0            # Fortran: SLNX3
    leaf_nitrogen_max::Float64 = 1.675         # Fortran: SLNYMX
    leaf_nitrogen_min::Float64 = 0.825         # Fortran: SLNYMN
    leaf_nitrogen_sensitivity::Float64 = 0.004 # Fortran: SLNK

    # ----- Temperature damage (→ crop_step!) -----
    cold_damage_threshold::Float64 = -100.0    # Fortran: TCmin
    heat_damage_threshold::Float64 = 100.0     # Fortran: THcrit

    # ----- Harvest (→ crop_step!) -----
    harvest_index::Float64 = 0.83              # Fortran: HI
    harvest_temp_threshold::Float64 = 8.0      # Fortran: HVT_TAVE

    # ----- Leaf loss (→ crop_step!) -----
    LAI_threshold_grain::Float64 = 4.0         # Fortran: LLFst
    k_leaf_loss::Float64 = 1.0                 # Fortran: kLLF

    # ----- Planting offset -----
    planting_offset::Float64 = 0.0              # Fortran: PLTDIF

    # ----- LTCH (bulk transfer, used in PHSYN) -----
    bulk_transfer_coeff::Float64 = 0.06         # Fortran: LTCH

    # ----- Rubisco kinetic params (for C3 photosynthesis) -----
    ZKCA::Float64 = 404.9
    ZKCB::Float64 = 79430.0
    ZKOA::Float64 = 278.4
    ZKOB::Float64 = 36380.0
    GMMA::Float64 = 42.8
    GMMB::Float64 = 37830.0
end

# ============================================================
# DailyForcing — one day of meteorological data at one point
# All temperatures in [K] (matching forcing file)
# ============================================================
@kwdef mutable struct DailyForcing
    doy::Int
    tmax::Float64              # daily maximum temperature [K]
    tmin::Float64              # daily minimum temperature [K]
    radiation::Float64         # downward shortwave radiation [W/m²]
    precip::Float64            # precipitation [kg/m²/s]
    humidity::Float64          # specific humidity [kg/kg]
    wind::Float64              # wind speed [m/s]
    pressure::Float64          # surface pressure [Pa]
    ozone::Float64 = 0.0       # ozone concentration
end

# ============================================================
# read_config — parse TOML config file → Config
# ============================================================
function read_config(config_path::String)::Config
    toml = TOML.parsefile(config_path)

    # Get config directory for relative path resolution
    config_dir = dirname(config_path)

    t = toml["time"]
    loc = toml["location"]
    crp = toml["crop"]
    sl = toml["soil"]
    co2 = toml["co2"]
    inp = toml["input"]
    out = toml["output"]

    # Auto-detect input format: [input.csv] or [input.netcdf]
    nc_vars = Dict{String,Dict{String,Any}}()
    nthreads = 1
    lon_dim = "lon"
    lat_dim = "lat"
    time_dim = "time"

    if haskey(inp, "csv") && haskey(inp["csv"], "path")
        input_format = "csv"
        csv_path = isabspath(inp["csv"]["path"]) ? inp["csv"]["path"] : joinpath(config_dir, inp["csv"]["path"])
    elseif haskey(inp, "netcdf")
        input_format = "netcdf"
        csv_path = ""

        # NetCDF settings
        nc_cfg = inp["netcdf"]
        nthreads = get(nc_cfg, "nthreads", 1)
        lon_dim = get(nc_cfg, "lon_dim", "lon")
        lat_dim = get(nc_cfg, "lat_dim", "lat")
        time_dim = get(nc_cfg, "time_dim", "time")

        # Per-variable metadata (each var has its own section)
        for (varname, meta) in nc_cfg
            if isa(meta, Dict)
                # Resolve file path relative to config_dir
                if haskey(meta, "file") && haskey(meta, "variable")
                    file_path = isabspath(meta["file"]) ? meta["file"] : joinpath(config_dir, meta["file"])
                    nc_vars[varname] = Dict{String,Any}(
                        "file" => file_path,
                        "variable" => meta["variable"],
                        "height" => get(meta, "height", 10.0),
                        "scale_factor" => get(meta, "scale_factor", 1.0),
                        "add_offset" => get(meta, "add_offset", 0.0),
                        "default_value" => get(meta, "default_value", nothing),
                    )
                elseif haskey(meta, "default_value")
                    # Management param with default_value only (no NC file)
                    nc_vars[varname] = Dict{String,Any}(
                        "default_value" => meta["default_value"],
                    )
                end
            end
        end

        # Ensure all 5 management params have entries in nc_vars
        mgmt_defaults = Dict{String,Any}(
            "planting_doy" => 120,
            "is_irrigated" => 0,
            "soil_type" => 9,
            "n_fertilizer" => 100.0,
            "thermal_time_requirement" => 1500.0,
        )
        for (pname, pdefault) in mgmt_defaults
            if !haskey(nc_vars, pname)
                nc_vars[pname] = Dict{String,Any}("default_value" => pdefault)
            elseif nc_vars[pname]["default_value"] === nothing
                nc_vars[pname]["default_value"] = pdefault
            end
        end
    else
        error("Must specify either [input.csv] with path or [input.netcdf] with base_dir")
    end

    # Resolve relative paths relative to config directory
    crop_param_file = isabspath(crp["param_file"]) ? crp["param_file"] : joinpath(config_dir, crp["param_file"])
    co2_file_val = get(co2, "file", "")
    if co2_file_val != "" && !isabspath(co2_file_val)
        co2_file_val = joinpath(config_dir, co2_file_val)
    end
    output_dir = isabspath(out["directory"]) ? out["directory"] : joinpath(config_dir, out["directory"])

    return Config(;
        start_year        = t["start_year"],
        end_year          = t["end_year"],
        start_doy         = t["start_doy"],
        end_doy           = t["end_doy"],
        time_step         = t["time_step"],
        latitude          = loc["latitude"],
        crop_name         = crp["crop_name"],
        crop_param_file   = crop_param_file,
        planting_doy      = crp["planting_doy"],
        is_irrigated      = crp["is_irrigated"],
        soil_type         = sl["soil_type"],
        n_fertilizer      = sl["n_fertilizer"],
        thermal_time_requirement = sl["thermal_time_requirement"],
        co2_file          = co2_file_val,
        co2_fixed_ppm     = get(co2, "fixed_ppm", 400.0),
        input_format      = input_format,
        csv_path          = csv_path,
        netcdf_vars       = nc_vars,
        nthreads          = nthreads,
        lon_dim           = lon_dim,
        lat_dim           = lat_dim,
        time_dim          = time_dim,
        output_dir        = output_dir,
        output_format     = out["format"],
    )
end

# ============================================================
# read_crop_params — parse TOML crop parameter file → CropParameters
# ============================================================
function read_crop_params(toml_path::String)::CropParameters
    toml = TOML.parsefile(toml_path)
    crop_params = toml["crop_parameters"]

    CropParameters(;
        respiration_coeff         = crop_params["respiration_coeff"],
        quantum_efficiency        = crop_params["quantum_efficiency"],
        a_theta                   = crop_params["a_theta"],
        b_theta                   = crop_params["b_theta"],
        m_H2O                     = crop_params["m_H2O"],
        b_H2O                     = crop_params["b_H2O"],
        k_nitrogen                = crop_params["k_nitrogen"],
        leaf_PAR_reflectance      = crop_params["leaf_PAR_reflectance"],
        leaf_PAR_transmittance    = crop_params["leaf_PAR_transmittance"],
        leaf_NIR_reflectance      = crop_params["leaf_NIR_reflectance"],
        leaf_NIR_transmittance    = crop_params["leaf_NIR_transmittance"],
        half_progress            = crop_params["half_progress"],
        needs_vernalization      = crop_params["needs_vernalization"],
        base_temp                = crop_params["base_temp"],
        optimal_temp             = crop_params["optimal_temp"],
        ceiling_temp             = crop_params["ceiling_temp"],
        vernalization_saturation = crop_params["vernalization_saturation"],
        k_leaf_convert           = crop_params["k_leaf_convert"],
        k_stem_convert           = crop_params["k_stem_convert"],
        k_root_convert           = crop_params["k_root_convert"],
        k_grain_convert          = crop_params["k_grain_convert"],
        shoot_progress_1         = crop_params["shoot_progress_1"],
        shoot_alloc_ratio_1      = crop_params["shoot_alloc_ratio_1"],
        shoot_progress_2         = crop_params["shoot_progress_2"],
        leaf_alloc_ratio_0       = crop_params["leaf_alloc_ratio_0"],
        leaf_progress_1          = crop_params["leaf_progress_1"],
        leaf_alloc_ratio_1       = crop_params["leaf_alloc_ratio_1"],
        leaf_progress_2          = crop_params["leaf_progress_2"],
        leaf_alloc_ratio_2       = crop_params["leaf_alloc_ratio_2"],
        panicle_progress_1       = crop_params["panicle_progress_1"],
        panicle_alloc_ratio_1    = crop_params["panicle_alloc_ratio_1"],
        panicle_progress_2       = crop_params["panicle_progress_2"],
        panicle_alloc_ratio_2    = crop_params["panicle_alloc_ratio_2"],
        panicle_progress_3       = crop_params["panicle_progress_3"],
        panicle_alloc_ratio_3    = crop_params["panicle_alloc_ratio_3"],
        dead_progress_1          = crop_params["dead_progress_1"],
        dead_ratio_1             = crop_params["dead_ratio_1"],
        dead_progress_2          = crop_params["dead_progress_2"],
        dead_ratio_2             = crop_params["dead_ratio_2"],
        dead_progress_3          = crop_params["dead_progress_3"],
        dead_ratio_3             = crop_params["dead_ratio_3"],
        fraction_starch_reserve  = crop_params["fraction_starch_reserve"],
        leaf_weight_asymptote    = crop_params["leaf_weight_asymptote"],
        leaf_weight_intercept    = crop_params["leaf_weight_intercept"],
        leaf_weight_decay_rate   = crop_params["leaf_weight_decay_rate"],
        height_coeff_a1          = crop_params["height_coeff_a1"],
        height_coeff_a2          = crop_params["height_coeff_a2"],
        height_coeff_b1          = crop_params["height_coeff_b1"],
        height_coeff_b2          = crop_params["height_coeff_b2"],
        root_growth_rate         = crop_params["root_growth_rate"],
        max_root_length          = crop_params["max_root_length"],
        soil_stress_param        = crop_params["soil_stress_param"],
        leaf_nitrogen_x1         = crop_params["leaf_nitrogen_x1"],
        leaf_nitrogen_x2         = crop_params["leaf_nitrogen_x2"],
        leaf_nitrogen_x3         = crop_params["leaf_nitrogen_x3"],
        leaf_nitrogen_max        = crop_params["leaf_nitrogen_max"],
        leaf_nitrogen_min        = crop_params["leaf_nitrogen_min"],
        leaf_nitrogen_sensitivity = crop_params["leaf_nitrogen_sensitivity"],
        cold_damage_threshold    = crop_params["cold_damage_threshold"],
        heat_damage_threshold    = crop_params["heat_damage_threshold"],
        harvest_index            = crop_params["harvest_index"],
        harvest_temp_threshold   = crop_params["harvest_temp_threshold"],
        LAI_threshold_grain      = crop_params["lai_threshold_grain"],
        k_leaf_loss              = crop_params["k_leaf_loss"],
        planting_offset          = crop_params["planting_offset"],
        bulk_transfer_coeff      = get(crop_params, "bulk_transfer_coeff", 0.06),
        ZKCA                     = get(crop_params, "ZKCA", 404.9),
        ZKCB                     = get(crop_params, "ZKCB", 79430.0),
        ZKOA                     = get(crop_params, "ZKOA", 278.4),
        ZKOB                     = get(crop_params, "ZKOB", 36380.0),
        GMMA                     = get(crop_params, "GMMA", 42.8),
        GMMB                     = get(crop_params, "GMMB", 37830.0),
    )
end

# ============================================================
# read_co2 — read CO2 concentration for a given year
# File format: CSV with columns year,co2_ppm
# If file is not provided or year not found, returns fixed_ppm
# ============================================================
function read_co2(config::Config, year::Int)::Float64
    if isempty(config.co2_file) || !isfile(config.co2_file)
        return config.co2_fixed_ppm
    end

    # Try to read as CSV
    try
        lines = readlines(config.co2_file)
        for line in lines
            line = strip(line)
            isempty(line) && continue
            startswith(line, '#') && continue  # skip comments

            parts = split(line, ',')
            if length(parts) >= 2
                y = tryparse(Int, strip(parts[1]))
                if y !== nothing && y == year
                    return parse(Float64, strip(parts[2]))
                end
            end
        end
    catch
        # If CSV parsing fails, fall back to fixed_ppm
    end

    # Year not found or file parse failed → use fixed
    return config.co2_fixed_ppm
end

# ============================================================
# calc_int_sinb — integrated sin(solar elevation) over a day
# Used by tinterp to distribute daily radiation to hourly
# ============================================================
function calc_int_sinb(doy::Int, lat::Float64, Δt::Int)::Float64
    int_sinb = 0.0
    n_steps = 86400 ÷ Δt
    for ihour in 1:n_steps
        hour = (Float64(ihour) - 0.5) * Float64(Δt) / 3600.0
        int_sinb += sin_solar_elevation(doy, hour, lat) * Float64(Δt)
    end
    return int_sinb
end
