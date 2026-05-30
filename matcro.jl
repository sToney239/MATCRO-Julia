# MATCRO - Unified Julia Crop Model
# Supports both CSV and NetCDF input formats
# Usage: julia matcro.jl config.toml

# ============================================================
# Module includes (in dependency order)
# Note: 00_routine.jl must be included LAST as it depends on all types
# ============================================================
include("lib/01_constants.jl")
include("lib/02_io.jl")
include("lib/02_1_point_io.jl")
include("lib/02_2_spatial_io.jl")
include("lib/03_hour_interpolation.jl")
include("lib/04_radiation.jl")
include("lib/05_1_photosynthesis_C3.jl")
include("lib/05_2_photosynthesis_C4.jl")
include("lib/05_photosynthesis.jl")
include("lib/06_crop.jl")
include("lib/07_soil.jl")
include("lib/00_routine.jl")

# ============================================================
# run_simulation — main entry point (dispatcher)
# ============================================================
function run_simulation(config_path::String)
    println("=" ^ 60)
    println("                   MATCRO (Julia Version)")
    println("=" ^ 60)
    flush(stdout)

    # Read configuration
    config = read_config(config_path)
    println("  Config:")
    println("    Crop: ", config.crop_name)
    println("    Period: $(config.start_year)-$(config.end_year)")
    println("    Input format: ", config.input_format)
    println("    Parameters loaded from: ", config.crop_param_file)

    # Dispatch to the appropriate simulation mode
    if config.input_format == "point"
        println("    Loading CSV weather data: ", config.csv_path)
        return run_point_simulation(config)
    elseif config.input_format == "raster"
        return run_spatial_simulation(config)
    else
        error("Unknown input format: $(config.input_format)")
    end
end

# Command-line entry point (skip when included from other scripts)
const _RUN_AS_MAIN = length(ARGS) > 0 && endswith(lowercase(ARGS[1]), ".toml")
if _RUN_AS_MAIN
    run_simulation(ARGS[1])
end
