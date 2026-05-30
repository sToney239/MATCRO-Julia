# MATCRO-Julia

This is the Julia version of the MATCRO crop model, originally written in Fortran. For model details, please refer to:

- Masutomi, Y., Usui, K., Manaka, T., Nishimori, M., Shimizu, M., Takimoto, K., and Arai, A.: Development of a model (MATCRO) for simulating the effects of climate change on regional crop production, Geosci. Model Dev., 9, 4133–4150, https://doi.org/10.5194/gmd-9-4133-2016, 2016.
- Masutomi, Y.: Improvement of MATCRO model for climate change impact assessment on global crop production, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2025-1885, 2025.
- Masutomi, Y.: MATCRO-SOM: a new soil organic matter model considering priming effects for climate change impact assessment on global crop production, Geosci. Model Dev., 18, 8801–8826, https://doi.org/10.5194/gmd-18-8801-2025, 2025.

> **Warning:** The modules of Maize and Soybeans have been fully cross-validated against the original Fortran code with acceptable precision (error < 0.01%), but currently for rice and wheats, no test has been done.

## Features of this version

- **Readable variable names** replacing original Fortran ones, making the code easier to follow:

  | Fortran name | Julia name | Meaning |
  |---|---|---|
  | `WST` | `stem_biomass` | Stem biomass of crop |
  | `RLFV` | `leaf_PAR_reflectance` | Leaf PAR reflectance |
  | `RESPCP` | `respiration_coeff` | Leaf respiration coefficient |

- **Cross-platform**: Julia runs natively on Windows, macOS, and Linux, whereas Fortran tooling on Windows is limited.
- **Parallel computation**: Julia provides built-in multi-core parallelism, enabling large-ensemble or multi-site simulations.
- **Equation-referenced comments**: Code is annotated with the corresponding equation numbers from the source papers as much as possible, and combined with the readable variable names this makes it easier to understand.

## Julia Setup

### Install Julia

Download Julia from [julialang.org](https://julialang.org/downloads/) (1.10 or later). Follow the platform-specific instructions for Windows, macOS, or Linux.

### Install dependencies

All external packages are listed in `Project.toml`. Install them at once:

```julia
using Pkg
Pkg.activate(".")   # activate the project environment
Pkg.instantiate()   # install all dependencies from Project.toml
```

Or install individually:

| Package | Required for |
|---------|-------------|
| `CSV` | Point (CSV) mode |
| `DataFrames` | Point (CSV) mode |
| `NCDatasets` | Spatial (NetCDF) mode |
| `ArchGDAL` | Spatial (GeoTIFF output, TIF input) |
| `JSON` | Spatial (GeoJSON boundary) |

> If you only use point (CSV) input mode, `NCDatasets`, `ArchGDAL`, and `JSON` are not required. If you only use spatial mode, `CSV` and `DataFrames` are not required. However, using `Pkg.instantiate()` with the project environment will install all packages cleanly and is recommended.

## Quick start

```bash
julia matcro.jl <config_path>
```

This command consists of three parts separated by spaces: the `julia` interpreter, the main program script `matcro.jl`, and the path to a `.toml` configuration file.

- `<config_path>` can be any `.toml` file, located anywhere on your system — no need to place it in the project root
- All relative paths in the `.toml` file (e.g., `crop_param`, `file`) are resolved relative to the config file's own directory
- The model source code lives in the `lib/` directory, with `matcro.jl` as the main entry point. The `.toml` file is a configuration file that Julia can read — you are free to modify its contents. For detailed parameter descriptions, see the example configs in `example/` and the [README-CONFIG.md](README-CONFIG.md) reference.

For full configuration `toml` reference, see [README-CONFIG.md](README-CONFIG.md).

### Try it out

The project includes ready-to-run example datasets. Try the following commands:

```bash
# Point simulation with single year CSV data
julia matcro.jl example/csv/config.toml

# Point simulation with multi-year CSV data
julia matcro.jl example/csv_multi_year/config.toml

# Spatial simulation with NetCDF weather data and NetCDF management data
julia matcro.jl example/netcdf/config.toml

# Spatial simulation with NetCDF weather data and TIF management data
# with boundary filter
julia matcro.jl example/tif/config.toml
```

### Multi-threaded parallel simulation

**Raster spatial mode** supports multi-threaded parallel pixel simulation using Julia's built-in threading. 

> Please note that for point format only supports single thread simulation.

The thread count is controlled by the Julia runtime, **not** by the config file.

**How to enable multi-threading:**

```bash
# Use 4 threads
julia -t 4 matcro.jl example/tif/config.toml
```

You can also set the `JULIA_NUM_THREADS` environment variable instead of the `-t` flag.

**How many threads should you use?**

- Check available CPU threads:
  - Windows: `echo %NUMBER_OF_PROCESSORS%` or Task Manager → Performance → CPU → Logical processors
  - macOS: `sysctl -n hw.ncpu`
  - Linux: `nproc`
- **Recommendation**: Use physical threads minus 1 or 2, leaving for your system and other tasks. For example, on a 16-thread machine, use 14 threads.

## Configuration Structure

The config TOML file has two main sections:

### 1. `[general]`

Core simulation settings that apply to both point and raster modes:

```toml
[general]
crop_name = "Maize"        # Maize, Rice, Wheat, Soybean
crop_param = "params/maize.toml"
start_year = 2000
end_year = 2000
start_doy = 1
end_doy = 365
time_step = 3600           # how many seconds in 1 step, normal should be 1 hour [seconds]
co2_ppm_default = 400.0    # default CO2 concentration [ppm]
co2_yearly_file = "data/co2/co2.csv" 
```

### 2. `[point_simulation]` **or** `[spatial_simulation]`

#### Point mode — `[point_simulation]`

```toml
[point_simulation]
latitude = 40.0
output_directory = "output/"

[point_simulation.weather]
csv_path = "weather_data.csv"

[point_simulation.management]
planting_doy = 120
is_irrigated = 0
soil_type = 9
n_fertilizer = 100.0
thermal_time_requirement = 1500.0
```

#### Spatial mode — `[spatial_simulation]`

```toml
[spatial_simulation]
output_directory = "output/"

[spatial_simulation.weather]
lon_dim = "lon"
lat_dim = "lat"
time_dim = "time"

[spatial_simulation.weather.temperature_max]
file = "data/netcdf/tmx.nc4"
variable = "tmx"
# ... other weather variables ...

[spatial_simulation.management]
[spatial_simulation.management.planting_doy]
default_value = 120
# ... other management params (with optional spatial files) ...
```

For full details, see [README-CONFIG.md](README-CONFIG.md).
