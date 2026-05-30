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

### `[general]`

Core simulation settings that apply to both point and raster modes:

| Parameter | Description |
|-----------|-------------|
| `crop_name` | Crop type: `"Maize"`, `"Rice"`, `"Wheat"`, `"Soybeans"` |
| `crop_param` | Path to crop parameter TOML file |
| `start_year` / `end_year` | Simulation year range |
| `start_doy` / `end_doy` | Simulation day-of-year range |
| `time_step` | Time step in seconds (typically use 1 hour, which is 3600s) |
| `co2_ppm_default` | Default CO₂ concentration [ppm] |
| `co2_yearly_file` | Path to yearly CO₂ CSV file (optional) |

### Input: choose `[point_simulation]` **or** `[spatial_simulation]`

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

## Input Data Format

### 1. Point Data Format (CSV)

When running point simulation, the weather data should be CSV format (specified by `csv_path` in `[point_simulation.weather]`). Your CSV file must have the following columns:

| Column | Unit | Description |
|--------|------|-------------|
| year | — | Year |
| doy | — | Day of year (1–366) |
| tmax | K | Daily maximum temperature |
| tmin | K | Daily minimum temperature |
| radiation | W/m² | Downward shortwave radiation |
| precip | kg/m²/s | Precipitation rate |
| humidity | kg/kg | Specific humidity |
| wind | m/s | Wind speed |
| pressure | Pa | Surface air pressure |
| ozone | — | Ozone concentration (can be 0) |

The file must include a header row with these column names.

Management parameters are specified in `[point_simulation.management]`.

### 2. Raster Input Data Format (NetCDF & TIFF)

When using NetCDF or TIFF input (`[spatial_simulation]` in config.toml), your data files must follow these format requirements.


#### 2.1 Weather variables (NetCDF format required)

##### 2.1.1 Spatial dimensions

Dimension names are configurable via `lon_dim`, `lat_dim`, `time_dim` in `[spatial_simulation.weather]` (defaults: `lon`, `lat`, `time`). Common alternatives like `latitude`, `longitude` are also auto-detected as fallbacks.

##### 2.1.2 Time dimension

The time variable in your NetCDF file must follow one of these formats:

1. **Stored as dates**: time values are directly stored as calendar dates (e.g., `2021-01-01`, `2021-01-02`, ...)
2. **Stored as numbers with a `units` attribute**: numeric values (e.g., 0, 1, 2, ...) with a `units` attribute following the format `"days since YYYY-MM-DD"` (e.g., `"days since 2021-01-01"`). The program converts each number to a date using the reference date in `units`, then slices by year.

> If the `units` attribute is missing or not in the `"days since YYYY-MM-DD"` format, the program will report an error.

##### 2.1.3 Weather Variable dimension

Each meteorological variable is specified in `[spatial_simulation.weather.<name>]` with user-friendly names:

| Config key | Variable | Description |
|------------|----------|-------------|
| `temperature_max` | tmx | Daily maximum temperature |
| `temperature_min` | tmn | Daily minimum temperature |
| `precipitation` | prc | Precipitation rate |
| `radiation` | rsd | Downward shortwave radiation |
| `humidity` | shm | Specific humidity |
| `wind_speed` | wnd | Wind speed |
| `pressure` | prs | Surface air pressure |

Each weather variable section requires:
- `file`: path to the NetCDF file (relative to config file location, or absolute)
- `variable`: variable name inside the file

Optional keys:
- `height`: wind measurement height (default: 10.0 m), only for `wind_speed`
- `scale_factor`: scaling factor for data values
- `add_offset`: offset added to data values

#### 2.2 Management parameters (NetCDF or TIF format)

> Belows describes only the TIF format requirement, if you choose NetCDF file, the format should be like in 2.2

Management parameters are specified in `[spatial_simulation.management.<param_name>]` sections:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `planting_doy` | Planting day of year | 120 |
| `is_irrigated` | 0 = rainfed, 1 = irrigated | 0 |
| `soil_type` | Soil texture index (1-13) | 9 |
| `n_fertilizer` | Nitrogen fertilizer [kg N/ha] | 100.0 |
| `thermal_time_requirement` | GDH at maturity | 1500.0 |

Each parameter section supports:
1. **Uniform default only**: set `default_value` — used for all pixels when no file is provided.
2. **Spatial NC file**: add `file` and `variable` keys alongside `default_value`. The file can be 2D (lon, lat) for static parameters, or 3D (lon, lat, time) with a year dimension for time-varying parameters.
3. **Spatial TIF file**: add `file` key with a `.tif`/`.tiff` path alongside `default_value`. The TIF file should be single-band (static, used for all years). Nearest-neighbor resampling is used when TIF and NC grid resolutions differ. A bbox check is performed to ensure the TIF extent covers the simulation grid.

Priority: file (NC or TIF) > `default_value` > built-in default.

#### 2.3 Boundary filtering (optional)

Specify a GeoJSON boundary file to filter pixels:

```toml
[spatial_simulation.boundary]
file = "data/boundaries/region.geojson"
buffer_deg = 0.01    # buffer distance (degrees) for boundary contact detection
```

>Later `shp` and `TIF` mask will be supported.

For full details, see [README-CONFIG.md](README-CONFIG.md).

### Spatial Output Format

When running raster spatial simulation with `format = "raster"` (or `"geotiff"`), the model outputs **one GeoTIFF file per year**. Four files are generated for each simulation year:

- `yield_YYYY.tif`: Crop yield (kg/ha), Float64
- `harvest_doy_YYYY.tif`: Harvest day of year, Int32
- `LAI_max_YYYY.tif`: Seasonal maximum leaf area index (m²/m²), Float64
- `biomass_aboveground_YYYY.tif`: Aboveground biomass at harvest (kg/ha), Float64

The TIF files use WGS84 (EPSG:4326) coordinate reference system, with geotransform metadata for proper geospatial alignment.