# MATCRO-Julia

This is the Julia version of the MATCRO crop model, originally written in Fortran. For model details, please refer to:

- Masutomi, Y., Usui, K., Manaka, T., Nishimori, M., Shimizu, M., Takimoto, K., and Arai, A.: Development of a model (MATCRO) for simulating the effects of climate change on regional crop production, Geosci. Model Dev., 9, 4133–4150, https://doi.org/10.5194/gmd-9-4133-2016, 2016.
- Masutomi, Y.: Improvement of MATCRO model for climate change impact assessment on global crop production, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2025-1885, 2025.
- Masutomi, Y.: MATCRO-SOM: a new soil organic matter model considering priming effects for climate change impact assessment on global crop production, Geosci. Model Dev., 18, 8801–8826, https://doi.org/10.5194/gmd-18-8801-2025, 2025.

> **Warning:** Only the C4 crop (Maize) pathway has been cross-validated against the original Fortran code. The C3 photosynthesis module (`src/05_1_photosynthesis_C3.jl`) for Soybeans (as well as Rice) still requires verification. Please wait for future updates.

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

Download Julia from [julialang.org](https://julialang.org/downloads/) (1.9 or later). Follow the platform-specific instructions for Windows, macOS, or Linux.

### Install dependencies

`NCDatasets` is the only external package you need to install (Other packages required are Julia standard libraries):

```julia
using Pkg
Pkg.add("NCDatasets")
```

> If you only use CSV input mode, `NCDatasets` is not required.

## Quick start

```bash
julia src/00_matcro.jl <config_path>
```

This command consists of three parts separated by spaces: the `julia` interpreter, the main program script `src/00_matcro.jl`, and the path to a `.toml` configuration file.

- `<config_path>` can be any `.toml` file, located anywhere on your system — no need to place it in the project root
- All relative paths in the `.toml` file (e.g., `param_file`, `file`) are resolved relative to the config file's own directory
- The model source code lives in the `src/` directory, with `00_matcro.jl` as the main entry point. The `.toml` file is a configuration file that Julia can read — you are free to modify its contents. For detailed parameter descriptions, see the example configs in `example/` and the [README-CONFIG.md](README-CONFIG.md) reference.

For full configuration `toml` reference, see [README-CONFIG.md](README-CONFIG.md).

### Try it out

The project includes ready-to-run example datasets. Try the following commands:

```bash
# CSV single-year simulation
julia src/00_matcro.jl example/csv/config.toml

# CSV multi-year simulation
julia src/00_matcro.jl example/csv_multi_year/config.toml

# NetCDF spatial simulation
julia src/00_matcro.jl example/netcdf/config.toml
```

### Multi-threaded parallel simulation

NetCDF spatial mode supports multi-threaded parallel pixel simulation. Set the thread count via `nthreads` in the `[input.netcdf]` section of your config file.

**How many threads should you use?**
- Check available CPU cores:
  - Windows: `echo %NUMBER_OF_PROCESSORS%` or Task Manager → Performance → CPU → Logical processors
  - macOS: `sysctl -n hw.ncpu`
  - Linux: `nproc`
- **Recommendation**: Set `nthreads` to physical cores minus 1 or 2, leaving headroom for your system and other tasks. For example, on a 16-core machine, use 14 threads.

## Input Data Format

### CSV Input Data Format

When using CSV input (`[input.csv]` in config.toml), your data file must have the following columns:

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

### NetCDF Input Data Format

When using NetCDF input (`[input.netcdf]` in config.toml), your data files must follow these format requirements:

#### Spatial dimensions

Dimension names are configurable via `lon_dim`, `lat_dim`, `time_dim` in `[input.netcdf]` (defaults: `lon`, `lat`, `time`). Common alternatives like `latitude`, `longitude` are also auto-detected as fallbacks.

#### Time dimension

The time variable in your NetCDF file must follow one of these formats:

1. **Stored as dates**: time values are directly stored as calendar dates (e.g., `2021-01-01`, `2021-01-02`, ...)
2. **Stored as numbers with a `units` attribute**: numeric values (e.g., 0, 1, 2, ...) with a `units` attribute following the format `"days since YYYY-MM-DD"` (e.g., `"days since 2021-01-01"`). The program converts each number to a date using the reference date in `units`, then slices by year.

> If the `units` attribute is missing or not in the `"days since YYYY-MM-DD"` format, the program will report an error.

#### Variable structure

Each meteorological variable (tmx, tmn, prc, rsd, shm, wnd, prs) is specified individually in `[input.netcdf.<var>]` with:

- `file`: path to the NetCDF file (relative to config file location, or absolute)
- `variable`: variable name inside the file 
  - you should check the dimension name of the NetCDF file.
  - It may differ from the file name, e.g., key `tmx` → variable `tasmax`.


#### Management parameters

Management parameters are specified in `[input.netcdf.<param_name>]` sections with a `default_value` key. Optionally, you can also provide a spatial NC file.

1. **Uniform default only**: set `default_value` — used for all pixels when no NC file is provided, or the file/year is not found.
2. **Spatial NC file**: add `file` and `variable` keys alongside `default_value`. The file can be 2D (lon, lat) for static parameters, or 3D (lon, lat, time) with a year dimension for time-varying parameters.

Priority: NC file > `default_value` > built-in default.

For full details, see [README-CONFIG.md](README-CONFIG.md).
