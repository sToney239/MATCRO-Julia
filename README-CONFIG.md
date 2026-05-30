# Configuration Reference

This document describes all parameters in the MATCRO-Julia configuration file (`config.toml`).

For a quick start guide, see [README.md](README.md).

---

## 1. `[general]`

Core simulation settings that apply to both point and spatial modes.

| Parameter | Type | Default | Unit | Description |
|-----------|------|---------|------|-------------|
| `crop_name` | String | — | — | Crop type: `"Maize"`, `"Rice"`, `"Wheat"`, or `"Soybeans"` |
| `crop_param` | String | — | — | Path to crop parameter TOML file (relative to config file) |
| `start_year` | Int | — | year | Simulation start year |
| `end_year` | Int | — | year | Simulation end year |
| `start_doy` | Int | 1 | day of year | Simulation start day of year. Must be ≥ 1 and < `end_doy`. |
| `end_doy` | Int | 365 | day of year | Simulation end day of year. Must be ≥ 1 and > `start_doy`. |
| `time_step` | Int | 3600 | seconds | Time step per iteration, often for 1 hour which should be 3600s. |
| `co2_ppm_default` | Float64 | 400.0 | ppm | Default CO₂ concentration. |
| `co2_yearly_file` | String | "" | — | Path to CO₂ CSV file (relative to config file, or absolute). Takes priority over `co2_ppm_default` |

**`start_doy` / `end_doy` behavior**:
- **Point mode**: defines the DOY range for the first and last simulation year. Middle years run DOY 1–365.
- **Spatial mode**: the DOY range is used to select a subset of days from the NetCDF forcing data. If the NetCDF file contains fewer days than `end_doy`, the actual range is clamped to the available data. By default (`start_doy=1`, `end_doy=365`), all days in the NetCDF file are used.

**CO₂ file format**: CSV with header row and two columns:

```
year,co2_ppm
2000,369.4
2001,371.1
...
```

If the file does not exist or the requested year is not found, `co2_ppm_default` is used.

---

## 2. Simulation: choose `[point_simulation]` **or** `[spatial_simulation]`

You must specify exactly one simulation mode. If both are present, `[point_simulation]` takes priority.

---

## 2.1 `[point_simulation]`

Point simulation mode with CSV weather data.

### Global settings

| Parameter | Type | Default | Unit | Description |
|-----------|------|---------|------|-------------|
| `latitude` | Float64 | — | degrees | Site latitude (positive = North) |
| `output_directory` | String | "output/" | — | Output directory path (relative to config file, or absolute) |

### `[point_simulation.weather]`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `csv_path` | String | — | Path to CSV weather data file (relative to config file) |

**CSV file format**: Must include a header row with the following columns:

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| `year` | Int | — | Year |
| `doy` | Int | — | Day of year (1–366) |
| `tmax` | Float64 | K | Daily maximum air temperature |
| `tmin` | Float64 | K | Daily minimum air temperature |
| `radiation` | Float64 | W/m² | Downward shortwave radiation |
| `precip` | Float64 | kg/m²/s | Precipitation rate |
| `humidity` | Float64 | kg/kg | Specific humidity |
| `wind` | Float64 | m/s | Wind speed |
| `pressure` | Float64 | Pa | Surface air pressure |
| `ozone` | Float64 | — | Ozone concentration (currently not in the model) |

Example header and first row:

```
year,doy,tmax,tmin,radiation,precip,humidity,wind,pressure,ozone
2021,1,290.15,276.58,104.51,0.0,0.00582,2.57,98652.0,26.0
```

### `[point_simulation.management]`

Management parameters for point simulation.

| Parameter | Type | Default | Unit | Description |
|-----------|------|---------|------|-------------|
| `planting_doy` | Int | 120 | day of year | Day of year for planting |
| `is_irrigated` | Int | 0 | — | 0 = rainfed, 1 = irrigated |
| `soil_type` | Int | 9 | — | Soil texture index (1–13). 9 = loam |
| `n_fertilizer` | Float64 | 100.0 | kg N/ha | Nitrogen fertilizer application rate |
| `thermal_time_requirement` | Float64 | 1500.0 | °C·day | Growing degree days required for maturity |

Output: CSV files (`yield_summary.csv`, `daily_output.csv`) in `output_directory`.

---

## 2.2 `[spatial_simulation]`

Spatial simulation mode with NetCDF/TIFF data.

### Global settings

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `output_directory` | String | "output/" | Output directory path (relative to config file, or absolute) |

Foe this spatial simulation section, the model outputs **one GeoTIFF file per year**. Four files are generated for each simulation year:

- `yield_YYYY.tif`: Crop yield (kg/ha), Float64
- `harvest_doy_YYYY.tif`: Harvest day of year, Int32
- `LAI_max_YYYY.tif`: Seasonal maximum leaf area index (m²/m²), Float64
- `biomass_aboveground_YYYY.tif`: Aboveground biomass at harvest (kg/ha), Float64

The TIF files use WGS84 (EPSG:4326) coordinate reference system, with geotransform metadata for proper geospatial alignment.

### 2.3.1 Weather: `[spatial_simulation.weather]`

Dimension name settings for weather NetCDF files:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `lon_dim` | String | `"lon"` | Longitude dimension name in NetCDF files |
| `lat_dim` | String | `"lat"` | Latitude dimension name in NetCDF files |
| `time_dim` | String | `"time"` | Time dimension name in NetCDF files |

Dimension names are configurable via `lon_dim`, `lat_dim`, `time_dim` in `[spatial_simulation.weather]` (defaults: `lon`, `lat`, `time`). Common alternatives like `latitude`, `longitude` are also auto-detected as fallbacks.

The time variable in your NetCDF file must follow one of these formats:

1. **Stored as dates**: time values are directly stored as calendar dates (e.g., `2021-01-01`, `2021-01-02`, ...)
2. **Stored as numbers with a `units` attribute**: numeric values (e.g., 0, 1, 2, ...) with a `units` attribute following the format `"days since YYYY-MM-DD"` (e.g., `"days since 2021-01-01"`). The program converts each number to a date using the reference date in `units`, then slices by year.

> If the `units` attribute is missing or not in the `"days since YYYY-MM-DD"` format, the program will report an error.

#### Per-variable specification: `[spatial_simulation.weather.<var>]`

Each meteorological variable is specified in its own section under `[spatial_simulation.weather]`. Use the following user-friendly variable names:

| Config key | Internal name | Description |
|------------|---------------|-------------|
| `temperature_max` | tmx | Daily maximum temperature |
| `temperature_min` | tmn | Daily minimum temperature |
| `precipitation` | prc | Precipitation rate |
| `radiation` | rsd | Downward shortwave radiation |
| `humidity` | shm | Specific humidity |
| `wind_speed` | wnd | Wind speed |
| `pressure` | prs | Surface air pressure |

**Required parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `file` | String | — | Path to NetCDF file (relative to config file, or absolute) |
| `variable` | String | — | Variable name inside the NetCDF file (may differ from the section key, e.g., key `temperature_max` → variable `tasmax`) |

**Optional parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `scale_factor` | Float64 | 1.0 | Scale factor applied as `value * scale_factor + add_offset` |
| `add_offset` | Float64 | 0.0 | Offset applied as `value * scale_factor + add_offset` |
| `height` | Float64 | 10.0 | Wind measurement height in meters (only for `wind_speed`) |

Example:

```toml
[spatial_simulation.weather.temperature_max]
file = "data/netcdf/tmx.nc4"
variable = "tmx"

[spatial_simulation.weather.wind_speed]
file = "data/netcdf/wnd.nc4"
variable = "wnd"
height = 10.0
```

### 2.3.2 Management: `[spatial_simulation.management]`

Management parameters for spatial simulation.

#### Per-parameter specification: `[spatial_simulation.management.<param_name>]`

Supported `<param_name>`: `planting_doy`, `is_irrigated`, `soil_type`, `n_fertilizer`, `thermal_time_requirement`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `default_value` | Int/Float64 | see below | Default value used when no file is provided, or the file/year is not found |
| `file` | String | — | (Optional) Path to spatial file (NetCDF `.nc` or GeoTIFF `.tif`/`.tiff`) |
| `variable` | String | — | (Required if `file` is set and is a NetCDF file) Variable name inside the NetCDF file |

**Built-in defaults** (used when `[spatial_simulation.management.<param_name>]` section is not present at all):

| Parameter | Default |
|-----------|---------|
| `planting_doy` | 120 |
| `is_irrigated` | 0 |
| `soil_type` | 9 |
| `n_fertilizer` | 100.0 |
| `thermal_time_requirement` | 1500.0 |

The file can be:
- **NetCDF 2D** (lon, lat): static parameter, same value for all years
- **NetCDF 3D** (lon, lat, time or year): time-varying parameter, sliced by year
- **GeoTIFF 2D** (single band): static parameter, same value for all years

When using a GeoTIFF file:
- Nearest-neighbor resampling is used when TIF and simulation grid resolutions differ
- A bbox check is performed to ensure the TIF extent covers the simulation grid — if not, a warning is issued and pixels outside the TIF coverage use `default_value`

The time dimension in management parameter NetCDF files can be named either `time` or `year`. The program will try `time` first, then fall back to `year`.

**Year matching logic** (for 3D management parameter files):

1. **Exact match**: if the requested simulation year exists in the file, use that year's data
2. **Backward fill**: if the exact year is not found but earlier years exist, use the nearest earlier year (e.g., file has 2020 and 2022, simulation requests 2021 → use 2020's data)
3. **Single value fallback**: if the file has only one time value, use it regardless of the simulation year (with a console warning)
4. **default_value**: if none of the above works, fall back to `default_value`

Priority: file (NC or TIF, with matching logic above) > `default_value` > built-in default.

Example — uniform value only (no spatial file):

```toml
[spatial_simulation.management.planting_doy]
default_value = 120
```

Example — spatial NC file with fallback:

```toml
[spatial_simulation.management.soil_type]
default_value = 9
file = "data/netcdf/soil_type.nc"
variable = "soil_type"
```

Example — spatial TIF file with fallback:

```toml
[spatial_simulation.management.soil_type]
default_value = 9
file = "data/tif/soil_type.tif"
```

> For TIF files, the `variable` key is not needed — the program reads the first band directly.

### 2.3.3 Boundary: `[spatial_simulation.boundary]` (optional)

Boundary file for spatial filtering — only pixels within (or contacting) the boundary are simulated.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `file` | String | — | Path to boundary file (GeoJSON `.geojson`) |
| `buffer_deg` | Float64 | 0.0 | Buffer distance (in degrees) for contact detection. Pixels within this distance of the boundary are also included. Default 0 means strict inside-polygon check. |

**Example**:

```toml
[spatial_simulation.boundary]
file = "data/boundaries/cornbelt_states.geojson"
buffer_deg = 0.01   # ~1km buffer at mid-latitudes
```

The boundary CRS is automatically converted to WGS84 (EPSG:4326) if it differs from the input data.

Output: GeoTIFF files per year in `output_directory`:
- `yield_YYYY.tif`: Crop yield (kg/ha), Float64
- `harvest_doy_YYYY.tif`: Harvest day of year, Int32
- `LAI_max_YYYY.tif`: Seasonal maximum leaf area index (m²/m²), Float64
- `biomass_aboveground_YYYY.tif`: Aboveground biomass at harvest (kg/ha), Float64

All files use WGS84 (EPSG:4326) CRS with proper geotransform metadata.

