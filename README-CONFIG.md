# Configuration Reference

This document describes all parameters in the MATCRO-Julia configuration file (`config.toml`).

For a quick start guide, see [README.md](README.md).

---

## `[time]`

| Parameter | Type | Default | Unit | Description |
|-----------|------|---------|------|-------------|
| `start_year` | Int | — | year | Simulation start year |
| `end_year` | Int | — | year | Simulation end year |
| `start_doy` | Int | 1 | day of year | Simulation start day of year |
| `end_doy` | Int | 365 | day of year | Simulation end day of year |
| `time_step` | Int | 3600 | seconds | Time step per iteration |

---

## `[location]`

| Parameter | Type | Default | Unit | Description |
|-----------|------|---------|------|-------------|
| `latitude` | Float64 | — | degrees | Site latitude (positive = North) |

> Longitude is not needed here. In NetCDF mode, it is read from the input data. In CSV mode, longitude is not required.

---

## `[crop]`

| Parameter | Type | Default | Unit | Description |
|-----------|------|---------|------|-------------|
| `crop_name` | String | — | — | Crop type: `"Maize"`, `"Rice"`, `"Wheat"`, or `"Soybeans"` |
| `param_file` | String | — | — | Path to crop parameter TOML file (relative to config file) |
| `planting_doy` | Int | — | day of year | Day of year for planting |
| `is_irrigated` | Int | 0 | — | 0 = rainfed, 1 = irrigated |

---

## `[soil]`

| Parameter | Type | Default | Unit | Description |
|-----------|------|---------|------|-------------|
| `soil_type` | Int | 9 | — | Soil texture index (1–13). 9 = loam |
| `n_fertilizer` | Float64 | 100.0 | kg N/ha | Nitrogen fertilizer application rate |
| `thermal_time_requirement` | Float64 | 1500.0 | °C·day | Growing degree days required for maturity |

---

## `[co2]`

| Parameter | Type | Default | Unit | Description |
|-----------|------|---------|------|-------------|
| `fixed_ppm` | Float64 | 400.0 | ppm | Fallback CO₂ concentration |
| `file` | String | "" | — | Path to CO₂ CSV file (relative to config file). Takes priority over `fixed_ppm` |

**CO₂ file format**: CSV with header row and two columns:

```
year,co2_ppm
2000,369.4
2001,371.1
...
```

If the file does not exist or the requested year is not found, `fixed_ppm` is used.

---

## `[input.csv]`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `path` | String | — | Path to CSV forcing data file (relative to config file) |

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
| `ozone` | Float64 | — | Ozone concentration (can be 0) |

Example header and first row:

```
year,doy,tmax,tmin,radiation,precip,humidity,wind,pressure,ozone
2021,1,290.15,276.58,104.51,0.0,0.00582,2.57,98652.0,26.0
```

---

## `[input.netcdf]`

Choose either `[input.csv]` **or** `[input.netcdf]` — if you have chosen one, you can delete the lines about the other or just comment them out.

### Global settings

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `nthreads` | Int | 1 | Number of parallel threads for spatial simulation |
| `lon_dim` | String | `"lon"` | Longitude dimension name in NetCDF files |
| `lat_dim` | String | `"lat"` | Latitude dimension name in NetCDF files |
| `time_dim` | String | `"time"` | Time dimension name in NetCDF files |

### Per-variable specification: `[input.netcdf.<var>]`

Each meteorological variable is specified in its own section. The key name (e.g., `tmx`) is the internal variable name used by MATCRO.

Supported `<var>` variables: `tmx`, `tmn`, `prc`, `rsd`, `shm`, `wnd`, `prs`

**Required parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `file` | String | — | Path to NetCDF file (relative to config file, or absolute) |
| `variable` | String | — | Variable name inside the NetCDF file (may differ from the section key, e.g., key `tmx` → variable `tasmax`) |

**Optional parameters**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `scale_factor` | Float64 | 1.0 | Scale factor applied as `value * scale_factor + add_offset` |
| `add_offset` | Float64 | 0.0 | Offset applied as `value * scale_factor + add_offset` |
| `height` | Float64 | 10.0 | Wind measurement height in meters (only for `wnd`) |

### Management parameters: `[input.netcdf.<param_name>]`

Management parameters are also specified in their own sections under `[input.netcdf]`.

Supported `<param_name>`: `planting_doy`, `is_irrigated`, `soil_type`, `n_fertilizer`, `thermal_time_requirement`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `default_value` | Int/Float64 | see below | Default value used when no file is provided, or the file/year is not found |
| `file` | String | — | (Optional) Path to spatial file (NetCDF `.nc` or GeoTIFF `.tif`/`.tiff`) |
| `variable` | String | — | (Required if `file` is set and is a NetCDF file) Variable name inside the NetCDF file |

**Built-in defaults** (used when `[input.netcdf.<param_name>]` section is not present at all):

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

Example — uniform value only (no spatial NC file):

```toml
[input.netcdf.planting_doy]
default_value = 120
```

Example — spatial NC file with fallback:

```toml
[input.netcdf.soil_type]
default_value = 9
file = "data/netcdf/soil_type.nc"
variable = "soil_type"
```

Example — spatial TIF file with fallback:

```toml
[input.netcdf.soil_type]
default_value = 9
file = "data/tif/soil_type.tif"
```

> For TIF files, the `variable` key is not needed — the program reads the first band directly.

---

## `[output]`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `directory` | String | — | Output directory path (relative to config file, or absolute) |
| `format` | String | — | Output format: `"csv"` or `"raster"` (or `"geotiff"`) |

When `format = "raster"` (or `"geotiff"`), spatial simulation outputs GeoTIFF files per year:
- `yield_YYYY.tif` — Crop yield (Float64, kg/ha)
- `harvest_doy_YYYY.tif` — Harvest day of year (Int32)

Both files use WGS84 (EPSG:4326) CRS with proper geotransform metadata.
