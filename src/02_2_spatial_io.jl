# IO_SPATIAL — Spatial I/O for MATCRO simulation (NetCDF + GeoTIFF)
# Used when config.input_format = "raster"
# Supports: single-point, spatial parallel, multi-year single file
# Companion: io.jl (core structs), io_csv.jl (CSV I/O)
#
# Each variable's file path is specified individually in config.toml [input.raster.weather] or [input.raster.management]
# No base_dir — paths resolved relative to config file directory by read_config

using NCDatasets
using Dates
using ArchGDAL
using JSON
import ArchGDAL as AG

# ============================================================
# read_forcing_netcdf — load one year of daily forcing for a
# single lat/lon point from NetCDF files, returns Dict{doy => DailyForcing}
# ============================================================
function read_forcing_netcdf(config::Config, year::Int)::Dict{Int,DailyForcing}
    nc = config.raster_vars
    lat = config.latitude
    lon = Float64(get(nc, "longitude", -100.0))

    result = Dict{Int,DailyForcing}()

    tmax_yr  = _read_nc_point(nc["tmx"], year, lat, lon; config=config)
    tmin_yr  = _read_nc_point(nc["tmn"], year, lat, lon; config=config)
    rsd_yr   = _read_nc_point(nc["rsd"], year, lat, lon; config=config)
    prc_yr   = _read_nc_point(nc["prc"], year, lat, lon; config=config)
    shm_yr   = _read_nc_point(nc["shm"], year, lat, lon; config=config)
    wnd_yr   = _read_nc_point(nc["wnd"], year, lat, lon; config=config)
    prs_yr   = _read_nc_point(nc["prs"], year, lat, lon; config=config)

    n_days = minimum([length(tmax_yr), length(tmin_yr), length(rsd_yr),
                      length(prc_yr), length(shm_yr), length(wnd_yr), length(prs_yr)])

    for doy in 1:n_days
        result[doy] = DailyForcing(;
            doy       = doy,
            tmax      = tmax_yr[doy],
            tmin      = tmin_yr[doy],
            radiation = rsd_yr[doy],
            precip    = prc_yr[doy],
            humidity  = shm_yr[doy],
            wind      = wnd_yr[doy],
            pressure  = prs_yr[doy],
        )
    end

    return result
end

# ============================================================
# _read_nc_point — read a 1-year time series at a single lat/lon
# Slices by time dimension using _find_year_indices
# ============================================================
function _read_nc_point(var_meta::Dict{String,Any}, year::Int,
                        target_lat::Float64, target_lon::Float64;
                        config::Union{Config,Nothing}=nothing)::Vector{Float64}
    filepath = var_meta["file"]

    ds = Dataset(filepath)

    lon_dim = config !== nothing ? config.lon_dim : "lon"
    lat_dim = config !== nothing ? config.lat_dim : "lat"
    time_dim = config !== nothing ? config.time_dim : "time"

    lat_var = _find_var(ds, lat_dim, ["lat", "latitude", "latitute"])
    lon_var = _find_var(ds, lon_dim, ["lon", "longitude"])

    lat_vals = Float64.(lat_var[:])
    lon_vals = Float64.(lon_var[:])

    ilat = argmin(abs.(lat_vals .- target_lat))
    ilon = argmin(abs.(lon_vals .- target_lon))

    vname = var_meta["variable"]
    var_data = ds[vname]

    scale_factor = Float64(get(var_meta, "scale_factor", 1.0))
    add_offset = Float64(get(var_meta, "add_offset", 0.0))

    result = Float64[]

    # Slice by year using time dimension
    year_indices = _find_year_indices(ds, time_dim, year)
    for t in year_indices
        if ndims(var_data) == 3
            val = Float64(var_data[ilon, ilat, t])
        else
            val = Float64(var_data[t])
        end
        push!(result, val * scale_factor + add_offset)
    end

    close(ds)
    return result
end

# ============================================================
# _find_year_indices — find time indices for a given year in a multi-year NC file
# ============================================================
function _find_year_indices(ds::NCDataset, time_dim::String, target_year::Int)::Vector{Int}
    time_var = ds[time_dim]
    time_vals = time_var[:]

    if eltype(time_vals) <: DateTime
        return findall(t -> year(t) == target_year, time_vals)
    else
        # Numeric time with units like "days since 2021-01-01"
        time_units = get(time_var.attrib, "units", "")
        if occursin("since", time_units)
            ref_str = strip(split(time_units, "since")[2])
            ref_date = DateTime(ref_str, "yyyy-mm-dd")
            indices = Int[]
            for (i, t) in enumerate(time_vals)
                dt = ref_date + Day(round(Int, Float64(t)))
                if year(dt) == target_year
                    push!(indices, i)
                end
            end
            return indices
        else
            error("Cannot parse time units: $time_units")
        end
    end
end

# ============================================================
# _find_var — find a variable in a Dataset by name, with fallbacks
# ============================================================
function _find_var(ds::NCDataset, preferred::String, fallbacks::Vector{String})
    if haskey(ds, preferred)
        return ds[preferred]
    end
    for name in fallbacks
        if haskey(ds, name)
            return ds[name]
        end
    end
    error("Cannot find variable '$preferred' or any of $fallbacks in dataset")
end

# ============================================================
# get_grid_info — read lon/lat arrays from a NetCDF file
# Uses the first variable's file to get grid info
# ============================================================
function get_grid_info(config::Config)
    nc = config.raster_vars
    # Always read grid from tmx forcing file (all forcing files share the same grid)
    # This ensures consistent lat/lon ordering regardless of management param format
    tmx_meta = get(nc, "tmx", nothing)
    if tmx_meta === nothing || !haskey(tmx_meta, "file")
        error("Cannot find tmx forcing file for grid info")
    end
    filepath = tmx_meta["file"]

    ds = Dataset(filepath)
    lat_var = _find_var(ds, config.lat_dim, ["lat", "latitude", "latitute"])
    lon_var = _find_var(ds, config.lon_dim, ["lon", "longitude"])
    lats = Float64.(lat_var[:])
    lons = Float64.(lon_var[:])
    close(ds)
    return lats, lons
end

# ============================================================
# read_forcing_netcdf_spatial — load one year of spatial forcing
# Returns Dict{doy => (tmax, tmin, rsd, prc, shm, wnd, prs)} where
# each field is Matrix{Float64} of size (nlon, nlat)
# ============================================================
function read_forcing_netcdf_spatial(config::Config, year::Int)
    nc = config.raster_vars

    tmax_3d = _read_nc_spatial(nc["tmx"], year; config=config)
    tmin_3d = _read_nc_spatial(nc["tmn"], year; config=config)
    rsd_3d  = _read_nc_spatial(nc["rsd"], year; config=config)
    prc_3d  = _read_nc_spatial(nc["prc"], year; config=config)
    shm_3d  = _read_nc_spatial(nc["shm"], year; config=config)
    wnd_3d  = _read_nc_spatial(nc["wnd"], year; config=config)
    prs_3d  = _read_nc_spatial(nc["prs"], year; config=config)

    n_lon, n_lat, n_days = size(tmax_3d)

    result = Dict{Int, @NamedTuple{tmax::Matrix{Float64}, tmin::Matrix{Float64},
                                    radiation::Matrix{Float64}, precip::Matrix{Float64},
                                    humidity::Matrix{Float64}, wind::Matrix{Float64},
                                    pressure::Matrix{Float64}}}()

    for doy in 1:n_days
        result[doy] = (
            tmax      = tmax_3d[:, :, doy],
            tmin      = tmin_3d[:, :, doy],
            radiation = rsd_3d[:, :, doy],
            precip    = prc_3d[:, :, doy],
            humidity  = shm_3d[:, :, doy],
            wind      = wnd_3d[:, :, doy],
            pressure  = prs_3d[:, :, doy],
        )
    end

    return result
end

# ============================================================
# _read_nc_spatial — read one year of spatial data from a NetCDF variable
# Returns Array{Float64, 3} of size (nlon, nlat, ndays)
# ============================================================
function _read_nc_spatial(var_meta::Dict{String,Any}, year::Int;
                          config::Union{Config,Nothing}=nothing)
    filepath = var_meta["file"]

    ds = Dataset(filepath)

    vname = var_meta["variable"]
    var_data = ds[vname]

    scale_factor = Float64(get(var_meta, "scale_factor", 1.0))
    add_offset = Float64(get(var_meta, "add_offset", 0.0))

    # Slice by year using time dimension
    time_dim = config !== nothing ? config.time_dim : "time"
    year_indices = _find_year_indices(ds, time_dim, year)
    result = Float64.(var_data[:, :, year_indices])

    close(ds)

    result .= result .* scale_factor .+ add_offset
    return result
end

# ============================================================
# load_management_param — load a management parameter for spatial mode
# Reads from NC file if file/variable keys exist in nc_vars.
# Otherwise returns a uniform matrix filled with default_value.
#
# Time matching logic:
# 1. Try exact year match
# 2. If no match but earlier years exist, use nearest earlier year (backward fill)
# 3. If only one time value exists, use it with a warning
# 4. If neither works, fall back to default_value
# ============================================================
function load_management_param(config::Config, param_name::String, year::Int,
                               n_lon::Int, n_lat::Int; lats=nothing, lons=nothing)
    var_meta = get(config.raster_vars, param_name, nothing)
    default_val = var_meta !== nothing ? Float64(get(var_meta, "default_value", 0)) : 0.0

    if var_meta !== nothing && haskey(var_meta, "file")
        filepath = var_meta["file"]

        if !isfile(filepath)
            @warn "Management param file not found: $filepath, using default=$default_val"
            return fill(default_val, n_lon, n_lat)
        end

        # Check if file is TIF based on extension
        if endswith(lowercase(filepath), ".tif") || endswith(lowercase(filepath), ".tiff")
            if lats === nothing || lons === nothing
                @warn "TIF management param requires lat/lon info, using default=$default_val"
                return fill(default_val, n_lon, n_lat)
            end
            return _load_management_tif(filepath, lats, lons, default_val)
        end

        ds = Dataset(filepath)
        vname = haskey(var_meta, "variable") ? var_meta["variable"] : param_name

        # Check lat/lon ordering in this NC file (before close)
        nc_lat_descending = false
        nc_lat_vals = Float64[]
        nc_lon_vals = Float64[]
        if haskey(ds, config.lat_dim) || haskey(ds, "lat") || haskey(ds, "latitude")
            nc_lat_var = _find_var(ds, config.lat_dim, ["lat", "latitude", "latitute"])
            nc_lat_vals = Float64.(nc_lat_var[:])
            if length(nc_lat_vals) >= 2
                nc_lat_descending = nc_lat_vals[1] > nc_lat_vals[2]
            end
        end
        if haskey(ds, config.lon_dim) || haskey(ds, "lon") || haskey(ds, "longitude")
            nc_lon_var = _find_var(ds, config.lon_dim, ["lon", "longitude"])
            nc_lon_vals = Float64.(nc_lon_var[:])
        end

        # Try "time" first, then fallback to "year" (management params often use "year")
        time_dim = haskey(ds, config.time_dim) ? config.time_dim : (haskey(ds, "year") ? "year" : nothing)

        if time_dim !== nothing
            time_var = ds[time_dim]
            time_vals = time_var[:]

            if eltype(time_vals) <: DateTime
                available_years = year.(time_vals)
            else
                # Numeric time - try to extract years
                time_units = get(time_var.attrib, "units", "")
                if occursin("since", time_units)
                    ref_str = strip(split(time_units, "since")[2])
                    ref_date = DateTime(ref_str, "yyyy-mm-dd")
                    available_years = [year(ref_date + Day(round(Int, t))) for t in time_vals]
                else
                    # Assume it's already years if not "days since..."
                    available_years = Int.(time_vals)
                end
            end

            # Case 1: Exact year match
            exact_matches = findall(y -> y == year, available_years)
            if !isempty(exact_matches)
                data = Float64.(ds[vname][:, :, exact_matches[1]])
            elseif length(available_years) == 1
                # Case 3: Only one time value - use it with warning
                @warn "Management param '$param_name' in $filepath has only one time value (year $(available_years[1])), using it regardless of simulation year $year"
                data = Float64.(ds[vname][:, :])
            else
                # Case 2: No exact match - try backward fill (nearest earlier year)
                earlier_years = available_years[available_years .< year]
                if !isempty(earlier_years)
                    nearest_idx = argmax(earlier_years)  # Get the closest (largest) earlier year
                    nearest_year = earlier_years[nearest_idx]
                    actual_idx = findall(y -> y == nearest_year, available_years)[1]
                    @warn "Year $year not found in management param '$param_name' ($filepath), using nearest earlier year $nearest_year (backward fill)"
                    data = Float64.(ds[vname][:, :, actual_idx])
                else
                    @warn "No suitable year found in management param '$param_name' ($filepath) for simulation year $year, using default=$default_val"
                    close(ds)
                    return fill(default_val, n_lon, n_lat)
                end
            end
        else
            # No time/year dimension - treat as 2D static file
            data = Float64.(ds[vname][:, :])
        end
        close(ds)

        scale_factor = Float64(get(var_meta, "scale_factor", 1.0))
        add_offset = Float64(get(var_meta, "add_offset", 0.0))
        data .= data .* scale_factor .+ add_offset

        # Reindex lat dimension to match grid (forcing) lat ordering
        # Management NC files may have lat descending while forcing/grid has ascending
        # Flip when: NC is descending AND grid is ascending (or vice versa)
        if lats !== nothing && size(data, 2) == length(lats)
            grid_lat_ascending = lats[1] < lats[end]
            # Flip if: (NC descending AND grid ascending) OR (NC ascending AND grid descending)
            should_flip = (nc_lat_descending && grid_lat_ascending) || (!nc_lat_descending && !grid_lat_ascending)
            if should_flip
                data = data[:, end:-1:1]
            end
        end

        # Bbox coverage check: warn if management NC does not fully cover weather grid
        if lats !== nothing && lons !== nothing && length(nc_lat_vals) >= 2
            nc_lon_min, nc_lon_max = extrema(nc_lon_vals)
            nc_lat_min, nc_lat_max = nc_lat_descending ? (nc_lat_vals[end], nc_lat_vals[1]) : (nc_lat_vals[1], nc_lat_vals[end])
            grid_lon_min, grid_lon_max = extrema(lons)
            grid_lat_min, grid_lat_max = extrema(lats)
            if nc_lon_min > grid_lon_min || nc_lon_max < grid_lon_max ||
               nc_lat_min > grid_lat_min || nc_lat_max < grid_lat_max
                @warn "Management param '$param_name' ($filepath) extent does not fully cover weather grid — uncovered pixels using default values"
            end
        end

        return data
    else
        return fill(default_val, n_lon, n_lat)
    end
end

# ============================================================
# _load_management_tif — load management parameter from GeoTIFF
# Only single-band (2D) TIF is supported; the value is used for all years
# ============================================================
function _load_management_tif(filepath::String,
                              target_lats::Vector{Float64}, target_lons::Vector{Float64},
                              default_val::Float64)::Matrix{Float64}
    AG.read(filepath) do dataset
        _extract_tif_band(dataset, 1, target_lons, target_lats, default_val)
    end
end

# ============================================================
# _extract_tif_band — extract one band from a TIF, reproject to NC grid
# ============================================================
function _extract_tif_band(dataset, band_idx::Int,
                            target_lons::Vector{Float64},
                            target_lats::Vector{Float64},
                            default_val::Float64)::Matrix{Float64}
    n_lon = length(target_lons)
    n_lat = length(target_lats)
    band = AG.getband(dataset, band_idx)
    tif_width = AG.width(dataset)
    tif_height = AG.height(dataset)
    gt = AG.getgeotransform(dataset)
    x_origin = gt[1]; pixel_w = gt[2]
    y_origin = gt[4]; pixel_h = gt[6]

    # Bbox check
    tif_xmin = x_origin
    tif_xmax = x_origin + tif_width * pixel_w
    tif_ymin = y_origin + tif_height * pixel_h
    tif_ymax = y_origin

    nc_xmin = minimum(target_lons) - (target_lons[2] - target_lons[1]) / 2
    nc_xmax = maximum(target_lons) + (target_lons[2] - target_lons[1]) / 2
    nc_ymin = minimum(target_lats) - (target_lats[1] - target_lats[2]) / 2
    nc_ymax = maximum(target_lats) + (target_lats[1] - target_lats[2]) / 2

    if tif_xmin > nc_xmin || tif_xmax < nc_xmax || tif_ymin > nc_ymin || tif_ymax < nc_ymax
        @warn "Management TIF bbox does not fully cover NC grid — some pixels may use defaults"
    end

    result = fill(default_val, n_lon, n_lat)
    # Read entire band first, then index
    tif_data = Float64.(AG.read(band))
    for i_lon in 1:n_lon, i_lat in 1:n_lat
        # GeoTIFF: pixel center at (x_origin + pixel_w*(col-0.5), y_origin + pixel_h*(row-0.5))
        col = round(Int, (target_lons[i_lon] - x_origin) / pixel_w + 0.5)
        row = round(Int, (target_lats[i_lat] - y_origin) / pixel_h + 0.5)
        if 1 <= col <= tif_width && 1 <= row <= tif_height
            # tif_data is [col, row] = [lon, lat]
            result[i_lon, i_lat] = tif_data[col, row]
        end
    end
    return result
end

# ============================================================
# write_output_netcdf — write single-point results to NetCDF
# ============================================================
function write_output_netcdf(results::Vector{NamedTuple}, output_path::String)
    if isempty(results)
        @warn "No results to write"
        return
    end

    n = length(results)
    years = [r.year for r in results]
    doys  = [r.doy for r in results]

    ds = Dataset(output_path, "c")

    defDim(ds, "time", n)

    defVar(ds, "year", years, ("time",))
    defVar(ds, "doy", doys, ("time",))

    if haskey(results[1], :yield)
        defVar(ds, "yield", Float64, ("time",), attrib = Dict("units" => "kg/ha"))
        ds["yield"][:] = [Float64(r.yield) for r in results]
    end
    if haskey(results[1], :LAI)
        defVar(ds, "LAI", Float64, ("time",), attrib = Dict("units" => "m2/m2"))
        ds["LAI"][:] = [Float64(r.LAI) for r in results]
    end
    if haskey(results[1], :development_stage)
        defVar(ds, "DVS", Float64, ("time",), attrib = Dict("units" => "-"))
        ds["DVS"][:] = [Float64(r.development_stage) for r in results]
    end
    if haskey(results[1], :biomass_aboveground)
        defVar(ds, "biomass_aboveground", Float64, ("time",), attrib = Dict("units" => "kg/ha"))
        ds["biomass_aboveground"][:] = [Float64(r.biomass_aboveground) for r in results]
    end

    close(ds)
    println("Output written to ", output_path)
end

# ============================================================
# write_output_netcdf_spatial — write spatial yield to 3D NetCDF
# (lon, lat, year) — appends year dimension incrementally
# ============================================================
function write_output_netcdf_spatial(yield_3d::Array{Float64,3},
                                     lats::Vector{Float64},
                                     lons::Vector{Float64},
                                     years::Vector{Int},
                                     output_path::String)
    n_lon = length(lons)
    n_lat = length(lats)
    n_year = length(years)

    ds = Dataset(output_path, "c")
    defDim(ds, "lon", n_lon)
    defDim(ds, "lat", n_lat)
    defDim(ds, "year", n_year)

    defVar(ds, "lon", lons, ("lon",), attrib = Dict("units" => "degrees_east"))
    defVar(ds, "lat", lats, ("lat",), attrib = Dict("units" => "degrees_north"))
    defVar(ds, "year", years, ("year",), attrib = Dict("units" => "year"))
    defVar(ds, "yield", Float64, ("lon", "lat", "year"),
           attrib = Dict("long_name" => "crop yield", "units" => "kg/ha"))

    # Write the actual 3D data
    ds["yield"][:, :, :] = yield_3d

    close(ds)
    println("Output written to ", output_path)
end

# ============================================================
# write_yield_slice — write one year's yield slice to an existing NC
# ============================================================
function write_yield_slice(output_path::String, year::Int, years::Vector{Int},
                           yield_map::Matrix{Float64})
    iy = findfirst(==(year), years)
    if iy === nothing
        error("Year $year not found in years vector")
    end

    ds = Dataset(output_path, "a")
    ds["yield"][:, :, iy] = yield_map
    close(ds)
end

# ============================================================
# write_float64_tif — write a 2D Float64 GeoTIFF for one year
# ============================================================
function write_float64_tif(data_map::Matrix{Float64}, lats::Vector{Float64},
                           lons::Vector{Float64}, year::Int, output_path::String;
                           unit::String="kg/ha")
    n_lon = length(lons)
    n_lat = length(lats)

    pixel_width = lons[2] - lons[1]    # positive
    pixel_height = lats[2] - lats[1]   # negative for north-up

    AG.create(output_path; driver=AG.getdriver("GTiff"), width=n_lon, height=n_lat, nbands=1, dtype=Float64) do dataset
        band = AG.getband(dataset, 1)
        AG.write!(band, data_map)

        y_origin = lats[1] - pixel_height / 2
        geo_transform = [lons[1] - pixel_width / 2, pixel_width, 0.0,
                         y_origin, 0.0, pixel_height]
        AG.setgeotransform!(dataset, geo_transform)
        AG.setproj!(dataset, "EPSG:4326")

        AG.setunittype!(band, unit)
        AG.setnodatavalue!(band, NaN)
    end
end

# ============================================================
# write_harvest_doy_tif — write a 2D harvest_doy GeoTIFF for one year
# ============================================================
function write_harvest_doy_tif(harvest_map::Matrix{Float64}, lats::Vector{Float64},
                               lons::Vector{Float64}, year::Int, output_path::String)
    n_lon = length(lons)
    n_lat = length(lats)

    pixel_width = lons[2] - lons[1]
    pixel_height = lats[2] - lats[1]

    AG.create(output_path; driver=AG.getdriver("GTiff"), width=n_lon, height=n_lat, nbands=1, dtype=Int32) do dataset
        band = AG.getband(dataset, 1)
        AG.write!(band, Int32.(round.(replace(harvest_map, NaN => -999))))

        y_origin = lats[1] - pixel_height / 2
        geo_transform = [lons[1] - pixel_width / 2, pixel_width, 0.0,
                         y_origin, 0.0, pixel_height]
        AG.setgeotransform!(dataset, geo_transform)
        AG.setproj!(dataset, "EPSG:4326")

        AG.setunittype!(band, "day of year")
        AG.setnodatavalue!(band, -999)
    end
end

# ============================================================
# load_boundary — load boundary polygon from GeoJSON/Shapefile
# Returns Vector of polygons, each polygon is Vector of (lon, lat) tuples
# CRS is automatically converted to WGS84 (EPSG:4326) if needed
# ============================================================
# load_boundary — load boundary polygon from GeoJSON file
# Returns Vector of polygons, each polygon is Vector of (lon, lat) tuples
# ============================================================
function load_boundary(filepath::String)::Vector{Vector{Tuple{Float64,Float64}}}
    polygons = Vector{Vector{Tuple{Float64,Float64}}}()

    # Read and parse GeoJSON manually (more reliable than GDAL for simple GeoJSON)
    json_str = read(filepath, String)
    json = JSON.parse(json_str)

    features = get(json, "features", [])
    for feature in features
        geom = get(feature, "geometry", nothing)
        if geom === nothing
            continue
        end

        geom_type = get(geom, "type", "")
        coords = get(geom, "coordinates", [])

        if geom_type == "Polygon"
            # coords[1] is the outer ring
            ring = coords[1]
            polygon = [(c[1], c[2]) for c in ring]
            push!(polygons, polygon)
        elseif geom_type == "MultiPolygon"
            # coords is array of polygons, each polygon is array of rings
            for poly_coords in coords
                ring = poly_coords[1]  # outer ring
                polygon = [(c[1], c[2]) for c in ring]
                push!(polygons, polygon)
            end
        end
    end

    return polygons
end

# ============================================================
# point_in_polygon — ray-casting algorithm to test if a point is inside a polygon
# ============================================================
function point_in_polygon(lon::Float64, lat::Float64, polygon::Vector{Tuple{Float64,Float64}})::Bool
    n = length(polygon)
    if n < 3
        return false
    end

    inside = false
    j = n
    for i in 1:n
        xi, yi = polygon[i]
        xj, yj = polygon[j]

        # Check if point is on the edge
        if ((yi > lat) != (yj > lat)) &&
           (lon < (xj - xi) * (lat - yi) / (yj - yi) + xi)
            inside = !inside
        end
        j = i
    end

    return inside
end

# ============================================================
# point_near_polygon — check if a point is within buffer distance of polygon
# Uses simple bounding box + distance to line segments
# ============================================================
function point_near_polygon(lon::Float64, lat::Float64,
                           polygon::Vector{Tuple{Float64,Float64}},
                           buffer_deg::Float64)::Bool
    n = length(polygon)
    if n < 3
        return false
    end

    # First check if inside (with original polygon)
    if point_in_polygon(lon, lat, polygon)
        return true
    end

    # Check distance to each edge
    for i in 1:n
        j = i % n + 1
        x1, y1 = polygon[i]
        x2, y2 = polygon[j]

        dist = _point_to_segment_distance(lon, lat, x1, y1, x2, y2)
        if dist <= buffer_deg
            return true
        end
    end

    return false
end

# Helper: perpendicular distance from point to line segment
function _point_to_segment_distance(px::Float64, py::Float64,
                                    x1::Float64, y1::Float64,
                                    x2::Float64, y2::Float64)::Float64
    dx = x2 - x1
    dy = y2 - y1

    # If segment is a point
    if dx == 0 && dy == 0
        return sqrt((px - x1)^2 + (py - y1)^2)
    end

    # Project point onto line, clamped to segment
    t = max(0, min(1, ((px - x1) * dx + (py - y1) * dy) / (dx^2 + dy^2)))

    proj_x = x1 + t * dx
    proj_y = y1 + t * dy

    return sqrt((px - proj_x)^2 + (py - proj_y)^2)
end

# ============================================================
# create_boundary_mask — create boolean mask for pixels within/contacting boundary
# Returns n_lon x n_lat matrix: true = pixel is within/contacting boundary
# ============================================================
function create_boundary_mask(lons::Vector{Float64}, lats::Vector{Float64},
                              boundary_file::String; buffer_deg::Float64=0.0)::Matrix{Bool}
    n_lon = length(lons)
    n_lat = length(lats)

    # Load boundary polygons
    polygons = load_boundary(boundary_file)

    if isempty(polygons)
        @warn "No polygons found in boundary file: $boundary_file"
        return fill(true, n_lon, n_lat)
    end

    mask = Matrix{Bool}(undef, n_lon, n_lat)

    for i_lon in 1:n_lon
        for i_lat in 1:n_lat
            lon = lons[i_lon]
            lat = lats[i_lat]

            # Check if point is in or near any polygon
            in_boundary = false
            for poly in polygons
                if buffer_deg > 0
                    if point_near_polygon(lon, lat, poly, buffer_deg)
                        in_boundary = true
                        break
                    end
                else
                    if point_in_polygon(lon, lat, poly)
                        in_boundary = true
                        break
                    end
                end
            end

            mask[i_lon, i_lat] = in_boundary
        end
    end

    return mask
end
