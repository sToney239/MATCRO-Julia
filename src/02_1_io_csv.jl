# IO_CSV — CSV-based I/O for single-point MATCRO simulation
# Used when config.input_format = "csv"
# Companion: io.jl (core structs), io_netcdf.jl (NetCDF I/O)

using CSV
using DataFrames

# ============================================================
# CSV forcing format:
#   year,doy,tmax,tmin,radiation,precip,humidity,wind,pressure[,ozone]
#   2000,1,271.3,257.8,86.6,0.0,0.005992,4.65,101384.1,26.5
# Units: tmax/tmin [K], radiation [W/m²], precip [kg/m²/s],
#        humidity [kg/kg], wind [m/s], pressure [Pa], ozone [ppb]
# ============================================================

# ============================================================
# read_forcing_csv — load all daily forcing from CSV into a
# Dict{year => Dict{doy => DailyForcing}}
# ============================================================
function read_forcing_csv(csv_path::String)::Dict{Int,Dict{Int,DailyForcing}}
    df = CSV.read(csv_path, DataFrame)

    # Normalize column names to symbols
    rename!(df, Symbol.(lowercase.(string.(names(df)))))

    result = Dict{Int,Dict{Int,DailyForcing}}()

    for row in eachrow(df)
        yr = Int(row[:year])
        dy = Int(row[:doy])

        forcing = DailyForcing(;
            doy       = dy,
            tmax      = Float64(row[:tmax]),
            tmin      = Float64(row[:tmin]),
            radiation = Float64(row[:radiation]),
            precip    = Float64(row[:precip]),
            humidity  = Float64(row[:humidity]),
            wind      = Float64(row[:wind]),
            pressure  = Float64(row[:pressure]),
            ozone     = hasproperty(row, :ozone) ? Float64(row[:ozone]) : 0.0,
        )

        if !haskey(result, yr)
            result[yr] = Dict{Int,DailyForcing}()
        end
        result[yr][dy] = forcing
    end

    return result
end

# ============================================================
# get_forcing — retrieve DailyForcing for a specific year/doy
# Returns nothing if data is missing
# ============================================================
function get_forcing(forcing_data::Dict{Int,Dict{Int,DailyForcing}},
                     year::Int, doy::Int)::Union{DailyForcing,Nothing}
    if haskey(forcing_data, year) && haskey(forcing_data[year], doy)
        return forcing_data[year][doy]
    end
    return nothing
end

# ============================================================
# write_output_csv — write simulation results to CSV
# ============================================================
function write_output_csv(results::Vector{NamedTuple}, output_path::String)
    if isempty(results)
        @warn "No results to write"
        return
    end

    df = DataFrame(results)
    CSV.write(output_path, df)
    println("Output written to ", output_path)
end

# ============================================================
# write_daily_csv — write detailed daily state variables to CSV
# ============================================================
function write_daily_csv(daily_records::Vector{NamedTuple}, output_path::String)
    if isempty(daily_records)
        return
    end

    df = DataFrame(daily_records)
    CSV.write(output_path, df)
    println("Daily output written to ", output_path)
end
