# TINTERP - Meteorological data time interpolation
# Interpolates daily weather data to hourly values
# Temperature: Goudriaan (1994) sine/exponential method
# Radiation: proportional to sin(solar elevation)
# Wind: logarithmic profile adjustment



const HOURTMX = 1.5         # hour of maximum temperature
const NocT = 4.0            # nocturnal time coefficient

# ============ DAYL: day length [hours] ============
function day_length(doy::Int, lat::Float64)::Float64
    # Solar declination
    del = -asin(sin(23.45 * 2π / 360) * cos(2π * (Float64(doy) + 10.0) / 365.0))
    latrad = 2π * lat / 360.0

    dd = sin(latrad) * sin(del) / (cos(latrad) * cos(del))

    if dd > 1.0
        return 24.0
    elseif dd < -1.0
        return 0.0
    else
        return 12.0 * (1.0 + (2.0 / π) * asin(dd))
    end
end

# ============ SINB: sin(solar elevation) ============
# Goudriaan and van Laar (1994)
function sin_solar_elevation(doy::Int, hour::Float64, lat::Float64)::Float64
    # Solar declination
    del = -asin(sin(23.45 * 2π / 360) * cos(2π * (Float64(doy) + 10.0) / 365.0))
    # Hour angle
    h = π * (hour - 12.0) / 12.0
    latrad = 2π * lat / 360.0

    return max(0.0, sin(latrad) * sin(del) + cos(latrad) * cos(del) * cos(h))
end

# ============ INTERPOLATE_TIME: main interpolation function ============
function interpolate_time(; doy::Int, prev_doy::Int, next_doy::Int, hour::Float64,
                 lat::Float64, Δt::Int,
                 # Daily inputs
                 tmax_prev::Float64,    # INTMXp: yesterday's max temp [K]
                 tmax::Float64,         # INTMX: today's max temp [K]
                 tmax_next::Float64,    # INTMXn: tomorrow's max temp [K]
                 tmin_prev::Float64,    # INTMNp: yesterday's min temp [K]
                 tmin::Float64,         # INTMN: today's min temp [K]
                 tmin_next::Float64,    # INTMNn: tomorrow's min temp [K]
                 radiation::Float64,    # INRSD: daily solar radiation [W/m²]
                 precip::Float64,       # INPRC: daily precipitation [mm/day or kg/m²/s]
                 humidity::Float64,     # INSHM: daily specific humidity [kg/kg]
                 wind::Float64,         # INWND: daily wind speed [m/s]
                 pressure::Float64,     # INPRS: daily air pressure [Pa]
                 ozone::Float64,        # INOZN: daily ozone
                 int_sinb::Float64,     # INTSINB: integrated sin(solar elevation) over day
                 wind_height::Float64   # WND_HGT: wind measurement height [m]
                 )::NamedTuple

    # Day length
    dl = day_length(doy, lat)
    dlp = day_length(prev_doy, lat)

    # Convert K → °C for Goudriaan temperature interpolation
    tmax_prev_c = tmax_prev - T_ice
    tmax_c      = tmax - T_ice
    tmax_next_c = tmax_next - T_ice
    tmin_prev_c = tmin_prev - T_ice
    tmin_c      = tmin - T_ice
    tmin_next_c = tmin_next - T_ice

    # Wind: logarithmic profile from measurement height to 2m
    wnd = wind * log(2.0 / 0.05) / log(wind_height / 0.05)

    # Pressure and precipitation: direct pass-through
    prs = pressure
    prc = precip
    shm = humidity

    # Radiation: distribute daily total proportional to sin(solar elevation)
    if int_sinb > 0.0
        rsd = radiation * 86400.0 * sin_solar_elevation(doy, hour, lat) / int_sinb
    else
        rsd = 0.0
    end

    # Temperature interpolation: Goudriaan (1994) P34
    # Three periods: before sunrise, daytime, after sunset
    sunrise = 12.0 - dl * 0.5
    sunset = 12.0 + dl * 0.5

    if hour < sunrise
        # Before sunrise: exponential decay from previous day's sunset
        tset = tmin_prev_c + (tmax_prev_c - tmin_prev_c) * sin(π * dlp / (dlp + 2.0 * HOURTMX))
        nt = 24.0 - dl  # night time [hours]
        tmp = (tmin_c - tset * exp(-nt / NocT) +
               (tset - tmin_c) * exp(-(hour + 24.0 - (12.0 + dlp * 0.5)) / NocT)) /
              (1.0 - exp(-nt / NocT))
    elseif hour < sunset
        # Daytime: sine curve from sunrise to sunset
        tmp = tmin_c + (tmax_c - tmin_c) * sin(π * (hour - 12.0 + dl * 0.5) / (dl + 2.0 * HOURTMX))
    else
        # After sunset: exponential decay
        tset = tmin_c + (tmax_c - tmin_c) * sin(π * dl / (dl + 2.0 * HOURTMX))
        nt = 24.0 - dl
        tmp = (tmin_next_c - tset * exp(-nt / NocT) +
               (tset - tmin_next_c) * exp(-(hour - (12.0 + dl * 0.5)) / NocT)) /
              (1.0 - exp(-nt / NocT))
    end

    # Convert back to Kelvin
    tmp += T_ice

    # Ozone: not implemented
    ozn = 0.0

    # Safety: minimum humidity
    if shm <= 0.0
        shm = 1.0e-5
    end

    return (temperature=tmp, precipitation=prc, radiation=rsd,
            humidity=shm, wind=wnd, pressure=prs, ozone=ozn)
end
