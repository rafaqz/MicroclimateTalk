# Point microclimate at Montpellier, July 2000.
#
# Driven by TerraClimate monthly weather; demonstrates the two-step API:
#   1. get_weather(TerraClimate, ...) → environment structs
#   2. simulate_microclimate(site, soil_thermal, soil_hydraulics, weather)
#
# Writes:
#   figures/point_temperatures.png — 2 m air, 1 cm air, soil surface, soil 10 cm

using BiophysicalGrids
using Microclimate
using FluidProperties
using FluidProperties: GoffGratch
using Unitful
using Printf
using Plots

const FIG_DIR = joinpath(@__DIR__, "..", "figures")

# Montpellier, southern France
lon, lat = 3.8772, 43.6108
elevation = 27.0u"m"

depths  = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm"
heights = [0.01, 2.0]u"m"
year    = 2000
month   = 7

println("Downloading TerraClimate weather for Montpellier, $year...")
weather = get_weather(TerraClimate, lon, lat;
    ystart = year,
    elevation,
    vapour_pressure_method = GoffGratch(),
    lapse_rate_type        = EnvironmentalLapseRate(),
)

# Slice to July (24 hourly entries per month).
function extract_month(ws, m)
    mm = ws.environment_minmax
    ed = ws.environment_daily
    eh = ws.environment_hourly
    new_mm = MonthlyMinMaxEnvironment(;
        reference_temperature_min = mm.reference_temperature_min[[m]],
        reference_temperature_max = mm.reference_temperature_max[[m]],
        reference_wind_min        = mm.reference_wind_min[[m]],
        reference_wind_max        = mm.reference_wind_max[[m]],
        reference_humidity_min    = mm.reference_humidity_min[[m]],
        reference_humidity_max    = mm.reference_humidity_max[[m]],
        cloud_min                 = mm.cloud_min[[m]],
        cloud_max                 = mm.cloud_max[[m]],
        minima_times              = mm.minima_times,
        maxima_times              = mm.maxima_times,
    )
    new_ed = DailyTimeseries(;
        shade                 = ed.shade[[m]],
        soil_wetness          = ed.soil_wetness[[m]],
        surface_emissivity    = ed.surface_emissivity[[m]],
        cloud_emissivity      = ed.cloud_emissivity[[m]],
        rainfall              = ed.rainfall[[m]],
        deep_soil_temperature = ed.deep_soil_temperature[[m]],
        leaf_area_index       = ed.leaf_area_index[[m]],
    )
    h_range = ((m - 1) * 24 + 1):(m * 24)
    new_eh = HourlyTimeseries(;
        pressure              = isnothing(eh.pressure)              ? nothing : eh.pressure[h_range],
        reference_temperature = isnothing(eh.reference_temperature) ? nothing : eh.reference_temperature[h_range],
        reference_humidity    = isnothing(eh.reference_humidity)    ? nothing : eh.reference_humidity[h_range],
        reference_wind_speed  = isnothing(eh.reference_wind_speed)  ? nothing : eh.reference_wind_speed[h_range],
        global_radiation      = isnothing(eh.global_radiation)      ? nothing : eh.global_radiation[h_range],
        longwave_radiation    = isnothing(eh.longwave_radiation)    ? nothing : eh.longwave_radiation[h_range],
        cloud_cover           = isnothing(eh.cloud_cover)           ? nothing : eh.cloud_cover[h_range],
        rainfall              = isnothing(eh.rainfall)              ? nothing : eh.rainfall[h_range],
        zenith_angle          = isnothing(eh.zenith_angle)          ? nothing : eh.zenith_angle[h_range],
    )
    return merge(ws, (;
        environment_minmax    = new_mm,
        environment_daily     = new_ed,
        environment_hourly    = new_eh,
        days                  = [ws.days[m]],
        soil_moisture_monthly = ws.soil_moisture_monthly[[m]],
    ))
end

weather_july = extract_month(weather, month)

site = Site(;
    latitude             = lat * u"°",
    longitude            = lon * u"°",
    elevation,
    slope                = 0.0u"°",
    aspect               = 0.0u"°",
    horizon_angles       = fill(0.0u"°", 24),
    sky_view_fraction    = 1.0,
    albedo               = 0.15,
    roughness_height     = 0.004u"m",
    atmospheric_pressure = atmospheric_pressure(elevation),
)

soil_thermal = CampbelldeVriesSoilThermal(;
    de_vries_shape_factor = 0.1,
    mineral_conductivity  = 2.5u"W/m/K",
    mineral_heat_capacity = 870.0u"J/kg/K",
    recirculation_power   = 4.0,
    return_flow_threshold = 0.162,
)

soil_hydraulics = example_soil_hydraulics(depths;
    bulk_density    = 1.3u"Mg/m^3",
    mineral_density = 2.56u"Mg/m^3",
    root_density    = fill(0.0, length(depths))u"m/m^3",
)

println("Running microclimate simulation...")
result = simulate_microclimate(
    site, soil_thermal, soil_hydraulics, weather_july;
    depths,
    heights,
    vapour_pressure_equation = GoffGratch(),
)

hours = 0:23
air_temperature_matrix = hcat([p.air_temperature for p in result.profile]...)'

T_air2m  = ustrip.(u"°C", u"K".(air_temperature_matrix[:, 2]))
T_air1cm = ustrip.(u"°C", u"K".(air_temperature_matrix[:, 1]))
T_soil0  = ustrip.(u"°C", u"K".(result.soil_temperature[:, 1]))
T_soil10 = ustrip.(u"°C", u"K".(result.soil_temperature[:, 4]))

p = plot(hours, T_air2m;
    label = "Air, 2 m",  color = :blue,  lw = 2,
    xlabel = "Hour of day",
    ylabel = "Temperature (°C)",
    title  = "Montpellier — representative July day",
    legend = :topleft,
    size   = (1000, 600),
)
plot!(p, hours, T_air1cm; label = "Air, 1 cm", color = :orange, lw = 2)
plot!(p, hours, T_soil0;  label = "Soil surface", color = :red,    lw = 2)
plot!(p, hours, T_soil10; label = "Soil, 10 cm",  color = :brown,  lw = 2)

mkpath(FIG_DIR)
savefig(p, joinpath(FIG_DIR, "point_temperatures.png"))
println("Wrote $(joinpath(FIG_DIR, "point_temperatures.png"))")
