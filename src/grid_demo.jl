# Gridded microclimate over Mont Aigoual (Cévennes, southern France).
#
# Pipeline:
#   1. Download SRTM DEM and reproject to UTM
#   2. Compute slope, aspect, horizon angles (Geomorphometry.jl)
#   3. Download TerraClimate weather for year 2000 at the centre pixel
#   4. Per-pixel simulate_microclimate with lapse-rate-corrected weather
#   5. Plot soil-surface T, air T at 1 cm and 2 m, soil T at 10 cm
#      at 6 representative hours; also save a 24-frame GIF.
#
# Threading: start Julia with `julia --threads auto`.
#
# Writes into figures/:
#   grid_Tair2m.png    grid_Tair1cm.png    grid_Tsoil0.png    grid_Tsoil10cm.png
#   grid_Tair1cm.gif   grid_Tsoil0.gif

using BiophysicalGrids
using RasterDataSources
using Rasters, ArchGDAL
using SolarRadiation
using FluidProperties
using FluidProperties: VPLookupTable
using Unitful
using Statistics: median
using Printf
using Plots
import Plots: heatmap, plot, savefig

const FIG_DIR = joinpath(@__DIR__, "..", "figures")

# Mont Aigoual, Cévennes — 1567 m peak, dramatic N/S aspect contrast,
# 30 km NW of Montpellier.
center_lon = 3.581    # °E
center_lat = 44.122   # °N
extent_lat = 0.060    # ~70 SRTM pixels N–S — keep CI build fast
extent_lon = 0.085    # ~70 SRTM pixels E–W

year             = 2000
july             = 7
n_horizon_angles = 24
vp_method        = VPLookupTable()

depths  = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm"
heights = [0.01, 2.0]u"m"

snapshot_hours = collect(0:23)
snapshot_steps = snapshot_hours .+ 1
hour_labels    = [@sprintf("%02d:00", h) for h in snapshot_hours]
nhours         = length(snapshot_hours)

panel_ks     = [1, 7, 10, 13, 16, 19]
panel_labels = ["Midnight", "Dawn", "Mid-morning", "Midday", "Mid-afternoon", "Dusk"]

println("Downloading SRTM DEM and reprojecting to UTM...")
(; utm_dem, x_coords_utm, y_coords_utm, nx_utm, ny_utm, cs) =
    load_utm_dem(center_lon, center_lat, extent_lon, extent_lat)
println("  UTM grid: $(nx_utm) × $(ny_utm) px, " *
        "cell ≈ $(round(cs[1]; digits=1)) × $(round(cs[2]; digits=1)) m")

println("Computing terrain grids (slope, aspect, horizons)...")
(; dem_data, data_is_xy, y_descending,
   elevation_m, slope_deg, aspect_deg,
   latitude_deg, longitude_deg, pressure_r,
   horizons_u) = compute_terrain_grids(utm_dem, x_coords_utm, y_coords_utm;
                                       n_horizon_angles)

println("Obtaining TerraClimate weather for year $year...")
valid_elev    = filter(!isnan, vec(dem_data))
center_elev_u = median(valid_elev) * u"m"

weather = get_weather(TerraClimate, center_lon, center_lat;
    ystart    = year,
    elevation = center_elev_u,
)

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

weather_july = extract_month(weather, july)

center_elev_m = round(ustrip(u"m",  center_elev_u);                                          digits = 0)
tmin_july_C   = round(ustrip(u"°C", weather_july.environment_minmax.reference_temperature_min[1]); digits = 1)
tmax_july_C   = round(ustrip(u"°C", weather_july.environment_minmax.reference_temperature_max[1]); digits = 1)
tdeep_C       = round(ustrip(u"°C", weather_july.environment_daily.deep_soil_temperature[1]);     digits = 1)
println("  Centre elevation: $center_elev_m m")
println("  July Tmin $tmin_july_C °C,  Tmax $tmax_july_C °C,  deep soil T $tdeep_C °C")

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

solar_model = SolarProblem()

function lapse_correct_weather(ws, elev_diff; method = vp_method)
    mm = ws.environment_minmax
    ed = ws.environment_daily
    T_min_new = lapse_adjust_temperature(mm.reference_temperature_min, elev_diff, EnvironmentalLapseRate())
    T_max_new = lapse_adjust_temperature(mm.reference_temperature_max, elev_diff, EnvironmentalLapseRate())
    new_mm = MonthlyMinMaxEnvironment(;
        reference_temperature_min = T_min_new,
        reference_temperature_max = T_max_new,
        reference_wind_min        = mm.reference_wind_min,
        reference_wind_max        = mm.reference_wind_max,
        reference_humidity_min    = rh_at_temperature(mm.reference_humidity_min, mm.reference_temperature_max, T_max_new, method),
        reference_humidity_max    = rh_at_temperature(mm.reference_humidity_max, mm.reference_temperature_min, T_min_new, method),
        cloud_min                 = mm.cloud_min,
        cloud_max                 = mm.cloud_max,
        minima_times              = mm.minima_times,
        maxima_times              = mm.maxima_times,
    )
    new_ed = DailyTimeseries(;
        shade                 = ed.shade,
        soil_wetness          = ed.soil_wetness,
        surface_emissivity    = ed.surface_emissivity,
        cloud_emissivity      = ed.cloud_emissivity,
        rainfall              = ed.rainfall,
        deep_soil_temperature = lapse_adjust_temperature(ed.deep_soil_temperature, elev_diff, EnvironmentalLapseRate()),
        leaf_area_index       = ed.leaf_area_index,
    )
    return merge(ws, (; environment_minmax = new_mm, environment_daily = new_ed))
end

println("Running per-pixel microclimate ($(nx_utm) × $(ny_utm), " *
        "$(Threads.nthreads()) thread(s))...")

wp_grid    = Array{Any}(undef, ny_utm, nx_utm)
Tmean_grid = fill(NaN * u"K", ny_utm, nx_utm)
Tdeep_grid = fill(NaN * u"K", ny_utm, nx_utm)

for I in CartesianIndices((ny_utm, nx_utm))
    i, j   = I[1], I[2]
    ri, rj = data_is_xy ? (j, i) : (i, j)
    elev   = elevation_m[ri, rj]
    ismissing(elev) && continue
    wp = lapse_correct_weather(weather_july, elev - center_elev_u)
    wp_grid[i, j]    = wp
    Tmean_grid[i, j] = (wp.environment_minmax.reference_temperature_min[1] +
                        wp.environment_minmax.reference_temperature_max[1]) / 2
    Tdeep_grid[i, j] = wp.environment_daily.deep_soil_temperature[1]
end

T_soil0  = fill(NaN, ny_utm, nx_utm, nhours)
T_air1   = fill(NaN, ny_utm, nx_utm, nhours)
T_air2   = fill(NaN, ny_utm, nx_utm, nhours)
T_soil10 = fill(NaN, ny_utm, nx_utm, nhours)

n_total = nx_utm * ny_utm
n_done  = Threads.Atomic{Int}(0)
T_converged = fill!(Matrix{Any}(undef, ny_utm, nx_utm), nothing)

@time for i in 1:ny_utm
    Threads.@threads :static for j in 1:nx_utm
        ri, rj = data_is_xy ? (j, i) : (i, j)

        elev = elevation_m[ri, rj]
        lat  = latitude_deg[ri, rj]
        lon  = longitude_deg[ri, rj]
        slp  = slope_deg[ri,    rj]
        asp  = aspect_deg[ri,   rj]
        pres = pressure_r[ri,   rj]

        if ismissing(elev) || ismissing(lat) || ismissing(lon) ||
           ismissing(slp)  || ismissing(asp) || ismissing(pres)
            continue
        end

        wp         = wp_grid[i, j]
        Tmean_july = Tmean_grid[i, j]
        Tdeep      = Tdeep_grid[i, j]

        T_init = if i > 1 && !isnothing(T_converged[i - 1, j])
            T_converged[i - 1, j]
        else
            T_tmp      = Vector{typeof(1.0u"K")}(undef, 10)
            T_tmp[1:8] .= Tmean_july
            T_tmp[9]    = (Tmean_july + Tdeep) / 2
            T_tmp[10]   = Tdeep
            T_tmp
        end

        site = Site(;
            latitude             = lat,
            longitude            = lon,
            elevation            = elev,
            slope                = slp,
            aspect               = asp,
            horizon_angles       = @view(horizons_u[i, j, :]),
            sky_view_fraction    = 1.0,
            albedo               = 0.15,
            roughness_height     = 0.004u"m",
            atmospheric_pressure = pres,
        )

        result = simulate_microclimate(
            site, soil_thermal, soil_hydraulics, wp;
            depths, heights, solar_model,
            initial_soil_temperature = T_init,
            vapour_pressure_equation = vp_method,
            iterate_day = 5,
        )

        T_converged[i, j] = collect(result.soil_temperature[1, :])

        @inbounds for (k, s) in enumerate(snapshot_steps)
            T_soil0[i,  j, k] = ustrip(u"°C", result.soil_temperature[s, 1])
            T_air1[i,   j, k] = ustrip(u"°C", result.profile[s].air_temperature[1])
            T_air2[i,   j, k] = ustrip(u"°C", result.profile[s].air_temperature[2])
            T_soil10[i, j, k] = ustrip(u"°C", result.soil_temperature[s, 4])
        end

        d = Threads.atomic_add!(n_done, 1)
        if Threads.threadid() == 1 && (d % nx_utm == 0 || d == n_total - 1)
            pct = round(100 * (d + 1) / n_total; digits = 1)
            print("  row $i/$ny_utm  ($(d+1)/$n_total, $pct%)   \r")
        end
    end
end
println("\nSimulation complete.")

y_plt = ascending_y(y_coords_utm, zeros(ny_utm, 1))[1]
common_kw = (; aspect_ratio = :equal, xlabel = "Easting (m)", ylabel = "Northing (m)")
mkpath(FIG_DIR)

function plot_variable(data4d, var_label, fname)
    all_vals = filter(!isnan, vec(data4d))
    isempty(all_vals) && return
    clims = (minimum(all_vals), maximum(all_vals))
    nframes = size(data4d, 3)
    ks = nframes <= 6 ? (1:nframes) : panel_ks
    ls = nframes <= 6 ? hour_labels[1:nframes] : panel_labels

    panels = [heatmap(x_coords_utm, y_plt, ascending_y(y_coords_utm, data4d[:, :, ks[n]])[2];
        color = cgrad(:RdYlBu, rev = true), clims = clims,
        title = ls[n], colorbar_title = "°C",
        titlefontsize = 9, common_kw...) for n in eachindex(ks)]

    fig = plot(panels...; layout = (2, 3), size = (1400, 900),
        left_margin = 5Plots.mm,
        plot_title = "$var_label — Mont Aigoual, July $year")
    savefig(fig, fname)
    println("  Saved $fname")
end

function animate_variable(data4d, var_label, fname; framerate = 4)
    all_vals = filter(!isnan, vec(data4d))
    isempty(all_vals) && return
    clims = (minimum(all_vals), maximum(all_vals))
    nframes = size(data4d, 3)
    labels_here = nframes == length(snapshot_hours) ?
        [@sprintf("%02d:00", snapshot_hours[k]) for k in 1:nframes] :
        [@sprintf("frame %d", k) for k in 1:nframes]

    anim = @animate for k in 1:nframes
        heatmap(x_coords_utm, y_plt, ascending_y(y_coords_utm, data4d[:, :, k])[2];
            color = cgrad(:RdYlBu, rev = true), clims = clims,
            title = "$var_label\n$(labels_here[k]) — Mont Aigoual, July $year",
            xlabel = "Easting (m)", ylabel = "Northing (m)",
            colorbar_title = "°C", aspect_ratio = :equal,
            titlefontsize = 9, size = (700, 600),
            left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)
    end
    gif(anim, fname; fps = framerate)
    println("  Saved $fname")
end

println("Plotting panels...")
plot_variable(T_air2,   "Air temperature at 2 m (°C)",    joinpath(FIG_DIR, "grid_Tair2m.png"))
plot_variable(T_air1,   "Air temperature at 1 cm (°C)",   joinpath(FIG_DIR, "grid_Tair1cm.png"))
plot_variable(T_soil0,  "Soil surface temperature (°C)",  joinpath(FIG_DIR, "grid_Tsoil0.png"))
plot_variable(T_soil10, "Soil temperature at 10 cm (°C)", joinpath(FIG_DIR, "grid_Tsoil10cm.png"))

println("Animating (1 cm air, soil surface)...")
animate_variable(T_air1,  "Air temperature at 1 cm (°C)",  joinpath(FIG_DIR, "grid_Tair1cm.gif"))
animate_variable(T_soil0, "Soil surface temperature (°C)", joinpath(FIG_DIR, "grid_Tsoil0.gif"))

println("\nDone. $(nx_utm)×$(ny_utm) px grid, July $year.")
