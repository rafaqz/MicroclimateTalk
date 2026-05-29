# Read cached NetCDF outputs and render every figure used in the slides.
# Each demo's `cached_solve` will write `output/*.nc` on first solve; this
# script touches neither MicroclimateMapper nor the solver.

include("shared.jl")

# ---------------------------------------------------------------------------
# Shared Makie helpers
# ---------------------------------------------------------------------------

# `raster` slices are 3-D (X, Y, Ti) already in plot units (°C, cm, ...).
function plot_panels(raster, panel_indices, panel_labels,
                     variable_label, filename;
                     cmap = Reverse(:RdYlBu), unit_label = "°C")
    valid = filter(!isnan, vec(parent(raster)))
    isempty(valid) && return
    clims = (minimum(valid), maximum(valid))
    figure = Figure(size = (1400, 900))
    Label(figure[0, 1:3], "$variable_label — Mont Aigoual";
        fontsize = 16, font = :bold)
    for k in eachindex(panel_indices)
        row, col = fldmod1(k, 3)
        ax = Axis(figure[row, col];
            title = panel_labels[k], titlesize = 11,
            xlabel = "Longitude", ylabel = "Latitude",
            aspect = DataAspect())
        plot!(ax, view(raster; Ti = panel_indices[k]);
            colormap = cmap, colorrange = clims)
    end
    Colorbar(figure[1:2, 4]; colormap = cmap, colorrange = clims, label = unit_label)
    save(filename, figure)
    println("  Saved $filename")
end

function animate_diurnal(raster, hour_labels, variable_label, filename;
                         cmap = Reverse(:RdYlBu),
                         unit_label = "°C", framerate = 4)
    valid = filter(!isnan, vec(parent(raster)))
    isempty(valid) && return
    clims = (minimum(valid), maximum(valid))
    n_frames = size(raster, Ti)
    figure = Figure(size = (700, 600))
    title_obs = Observable("$variable_label\n$(hour_labels[1]) — Mont Aigoual")
    ax = Axis(figure[1, 1];
        title = title_obs, titlesize = 11,
        xlabel = "Longitude", ylabel = "Latitude",
        aspect = DataAspect())
    frame_obs = Observable(view(raster; Ti = 1))
    plot!(ax, frame_obs; colormap = cmap, colorrange = clims)
    Colorbar(figure[1, 2]; colormap = cmap, colorrange = clims, label = unit_label)
    record(figure, filename, 1:n_frames; framerate) do k
        frame_obs[] = view(raster; Ti = k)
        title_obs[] = "$variable_label\n$(hour_labels[k]) — Mont Aigoual"
    end
    println("  Saved $filename")
end

# ---------------------------------------------------------------------------
# Vector demo plots (Cévennes points, winter 2010)
# ---------------------------------------------------------------------------

println("Vector demo plots...")
# Vector output uses a MergedLookup `point` dim that NCDatasets can't store,
# so we re-run vector_demo.jl rather than caching. It's only ~14 s.
include("vector_demo.jl")
vector_output = strip_to_canonical(output)
point_labels  = ["Valley (900 m)", "South slope (1200 m)", "Summit (1567 m)"]
vector_year   = Dates.year(first(lookup(vector_output, Ti)))

ti_vector       = lookup(vector_output, Ti)
noon_indices    = findall(t -> hour(t) == 12, ti_vector)
winter_indices  = filter(i -> month(ti_vector[i]) in (1, 2, 3, 11, 12), noon_indices)
winter_dates    = collect(ti_vector[winter_indices])

snow_figure = Figure(size = (1000, 500))
snow_axis = Axis(snow_figure[1, 1];
    title  = "Snow depth — winter $vector_year",
    xlabel = "Date", ylabel = "Snow depth (cm)",
)
for (i, name) in enumerate(point_labels)
    depths_cm = view(vector_output.snow_depth; point = i, Ti = winter_indices)
    lines!(snow_axis, winter_dates, parent(depths_cm); label = name, linewidth = 2)
end
axislegend(snow_axis; position = :lt)
save(joinpath(FIG_DIR, "point_snow_depth.png"), snow_figure)
println("  Saved $(joinpath(FIG_DIR, "point_snow_depth.png"))")

target_date  = Date(vector_year, 1, 15)
day_indices  = findall(t -> Date(t) == target_date, ti_vector)
hours_of_day = hour.(ti_vector[day_indices])

extract_day(layer; kw...) = parent(view(layer; Ti = day_indices, kw...))

variable_panels = [
    ("Air, 2 m",     i -> extract_day(vector_output.air_temperature;  point = i, height = 2)),
    ("Air, 1 cm",    i -> extract_day(vector_output.air_temperature;  point = i, height = 1)),
    ("Soil surface", i -> extract_day(vector_output.soil_temperature; point = i, depth  = 1)),
    ("Soil, 10 cm",  i -> extract_day(vector_output.soil_temperature; point = i, depth  = 4)),
]
colors = [:steelblue, :darkorange, :firebrick]

temperatures_figure = Figure(size = (1200, 800))
Label(temperatures_figure[0, 1:2], "Mont Aigoual — $target_date";
    fontsize = 16, font = :bold)
for (k, (label, extract)) in enumerate(variable_panels)
    row, col = fldmod1(k, 2)
    ax = Axis(temperatures_figure[row, col];
        title = label, xlabel = "Hour of day", ylabel = "°C")
    for (i, name) in enumerate(point_labels)
        lines!(ax, hours_of_day, extract(i);
            label = name, color = colors[i], linewidth = 2)
    end
    k == 1 && axislegend(ax; position = :lt)
end
save(joinpath(FIG_DIR, "point_temperatures.png"), temperatures_figure)
println("  Saved $(joinpath(FIG_DIR, "point_temperatures.png"))")

# ---------------------------------------------------------------------------
# Raster demo plots (small high-res Mont Aigoual box, July 2000)
# ---------------------------------------------------------------------------

println("Raster demo plots...")
raster_output = RasterStack(joinpath(OUTPUT_DIR, "raster_summer.nc"))

ti_raster = lookup(raster_output, Ti)
july = findall(d -> month(d) == 7, ti_raster)
hour_labels_summer = [@sprintf("%02d:00", hour(t)) for t in ti_raster[july]]

slice_raster(layer; kw...) = view(layer; Ti = july, kw...)

air_temperature_2m       = slice_raster(raster_output.air_temperature;  height = 2)
air_temperature_1cm      = slice_raster(raster_output.air_temperature;  height = 1)
soil_temperature_surface = slice_raster(raster_output.soil_temperature; depth  = 1)
soil_temperature_10cm    = slice_raster(raster_output.soil_temperature; depth  = 4)

panel_indices = [1, 7, 10, 13, 16, 19]
panel_labels  = ["Midnight", "Dawn", "Mid-morning", "Midday", "Mid-afternoon", "Dusk"]

plot_panels(air_temperature_2m,       panel_indices, panel_labels,
    "Air temperature at 2 m (°C)",    joinpath(FIG_DIR, "grid_air_temperature_2m.png"))
plot_panels(air_temperature_1cm,      panel_indices, panel_labels,
    "Air temperature at 1 cm (°C)",   joinpath(FIG_DIR, "grid_air_temperature_1cm.png"))
plot_panels(soil_temperature_surface, panel_indices, panel_labels,
    "Soil surface temperature (°C)",  joinpath(FIG_DIR, "grid_soil_temperature_surface.png"))
plot_panels(soil_temperature_10cm,    panel_indices, panel_labels,
    "Soil temperature at 10 cm (°C)", joinpath(FIG_DIR, "grid_soil_temperature_10cm.png"))

animate_diurnal(air_temperature_1cm,      hour_labels_summer,
    "Air temperature at 1 cm (°C)",  joinpath(FIG_DIR, "grid_air_temperature_1cm.gif"))
animate_diurnal(soil_temperature_surface, hour_labels_summer,
    "Soil surface temperature (°C)", joinpath(FIG_DIR, "grid_soil_temperature_surface.gif"))

# ---------------------------------------------------------------------------
# Snow demo plots
# ---------------------------------------------------------------------------

println("Spring snowmelt plots...")
spring_snow_output = RasterStack(joinpath(OUTPUT_DIR, "snow_spring.nc"))
ti_spring          = lookup(spring_snow_output, Ti)
april              = findall(d -> month(d) == 4, ti_spring)
spring_hour_labels = [@sprintf("%02d:00", hour(t)) for t in ti_spring[april]]
spring_snow_depth  = view(spring_snow_output.snow_depth; Ti = april)

plot_panels(spring_snow_depth, panel_indices, panel_labels,
    "Snow depth (cm) — April diurnal cycle",
    joinpath(FIG_DIR, "grid_snow_depth.png");
    cmap = Reverse(:ice), unit_label = "cm")
animate_diurnal(spring_snow_depth, spring_hour_labels,
    "Snow depth (cm) — April diurnal cycle",
    joinpath(FIG_DIR, "grid_snow_depth.gif");
    cmap = Reverse(:ice), unit_label = "cm")

println("Year-long snow plots...")
yearly_snow_output = RasterStack(joinpath(OUTPUT_DIR, "snow_yearly.nc"))
ti_yearly          = lookup(yearly_snow_output, Ti)
yearly_year        = year(first(ti_yearly))

weekly_noon   = findall(t -> hour(t) == 12 && (dayofyear(t) - 1) % 7 == 0, ti_yearly)
weekly_labels = [Dates.format(ti_yearly[k], "u dd") for k in weekly_noon]
snow_weekly   = view(yearly_snow_output.snow_depth; Ti = weekly_noon)

animate_diurnal(snow_weekly, weekly_labels,
    "Snow depth (cm) — seasonal cycle",
    joinpath(FIG_DIR, "grid_snow_depth_yearly.gif");
    cmap = Reverse(:ice), unit_label = "cm", framerate = 8)

monthly_idx  = [findfirst(t -> month(t) == m && day(t) == 15 && hour(t) == 12, ti_yearly) for m in 1:12]
month_labels = [Dates.format(ti_yearly[k], "U") for k in monthly_idx]
snow_monthly = view(yearly_snow_output.snow_depth; Ti = monthly_idx)

let valid = filter(!isnan, vec(parent(snow_monthly)))
    if !isempty(valid)
        clims = (minimum(valid), maximum(valid))
        cmap  = Reverse(:ice)
        figure = Figure(size = (1600, 1100))
        Label(figure[0, 1:4],
            "Snow depth (cm) — seasonal cycle at noon, Mont Aigoual ($yearly_year)";
            fontsize = 16, font = :bold)
        for k in eachindex(month_labels)
            row, col = fldmod1(k, 4)
            ax = Axis(figure[row, col];
                title = month_labels[k], titlesize = 11,
                xlabel = "Longitude", ylabel = "Latitude",
                aspect = DataAspect())
            plot!(ax, view(snow_monthly; Ti = k);
                colormap = cmap, colorrange = clims)
        end
        Colorbar(figure[1:3, 5]; colormap = cmap, colorrange = clims, label = "cm")
        filename = joinpath(FIG_DIR, "grid_snow_depth_yearly.png")
        save(filename, figure)
        println("  Saved $filename")
    end
end
