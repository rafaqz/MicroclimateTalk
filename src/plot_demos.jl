# Read cached NetCDF outputs and render every figure used in the slides.
# Each demo's `cached_solve` will write `output/*.nc` on first solve; this
# script touches neither MicroclimateMapper nor the solver.

include("shared.jl")

# ---------------------------------------------------------------------------
# Tufte-inspired Makie theme — matches the slide CSS:
#   - off-white #fffff8 background
#   - serif typography
#   - thin medium-grey axis spines, no top/right spines
#   - no bold titles, no gridlines
# ---------------------------------------------------------------------------

const TUFTE_BG    = RGBf(1.0, 1.0, 0.973)   # #fffff8
const TUFTE_INK   = RGBf(0.07, 0.07, 0.07)
const TUFTE_GREY  = RGBf(0.45, 0.45, 0.45)
const TUFTE_FAINT = RGBf(0.78, 0.76, 0.70)
const TUFTE_RED   = RGBf(0.63, 0.0,  0.0)

const TUFTE_SERIF = "TeX Gyre Pagella"   # Palatino-style serif, ships with Makie

set_theme!(Theme(
    backgroundcolor = TUFTE_BG,
    textcolor       = TUFTE_INK,
    fontsize        = 14,
    font            = TUFTE_SERIF,
    fonts           = (regular = TUFTE_SERIF, italic = TUFTE_SERIF,
                       bold = TUFTE_SERIF, bold_italic = TUFTE_SERIF),
    Figure = (; backgroundcolor = TUFTE_BG),
    Axis = (
        backgroundcolor    = TUFTE_BG,
        xgridvisible       = false,
        ygridvisible       = false,
        topspinevisible    = false,
        rightspinevisible  = false,
        leftspinecolor     = TUFTE_GREY,
        bottomspinecolor   = TUFTE_GREY,
        spinewidth         = 0.8,
        xtickcolor         = TUFTE_GREY,
        ytickcolor         = TUFTE_GREY,
        xticklabelcolor    = TUFTE_INK,
        yticklabelcolor    = TUFTE_INK,
        xticklabelsize     = 12,
        yticklabelsize     = 12,
        xlabelsize         = 13,
        ylabelsize         = 13,
        titlesize          = 14,
        titlefont          = TUFTE_SERIF,
        titlealign         = :left,
        titlegap           = 8,
        xlabelcolor        = TUFTE_INK,
        ylabelcolor        = TUFTE_INK,
    ),
    Axis3 = (
        xgridvisible    = false,
        ygridvisible    = false,
        zgridvisible    = false,
        xspinecolor_1   = TUFTE_FAINT,
        yspinecolor_1   = TUFTE_FAINT,
        zspinecolor_1   = TUFTE_FAINT,
        titlesize       = 14,
        titlefont       = TUFTE_SERIF,
        titlealign      = :left,
    ),
    Colorbar = (
        labelsize       = 12,
        ticklabelsize   = 11,
        spinewidth      = 0.5,
        tickcolor       = TUFTE_GREY,
        labelcolor      = TUFTE_INK,
        ticklabelcolor  = TUFTE_INK,
    ),
    Label = (; font = TUFTE_SERIF, color = TUFTE_INK),
    Lines = (; linewidth = 1.4),
    Legend = (
        framevisible    = false,
        labelsize       = 12,
        labelfont       = TUFTE_SERIF,
        backgroundcolor = :transparent,
        patchsize       = (18, 2),
    ),
))

# ---------------------------------------------------------------------------
# Shared Makie helpers
# ---------------------------------------------------------------------------

# `slice` is an X-Y Raster whose `missingval` marks cells to draw transparent.
# Returns a Float64 Matrix with those cells replaced by NaN so Makie's
# `nan_color = :transparent` punches them out of the surface colouring.
function drape_color(slice)
    mv = missingval(slice)
    A  = parent(slice)
    [_skip_cell(v, mv) ? NaN : Float64(v) for v in A]
end

# Returns a positive colour range from data that may be entirely missing/zero.
function drape_clims(raster)
    A = parent(raster)
    mv = missingval(raster)
    valid = Float64[]
    for v in A
        _skip_cell(v, mv) && continue
        push!(valid, Float64(v))
    end
    isempty(valid) && return nothing
    lo = minimum(valid); hi = maximum(valid)
    lo == hi ? (lo, lo + 1.0) : (lo, hi)
end

# A cell is "skip" if it's `missing`, equals a concrete missingval, or NaN.
# Guarding `v == mv` with `ismissing(mv)` avoids `missing > 0`-style propagation
# when the raster's missingval is itself `missing` (default for NCDF reads).
function _skip_cell(v, mv)
    ismissing(v) && return true
    !ismissing(mv) && v == mv && return true
    isnan(v) && return true
    return false
end

# Returns true if any cell in the X-Y slice has snow.
function _slice_has_snow(slice)
    A  = parent(slice)
    mv = missingval(slice)
    for v in A
        _skip_cell(v, mv) && continue
        v > 0 && return true
    end
    return false
end

# Scans `snow_depth` over the supplied `time_indices` (positions in its Ti dim),
# returns a sub-range covering the main snowpack: from `buffer_days` before
# first snow to `buffer_days` after the snowpack disappears for at least
# `min_bare_gap` consecutive samples. This excludes late sporadic summer/fall
# events that would otherwise stretch the animation indefinitely.
# Returns `nothing` if no snow is found.
function cold_season_window(snow_depth, time_indices;
                            buffer_days = 7, min_bare_gap = 14)
    has_snow  = [_slice_has_snow(view(snow_depth; Ti = k)) for k in time_indices]
    first_idx = findfirst(has_snow)
    first_idx === nothing && return nothing
    n        = length(has_snow)
    last_idx = first_idx
    bare_run = 0
    for k in first_idx:n
        if has_snow[k]
            last_idx = k
            bare_run = 0
        else
            bare_run += 1
            bare_run >= min_bare_gap && break
        end
    end
    lo = max(1, first_idx - buffer_days)
    hi = min(n, last_idx + buffer_days)
    time_indices[lo:hi]
end

# Pixel-thin z-offset so the snow overlay sits visibly above the DEM without
# z-fighting in CairoMakie's depth ordering.
const DRAPE_Z_OFFSET = 5.0

function add_terrain_base!(ax, xs, ys, zs)
    surface!(ax, xs, ys, zs;
        color = zs, colormap = :greys,
        colorrange = (minimum(zs), maximum(zs)),
        shading = NoShading)
end

# `dem` is an X-Y elevation Raster matching `snow_raster`'s X-Y grid; each
# `panel_indices[k]` is a Ti index whose snow slice is draped on top.
function plot_drape_panels(snow_raster, dem, panel_indices, panel_labels,
                            variable_label, filename;
                            cmap = Reverse(:ice), unit_label = "cm",
                            ncols = 3)
    clims = drape_clims(snow_raster)
    clims === nothing && return
    xs = collect(lookup(dem, X))
    ys = collect(lookup(dem, Y))
    zs = Float64.(parent(dem))
    zs_snow = zs .+ DRAPE_Z_OFFSET
    nrows = cld(length(panel_indices), ncols)
    figure = Figure(size = (ncols * 460, nrows * 380 + 80))
    Label(figure[0, 1:ncols], variable_label;
        fontsize = 15, halign = :left, padding = (10, 0, 0, 0))
    for k in eachindex(panel_indices)
        row, col = fldmod1(k, ncols)
        ax = Axis3(figure[row, col];
            title = panel_labels[k],
            xlabel = "Lon", ylabel = "Lat", zlabel = "m",
            azimuth = -π/4, elevation = π/8,
            aspect = (1, 1, 0.35))
        add_terrain_base!(ax, xs, ys, zs)
        surface!(ax, xs, ys, zs_snow;
            color = drape_color(view(snow_raster; Ti = panel_indices[k])),
            colormap = cmap, colorrange = clims,
            nan_color = :transparent, shading = NoShading)
    end
    Colorbar(figure[1:nrows, ncols + 1];
        colormap = cmap, colorrange = clims, label = unit_label)
    save(filename, figure)
    println("  Saved $filename")
end

function animate_drape(snow_raster, dem, frame_labels,
                       variable_label, filename;
                       cmap = Reverse(:ice), unit_label = "cm", framerate = 4)
    clims = drape_clims(snow_raster)
    clims === nothing && return
    xs = collect(lookup(dem, X))
    ys = collect(lookup(dem, Y))
    zs = Float64.(parent(dem))
    zs_snow = zs .+ DRAPE_Z_OFFSET
    n_frames = size(snow_raster, Ti)
    figure = Figure(size = (820, 720))
    title_obs = Observable("$variable_label · $(frame_labels[1])")
    ax = Axis3(figure[1, 1];
        title = title_obs,
        xlabel = "Lon", ylabel = "Lat", zlabel = "m",
        azimuth = -π/4, elevation = π/8,
        aspect = (1, 1, 0.35))
    add_terrain_base!(ax, xs, ys, zs)
    color_obs = Observable(drape_color(view(snow_raster; Ti = 1)))
    surface!(ax, xs, ys, zs_snow;
        color = color_obs, colormap = cmap, colorrange = clims,
        nan_color = :transparent, shading = NoShading)
    Colorbar(figure[1, 2]; colormap = cmap, colorrange = clims, label = unit_label)
    record(figure, filename, 1:n_frames; framerate) do k
        color_obs[] = drape_color(view(snow_raster; Ti = k))
        title_obs[] = "$variable_label · $(frame_labels[k])"
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
noon_dates      = collect(ti_vector[noon_indices])

point_colors = [
    RGBf(0.30, 0.30, 0.30),   # valley — graphite
    RGBf(0.55, 0.55, 0.55),   # mid    — mid grey
    RGBf(0.63, 0.0,  0.0),    # summit — brick red accent
]

snow_figure = Figure(size = (1000, 500))
snow_axis = Axis(snow_figure[1, 1];
    title  = "Snow depth, $vector_year",
    xlabel = "Date", ylabel = "Snow depth (cm)",
)
for (i, name) in enumerate(point_labels)
    depths_cm = view(vector_output.snow_depth; point = i, Ti = noon_indices)
    lines!(snow_axis, noon_dates, parent(depths_cm);
        label = name, color = point_colors[i], linewidth = 1.4)
end
axislegend(snow_axis; position = :rt, framevisible = false)
save(joinpath(FIG_DIR, "point_snow_depth.png"), snow_figure)
println("  Saved $(joinpath(FIG_DIR, "point_snow_depth.png"))")

target_date  = Date(vector_year, 3, 1)
day_indices  = findall(t -> Date(t) == target_date, ti_vector)
hours_of_day = hour.(ti_vector[day_indices])

extract_day(layer; kw...) = parent(view(layer; Ti = day_indices, kw...))

variable_panels = [
    ("Air, 2 m",     i -> extract_day(vector_output.air_temperature;  point = i, height = 2)),
    ("Air, 1 cm",    i -> extract_day(vector_output.air_temperature;  point = i, height = 1)),
    ("Soil surface", i -> extract_day(vector_output.soil_temperature; point = i, depth  = 1)),
    ("Soil, 10 cm",  i -> extract_day(vector_output.soil_temperature; point = i, depth  = 4)),
]
temperatures_figure = Figure(size = (1200, 800))
Label(temperatures_figure[0, 1:2], "$target_date";
    fontsize = 15, halign = :left, padding = (10, 0, 0, 0))
for (k, (label, extract)) in enumerate(variable_panels)
    row, col = fldmod1(k, 2)
    ax = Axis(temperatures_figure[row, col];
        title = label, xlabel = "Hour of day", ylabel = "°C")
    for (i, name) in enumerate(point_labels)
        lines!(ax, hours_of_day, extract(i);
            label = name, color = point_colors[i], linewidth = 1.4)
    end
    k == 1 && axislegend(ax; position = :lt, framevisible = false)
end
save(joinpath(FIG_DIR, "point_temperatures.png"), temperatures_figure)
println("  Saved $(joinpath(FIG_DIR, "point_temperatures.png"))")

# ---------------------------------------------------------------------------
# Raster demo plots (small high-res Mont Aigoual box, July 2000)
# ---------------------------------------------------------------------------

println("Raster demo plots...")
# Only the soil-temperature-5cm drape panel is used by the talk; load just
# that layer and slice to July to keep things compact.
raster_output = RasterStack(joinpath(OUTPUT_DIR, "raster_summer.nc"))

ti_raster = lookup(raster_output, Ti)
july = findall(d -> month(d) == 7, ti_raster)

soil_temperature_5cm = view(raster_output.soil_temperature; Ti = july, depth = 3)

panel_indices = [1, 7, 10, 13, 16, 19]
panel_labels  = ["Midnight", "Dawn", "Mid-morning", "Midday", "Mid-afternoon", "Dusk"]

# Load SRTM at the raster demo's extent and resample onto the run grid so
# DEM X/Y match the soil-temperature raster, then drape on the 3-D mesh.
raster_area = Extent(X = (3.554, 3.608), Y = (44.095, 44.149))
raster_dem  = Rasters.resample(
    read(crop(Raster(SRTM; extent = raster_area, lazy = true, missingval = Int16(0));
              to = raster_area, touches = true));
    to = first(raster_output),
)

plot_drape_panels(soil_temperature_5cm, raster_dem, panel_indices, panel_labels,
    "Soil temperature at 5 cm (°C)",
    joinpath(FIG_DIR, "grid_soil_temperature_5cm_drape.png");
    cmap = Reverse(:RdYlBu), unit_label = "°C")

# ---------------------------------------------------------------------------
# Snow demo plot — yearly cold-season GIF
# ---------------------------------------------------------------------------
# Snow depth is draped on a 3-D mesh built from the aggregated SRTM DEM.
# Zero-depth cells are mapped to NaN via `replace_missing(rast, 0)` + the
# `drape_color` helper, so bare ground shows the underlying terrain.

snow_area = Extent(X = (3.33, 3.83), Y = (43.87, 44.37))
snow_dem  = aggregate(mean,
    read(crop(Raster(SRTM; extent = snow_area, lazy = true, missingval = Int16(0));
              to = snow_area, touches = true)),
    6)

println("Year-long snow plot...")
yearly_snow_output = RasterStack(joinpath(OUTPUT_DIR, "snow_yearly.nc"))
ti_yearly          = lookup(yearly_snow_output, Ti)
yearly_year        = Dates.year(first(ti_yearly))
yearly_snow_depth  = yearly_snow_output.snow_depth
noon_indices_y     = findall(t -> hour(t) == 12, ti_yearly)

label_for(idx) = Dates.format(ti_yearly[idx], "u dd")

cold_idx   = cold_season_window(yearly_snow_depth, noon_indices_y; buffer_days = 7)
cold_idx === nothing && error("No snow found in $(yearly_year) yearly run")
cold_snow  = replace_missing(view(yearly_snow_depth; Ti = cold_idx), 0)
animate_drape(cold_snow, snow_dem, label_for.(cold_idx),
    "Snow depth (cm) — daily cold-season $(yearly_year)",
    joinpath(FIG_DIR, "grid_snow_depth_yearly.gif"); framerate = 4)
