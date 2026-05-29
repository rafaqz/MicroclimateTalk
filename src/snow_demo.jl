# Snow demos — two scenarios over a ~55 km × 55 km Cévennes window
# centred on Mont Aigoual. SRTM is loaded at native resolution and
# aggregated by 6 to a 100×100 template so the runs finish in minutes
# while still showing real spatial structure.
#
# Outputs cached to `output/snow_spring.nc` + `output/snow_yearly.nc`.

include("shared.jl")

snow_area = Extent(X = (3.33, 3.83), Y = (43.87, 44.37))
snow_template_native = read(crop(
    Raster(SRTM; extent = snow_area, lazy = true, missingval = Int16(0));
    to = snow_area, touches = true,
))
snow_template = aggregate(mean, snow_template_native, 6)
println("Snow template grid: $(size(snow_template))")

snow_micro_model = MicroModel(;
    depths,
    heights,
    soil_profile          = example_soil_profile(depths),
    soil_properties_model = example_soil_properties_model(),
    soil_hydraulic_model  = example_soil_hydraulic_model(depths),
    snow_model            = SnowModel(),
)

# ---------------------------------------------------------------------------
# Scenario 1 — spring snowmelt (TerraClimate, 30 cm initial snowpack)
# ---------------------------------------------------------------------------
# Monthly resolution triggers NonConsecutiveDayMode: each month is a
# representative day re-using the 30 cm initial snowpack.

spring_map_model = MicroMapModel(;
    micro_model             = snow_micro_model,
    dem_source              = SRTM,
    weather_source          = TerraClimate{Historical},
    surface_albedo_source   = 0.15,
    roughness_height_source = 0.004u"m",
)

spring_problem = MicroRasterProblem(;
    model    = spring_map_model,
    area     = snow_area,
    years    = 2000:2000,
    template = snow_template,
    init     = (; snow_depth = 30.0u"cm"),
)

spring_output = solve(spring_problem)

write(joinpath(OUTPUT_DIR, "snow_spring.nc"), strip_to_canonical(spring_output); force = true)

# ---------------------------------------------------------------------------
# Scenario 2 — year-long snow (NCEP daily, ConsecutiveDayMode)
# ---------------------------------------------------------------------------
# Daily resolution triggers ConsecutiveDayMode automatically. Snowpack
# accumulates and melts day to day across the year.

yearly_map_model = MicroMapModel(;
    micro_model             = snow_micro_model,
    dem_source              = SRTM,
    weather_source          = NCEP{SurfaceGauss},
    surface_albedo_source   = 0.15,
    roughness_height_source = 0.004u"m",
)

yearly_problem = MicroRasterProblem(;
    model    = yearly_map_model,
    area     = snow_area,
    years    = 2010:2010,
    template = snow_template,
    init     = (;
        soil_moisture = fill(0.25, length(depths)),
        snow_depth    = 0.0u"cm",
    ),
)

yearly_output = solve(yearly_problem)

write(joinpath(OUTPUT_DIR, "snow_yearly.nc"), strip_to_canonical(yearly_output); force = true)
