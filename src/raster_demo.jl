# Gridded microclimate over Mont Aigoual (Cévennes, southern France).
# Small high-resolution box at native SRTM resolution; TerraClimate forcing.
# Output cached to `output/raster_summer.nc`.

include("shared.jl")

# Mont Aigoual: 1567 m peak, dramatic N/S aspect contrast,
# 30 km NW of Montpellier. ~0.054° square centred on the summit.
area = Extent(X = (3.554, 3.608), Y = (44.095, 44.149))
year = 2000

micro_model = MicroModel(;
    depths,
    heights,
    soil_profile          = example_soil_profile(depths),
    soil_properties_model = example_soil_properties_model(),
    soil_hydraulic_model  = example_soil_hydraulic_model(depths),
    snow_model            = NoSnow(),
)

map_model = MicroMapModel(;
    micro_model,
    dem_source              = SRTM,
    weather_source          = TerraClimate{Historical},
    surface_albedo_source   = 0.15,
    roughness_height_source = 0.004u"m",
)

problem = MicroRasterProblem(;
    model    = map_model,
    area,
    years    = year:year,
    template = SRTM,
)

output = solve(problem)

write(joinpath(OUTPUT_DIR, "raster_summer.nc"), strip_to_canonical(output); force = true)
