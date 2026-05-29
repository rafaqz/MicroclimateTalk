# Point microclimate at three Cévennes sites, winter 2010.
#
# Daily NCEP reanalysis forcing; no snow model. The three sites share the
# same NCEP cell (~2.5° resolution) — terrain lapse-correction and aspect
# alone explain the differences. The solve is only ~14 s so we don't cache
# it; `plot_demos.jl` re-includes this file when plotting.

include("shared.jl")

# Valley floor → south-facing slope → summit. Same NCEP cell, very
# different microclimates from elevation + aspect.
points = [
    (3.520, 44.143),   # Camprieu (valley, ~900 m)
    (3.582, 44.115),   # South-facing slope (~1200 m)
    (3.581, 44.122),   # Mont Aigoual summit (~1567 m)
]

year = 2010

model = MicroMapModel(;
    micro_model = MicroModel(;
        depths,
        heights,
        soil_profile          = example_soil_profile(depths),
        soil_properties_model = example_soil_properties_model(),
        soil_hydraulic_model  = example_soil_hydraulic_model(depths),
        snow_model            = NoSnow(),
    ),
    dem_source              = SRTM,
    weather_source          = NCEP{SurfaceGauss},
    surface_albedo_source   = 0.15,
    roughness_height_source = 0.004u"m",
)

problem = MicroVectorProblem(;
    model,
    points,
    years = year:year,
    init  = (; soil_moisture = fill(0.25, length(depths))),
)

output = solve(problem)
