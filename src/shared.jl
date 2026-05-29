# Shared imports + run paths used by every demo and the plot script.

using Microclimate
using Microclimate: example_soil_profile,
    example_soil_properties_model, example_soil_hydraulic_model
using MicroclimateMapper
using RasterDataSources
using Rasters
using Rasters: Ti, X, Y, lookup, aggregate
using Rasters.Extents: Extent
using Statistics: mean
using Dates
using Printf
using Unitful
using CairoMakie

const FIG_DIR    = joinpath(@__DIR__, "..", "figures")
const OUTPUT_DIR = joinpath(@__DIR__, "..", "output")
mkpath(FIG_DIR)
mkpath(OUTPUT_DIR)

# Soil + air discretisation shared across every demo.
depths  = [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 200.0]u"cm"
heights = [0.01, 2.0]u"m"
