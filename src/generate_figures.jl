# Driver: solve every demo (or reuse cached NetCDF in `output/`), then
# render every figure used in the slides. Plotting is in a separate file
# so it can be iterated on without re-solving.

println("=== Vector demo ===")
include("vector_demo.jl")

println()
println("=== Raster demo (small box, summer) ===")
include("raster_demo.jl")

println()
println("=== Snow demos (Cévennes window) ===")
include("snow_demo.jl")

println()
println("=== Plots ===")
include("plot_demos.jl")
