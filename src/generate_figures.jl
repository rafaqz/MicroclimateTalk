# Driver: run both demos to regenerate every figure used in the slides.
# Used locally and by the GitHub Actions build.

println("=== Point demo (Montpellier) ===")
include("point_demo.jl")

println()
println("=== Grid demo (Mont Aigoual) ===")
include("grid_demo.jl")
