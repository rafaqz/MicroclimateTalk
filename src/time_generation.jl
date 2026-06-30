# Run every demo + plot script and record how long each took.
# Run with:
#     julia --project=. src/time_generation.jl 2>&1 | tee output/timings.txt
# Or just let it stream to stdout — the summary at the end is the same either way.

using Dates

const TIMING_DIR = joinpath(@__DIR__, "..", "output")
mkpath(TIMING_DIR)

parts = [
    "Vector demo" => "vector_demo.jl",
    "Raster demo" => "raster_demo.jl",
    "Snow demo"   => "snow_demo.jl",
    "Plot demos"  => "plot_demos.jl",   # re-includes vector_demo internally (~14 s)
]

summary = Pair{String,Float64}[]

println("\nMicroclimate talk data-generation timings")
println("Started: ", Dates.now(), "  Julia ", VERSION, "  threads ", Threads.nthreads())

for (name, file) in parts
    println("\n", "="^60)
    println("  ", name, " — ", file)
    println("="^60)
    flush(stdout)
    t = @elapsed include(file)
    push!(summary, name => t)
    println("==> ", name, ": ", round(t, digits=2), " s")
end

println("\n", "="^60)
println("  Summary")
println("="^60)
total = 0.0
for (name, t) in summary
    println(rpad(name * ":", 24), lpad(round(t, digits=2), 8), " s")
    global total += t
end
println("-"^36)
println(rpad("Total:", 24), lpad(round(total, digits=2), 8), " s")
println("Finished: ", Dates.now())
