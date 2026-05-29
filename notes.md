# Need to discuss
- performance/accuracy tradeoffs when making large/high res rasters
  - need to start with a fast language for this to be possible
  - ability to mix and match formulations
- units make writing accurate code easy
- ease of collaboration and sharing vs C++

- hard-coded datasets as one approach to performance mitigation

# Audience
- people doing high-res SDMs
- people are measuring microclimates (soil temp group)
  - how do we feed measurements into model improvements
  - fitting parameters from data, inverse fits

# Names ??
- BiophysicalGrids vs BiophysicalRasters vs MicroclimateRasters vs MicroclimateGrids
  - what level of abstraction is it at



- need high performance
- lots of models - lots of busywork to switch between models and datasets

We need 
- declarative models and datasets


Julia
- multliple dispatch and a type system
- saves busywork
- compiler handles your problems
- object describes the entire run
- DifferentialEquations/SciML

- R + Fortran / R + C

THe best bling:

- animation of snow melting at the end


Tutorial
- come with vscode/julia installed use juilaup
- show basic unit
- use SolarRadiation
- use FluidProperties ?


# Talk structure
Introduce the problems
Why julia
Demonstrate packages point -> spatial -> animation
Rest of the ecosystem/context
Question

