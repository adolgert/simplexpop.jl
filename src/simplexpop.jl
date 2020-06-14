module simplexpop
# These have no dependencies.
include("sort.jl")
include("hilbert.jl")

# Manipulate rasters.
include("assignweights.jl")

# Finite elements.
include("decplex.jl")
include("fibgrid.jl")
include("hexlattice.jl")
end # module
