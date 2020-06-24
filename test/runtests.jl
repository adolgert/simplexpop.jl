using simplexpop
using Test

# Make these the same order as their includes in simplexpop.jl.
@testset "sort.jl" begin include("test_sort.jl") end
@testset "hilbert.jl" begin include("test_hilbert.jl") end
@testset "assign_weights.jl" begin include("test_assign_weights.jl") end
