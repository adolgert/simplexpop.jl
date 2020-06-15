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

# From sort
export countedinsertionsort!

# From hilbert
export encode_hilbert_zero, decode_hilbert_zero
export encode_hilbert, decode_hilbert, hilbert_order

# From assignweights.jl
export PixelGrid, pixelgrid, DataGrid
export assignweights, intersection_area
export geo_to_transform, geo_to_xy_corner, geo_to_xy_center
export xy_bounds, pixel_containing, ij_cover_rect, block_cover_ij_rect
export pixels_of_block, corner_pixel, pixelgrid, load_single_block
export load_pixel_grid

end # module
