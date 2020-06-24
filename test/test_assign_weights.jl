using ArchGDAL
using Test
include("../src/hilbert.jl")

coarse_geo = [-180.0, 0.0083333333333333, 0.0, 89.99999999999929, 0.0, -0.0083333333333333]
coarse_width_height = (43200, 21600)
fine_geo = [29.5734769457, 0.00027777777999999997, 0.0, 4.22813426181, 0.0, -0.00027777778]
fine_width_height = (19536, 20540)

origin, transform = geo_to_transform(coarse_geo)

xy1 = geo_to_xy_corner(origin, transform, 1, 1)
@test xy1[1] ≈ -180
@test xy1[2] ≈ 90
xy2 = geo_to_xy_corner(origin, transform, coarse_width_height[1] + 1, coarse_width_height[2] + 1)
@test xy2[1] ≈ 180
@test xy2[2] ≈ -90

xy3 = geo_to_xy_corner(origin, transform, 24, 36)
xy4 = geo_to_xy_center(origin, transform, 24, 36)
@test xy4[1] ≈ xy3[1] + coarse_geo[2] / 2
@test xy4[2] ≈ xy3[2] + coarse_geo[6] / 2 # Note that the offset is negative, -0.008333.

coarse_pg = PixelGrid((1, 1), coarse_width_height, geo_to_transform(coarse_geo)...)
xyb1 = xy_bounds(coarse_pg)
@test xyb1[1][1] ≈ -180
@test xyb1[1][2] ≈ 90
@test xyb1[2][1] ≈ 180
@test xyb1[2][2] ≈ -90

xy5 = geo_to_xy_center(origin, transform, 24, 36)
ij1 = pixel_containing(origin, transform, xy5)
@test ij1[1] == 24
@test ij1[2] == 36

xyb2 = xy_bounds(coarse_pg)
# Nudge the corner just a hair so it isn't right on the boundary.
xyb3 = ([xyb2[1][1], xyb2[1][2]], [xyb2[2][1] - 1e-5, xyb2[2][2] + 1e-5])
ij_corners = ij_cover_rect(origin, transform, xyb3)
@test ij_corners[1][1] == 1
@test ij_corners[1][2] == 1
@test ij_corners[2][1] == coarse_width_height[1]
@test ij_corners[2][1] == coarse_width_height[1]

bci_trial = [
   (([1, 1], [10, 10]), (200, 200), ((1, 1), (1, 1))),
   (((1, 1), (200, 200)), (200, 200), ((1, 1), (1, 1))),
   (((1, 1), (201, 250)), (200, 300), ((1, 1), (2, 1))),
   (((1, 1), (200, 301)), (200, 300), ((1, 1), (1, 2))),
   (((204, 1), (250, 601)), (200, 300), ((2, 1), (2, 3))),
   (((5, 50), (5, 50)), (2, 10), ((3, 5), (3, 5))),
   (((4, 50), (5, 50)), (2, 10), ((2, 5), (3, 5))),
   (((5, 40), (5, 50)), (2, 10), ((3, 4), (3, 5))),
   (((5, 50), (6, 50)), (2, 10), ((3, 5), (3, 5))),
   (((5, 50), (7, 50)), (2, 10), ((3, 5), (4, 5))),
]
for (ij_rect, bs, result) in bci_trial[8:8]
   ij = (Int.(ij_rect[1]), Int.(ij_rect[2]))
   bl1 = block_cover_ij_rect(ij, Int.(bs))
   @test bl1 == result
end

rectish = [
   ((1, 1), ((0, 2), (2, 0)), true),
   ((3, 1), ((0, 4), (2, 0)), false),
   ((1, 3), ((0, 4), (2, 0)), true),
   ((-1, 1), ((0, 4), (2, 0)), false),
   ((1, -1), ((0, 4), (2, 0)), false),
   ((-1, -1), ((0, 4), (2, 0)), false),
   ((3, 5), ((0, 4), (2, 0)), false),
   ((1.5, 3.5), ((1, 4), (2, 3)), true),
   ((.5, 3.5), ((1, 4), (2, 3)), false),
   ((1.5, 2.5), ((1, 4), (2, 3)), false),
]
for (pt, rect, result) in rectish
   @test pointinrect(pt, rect) == result
end

boundsrect = ([-6.0, 12.0], [-5.0, 11.0])
pt = [-4.9, 11.9]
@test !pointinrect(pt, boundsrect)

# Here's our stub setup.
# Coarse grid from ul=(-10, 16) to lr=(10, 0) in xy
# Give it w=20, h=16.
# Fine grid from ul=(-5, 12) to lr=(5, 4) in xy
# Give it w=10*5=50, h=8*5=40, so 5x the resolution of the coarse.
function make_coarse(A)
   StubBand{Int32}(
      A,
      [6, 2],  # blocksize
      Float64[-10, 16] , # origin
      Float64[1 0.0; 0.0 -1],  # pixel size
   )
end
coarse_band = make_coarse(fill(Int32(200), 20, 16))
function make_fine(A)
   StubBand{Float64}(
      A,
      [5, 3],  # blocksize
      Float64[-5, 12] , # origin
      Float64[.2 0.0; 0.0 -.2],  # pixel size
   )
end
fine_band = make_fine(fill(2.5, 50, 40))

coarse_pg = pixelgrid(coarse_band)
fine_pg = pixelgrid(fine_band)
# This should be the center of the pixel.
pc1 = pixel_containing(coarse_pg.origin, coarse_pg.transform, [-4.5, 11.5])
@test pc1 == [6, 5]
pc2 = pixel_containing(coarse_pg.origin, coarse_pg.transform, [-4.001, 11.001])
@test pc2 == [6, 5]
r1 = ij_cover_rect(coarse_pg.origin, coarse_pg.transform, ([-5, 12], [-4, 11]))
@test r1[1] == [6, 5]
@test r1[2][1] <= 7  # The cover rect can bring in a neighboring pixel.
@test r1[2][2] <= 6
r2 = ij_cover_rect(coarse_pg.origin, coarse_pg.transform, ([-4.999, 11.99], [4.999, 4.01]))
@test r2 == ([6, 5], [15, 12])
large_crop = crop_to(coarse_pg, xy_bounds(fine_pg))
@test large_crop.ul == (6, 5)


sb1 = load_single_block(coarse_band, coarse_pg, (1, 1))
sb2 = load_single_block(coarse_band, coarse_pg, (3, 3))
@test sb2.A[13, 5] == 200
# This is past the last block in x. There are 4 blocks.
@test_throws AssertionError sb3 = load_single_block(coarse_band, coarse_pg, (5, 3))

# Create an array with unique values we can check.
AA = Int32[encode_hilbert(x[1], x[2]) for x in CartesianIndices((20,16))]
copya = copy(AA)
cb2 = make_coarse(AA)
dg1 = load_data_grid(cb2, coarse_pg)
@test size(dg1.A) == (6 * 4, 16)  # Must be integral blocks to cover the domain.
@test dg1.A == copya
tune_band(fine_band, coarse_band, typemax(Int64))
@test fine_band.A[10, 10] ≈ 200 / 25

# Now try some nonzero structure.
A2 = fill(-floatmax(Float64), 50, 40)
fb2 = make_fine(A2)
tune_band(fb2, coarse_band, typemax(Int64))
@test all(fb2.A .== -floatmax(Float64))

# Now try some nonzero structure.
A3 = fill(-floatmax(Float64), 50, 40)
A3[1:5, 1:5] .= 1  # first 1x1 all filled.
A3[3:3, 7:7] .= 1  # Only one nonzero in this 1x1.
A3[6:7, 6:7] .= 1  # Four nonzero here.
fb3 = make_fine(A3)
@test sum(fb3.A .> 0) == 25 + 1 + 4
tune_band(fb3, coarse_band, 2)

@test all(fb3.A[1:5, 1:5] .≈ 200 / 25)
@test all(fb3.A[3:3, 7:7] .≈ 100)  # This will hit the max.
@test all(fb3.A[6:7, 6:7].≈ 200 / 4)



## A live data test.
hrsl_path = "~/data/inputs/HRSL/uganda2018/hrsl_uga_pop.tif"
landscan_path = "~/data/inputs/landscan/LandScan Global 2018/lspop2018/w001001.adf"
hrsl_exists = isfile(expanduser(hrsl_path))
landscan_exists = isfile(expanduser(landscan_path))

if hrsl_exists && landscan_exists
   out_path = "hrsl_ls_adjusted.tif"
   if ispath(out_path)
      rm(out_path)
   end
   assignweights(hrsl_path, landscan_path, out_path; max_pixels = 1)
   @test ispath(out_path)
   @test stat(out_path).size > 1210320848
else
   println("Could not find HRSL or landscan")
   println("$(hrsl_path)=$(hrsl_exists)")
   println("$(landscan_path)=$(landscan_exists)")
end
