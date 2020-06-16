using ArchGDAL
using Test

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

bl1 = block_cover_ij_rect(([1, 1], [10, 10]), (200, 200))
@test bl1[1][1] == 1
@test bl1[1][2] == 1
@test bl1[2][1] == 1
@test bl1[2][2] == 1
bl2 = block_cover_ij_rect(((1, 1), (200, 200)), (200, 200))
@test bl2[1][1] == 1
@test bl2[1][2] == 1
@test bl2[2][1] == 1
@test bl2[2][2] == 1
bl3 = block_cover_ij_rect(((1, 1), (201, 250)), (200, 300))
@test bl3[1][1] == 1
@test bl3[1][2] == 1
@test bl3[2][1] == 2
@test bl3[2][2] == 1
bl4 = block_cover_ij_rect(((1, 1), (200, 301)), (200, 300))
@test bl4[1][1] == 1
@test bl4[1][2] == 1
@test bl4[2][1] == 1
@test bl4[2][2] == 2
bl5 = block_cover_ij_rect(((204, 1), (250, 601)), (200, 300))
@test bl5[1][1] == 2
@test bl5[1][2] == 1
@test bl5[2][1] == 2
@test bl5[2][2] == 3

rectish = [
   ((1, 1), ((0, 0), (2, 2)), true),
   ((3, 1), ((0, 0), (2, 4)), false),
   ((1, 3), ((0, 0), (2, 4)), true),
   ((-1, 1), ((0, 0), (2, 4)), false),
   ((1, -1), ((0, 0), (2, 4)), false),
   ((-1, -1), ((0, 0), (2, 4)), false),
   ((3, 5), ((0, 0), (2, 4)), false),
   ((1.5, 3.5), ((1, 3), (2, 4)), true),
   ((.5, 3.5), ((1, 3), (2, 4)), false),
   ((1.5, 2.5), ((1, 3), (2, 4)), false),
]
for (pt, rect, result) in rectish
   @test pointinrect(pt, rect) == result
end

# Here's our stub setup.
# Coarse grid from ul=(-10, 16) to lr=(10, 0) in xy
# Give it w=20, h=16.
# Fine grid from ul=(-5, 12) to lr=(5, 4) in xy
# Give it w=10*5=50, h=8*5=40, so 5x the resolution of the coarse.
coarse_band = StubBand{Int32}(
   fill(400, 20, 16),
   [6, 2],  # blocksize
   Float64[-10, 16] , # origin
   Float64[1 0.0; 0.0 -1],  # pixel size
)
fine_band = StubBand{Float64}(
   fill(2.5, 50, 40),
   [5, 3],  # blocksize
   Float64[-5, 12] , # origin
   Float64[.2 0.0; 0.0 -.2],  # pixel size
)

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
large_crop = crop_to(coarse_pg, fine_pg)

sb1 = load_single_block(coarse_band, pg1, (1, 1))
sb2 = load_single_block(coarse_band, pg1, (3, 3))
@test sb2.A[13, 5] == 400
# This is past the last block in x. There are 4 blocks.
@test_throws AssertionError sb3 = load_single_block(coarse_band, pg1, (5, 3))

# Load the whole thing.
dg1 = load_data_grid(coarse_band, pg1)
@test size(dg1.A) == (6 * 4, 16)  # Must be integral blocks to cover the domain.

tune_band(fine_band, coarse_band, typemax(Int64))

## A live data test.
hrsl_path = "~/data/inputs/HRSL/uganda2018/hrsl_uga_pop.tif"
landscan_path = "~/data/inputs/landscan/LandScan Global 2018/lspop2018/w001001.adf"
hrsl_exists = isfile(expanduser(hrsl_path))
landscan_exists = isfile(expanduser(landscan_path))

if hrsl_exists && landscan_exists && false
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
