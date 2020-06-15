using ArchGDAL

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


xyb1 = xy_bounds(origin, transform, coarse_width_height)
@test xyb1[1][1] ≈ -180
@test xyb1[1][2] ≈ 90
@test xyb1[2][1] ≈ 180
@test xyb1[2][2] ≈ -90

xy5 = geo_to_xy_center(origin, transform, 24, 36)
ij1 = pixel_containing(origin, transform, xy5)
@test ij1[1] == 24
@test ij1[2] == 36

xyb2 = xy_bounds(origin, transform, coarse_width_height)
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

pb1 = pixels_of_block(CartesianIndex(1, 1), (5, 3))
@test length(pb1) == 15
@test size(pb1) == (5, 3)
@test pb1[1][1] == 1
@test pb1[1][2] == 1
@test pb1[15][1] == 5
@test pb1[15][2] == 3
cp1 = corner_pixel(CartesianIndex(1, 1), (5, 3))
@test cp1 == pb1[1]
pb2 = pixels_of_block(CartesianIndex(3, 2), (5, 3))
@test length(pb2) == 15
@test size(pb2) == (5, 3)
@test pb2[1][1] == 11
@test pb2[1][2] == 4
@test pb2[15][1] == 15
@test pb2[15][2] == 6
cp2 = corner_pixel(CartesianIndex(3, 2), (5, 3))
@test cp2 == pb2[1]

hrsl_path = "~/data/inputs/HRSL/uganda2018/hrsl_uga_pop.tif"
landscan_path = "~/data/inputs/landscan/LandScan Global 2018/lspop2018/w001001.adf"
hrsl_exists = isfile(expanduser(hrsl_path))
landscan_exists = isfile(expanduser(landscan_path))

if hrsl_exists && landscan_exists
   ArchGDAL.read(expanduser(landscan_path)) do ls_ds
      band = ArchGDAL.getband(ls_ds, 1)
      pg1 = pixelgrid(band)

      sb1 = load_single_block(band, pg1, (1, 1))
      sb2 = load_single_block(band, pg1, (20, 20))

      # Load the whole thing.
      dg1 = load_pixel_grid(band, pg1)
   end

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
