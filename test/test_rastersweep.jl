using ArchGDAL
using Test

# Here's our stub setup.
# Coarse grid from ul=(-10, 16) to lr=(10, 0) in xy
# Give it w=20, h=16.
# Fine grid from ul=(-5, 12) to lr=(5, 4) in xy
# Give it w=10*5=50, h=8*5=40, so 5x the resolution of the coarse.
# CCCCCCCCCC|CCCCCCCCCC
# CCCCCCCCCC|CCCCCCCCCC
# CCCCCCCCCC|CCCCCCCCCC
# CCCCCCCCCC|CCCCCCCCCC
# CCCCCFFFFF|FFFFFCCCCC
# CCCCCFFFFF|FFFFFCCCCC
# CCCCCFFFFF|FFFFFCCCCC
# CCCCCFFFFF|FFFFFCCCCC
# CCCCCFFFFF|FFFFFCCCCC
# ----------|----------
# CCCCCFFFFF|FFFFFCCCCC
# CCCCCFFFFF|FFFFFCCCCC
# CCCCCFFFFF|FFFFFCCCCC
# CCCCCFFFFF|FFFFFCCCCC
# CCCCCFFFFF|FFFFFCCCCC
# CCCCCCCCCC|CCCCCCCCCC
# CCCCCCCCCC|CCCCCCCCCC
# CCCCCCCCCC|CCCCCCCCCC
# CCCCCCCCCC|CCCCCCCCCC
coarse_geo = Geo(
      Float64[-10, 16] , # origin
      Float64[1 0.0; 0.0 -1],  # pixel size
   )
fine_geo = Geo(
      Float64[-5, 12] , # origin
      Float64[.2 0.0; 0.0 -.2],  # pixel size
   )


rr = pixel_rect(coarse_geo, [1, 1])
@test rr ≈ [-10 -9; 16 15]

ci1 = cover_rect([-10 -9; 16 15], coarse_geo, (26, 10))
@test ci1[1, 1] == CartesianIndex(1, 1)

ci2 = cover_rect([-10 -8; 16 15], coarse_geo, (26, 10))
cci2 = collect(ci2)
@test length(cci2) == 2


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
   assignweights(hrsl_path, landscan_path, out_path)
   @test ispath(out_path)
   @test stat(out_path).size > 1210320848
else
   println("Could not find HRSL or landscan")
   println("$(hrsl_path)=$(hrsl_exists)")
   println("$(landscan_path)=$(landscan_exists)")
end
