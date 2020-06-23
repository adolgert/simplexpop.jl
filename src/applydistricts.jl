using ArchGDAL
using DelimitedFiles


function applydistricts(districts_path, hrsl_path)
end

function layerinfo(layer)
   @show ArchGDAL.getname(layer)
   @show ArchGDAL.nfeature(layer)
   featuredefn = ArchGDAL.layerdefn(layer)
   nfield = ArchGDAL.nfield(featuredefn)
   ngeom = ArchGDAL.nfield(featuredefn)
   @show ngeom
   for fieldidx in 1:nfield
      field = ArchGDAL.getfielddefn(featuredefn, fieldidx - 1)
      name = ArchGDAL.getname(field)
      type = ArchGDAL.gettype(field)
      @show (name, type, fieldidx)
   end
   for geomidx in 1:1
      geom = ArchGDAL.getgeomdefn(featuredefn, geomidx - 1)
      name = ArchGDAL.getname(geom)
      type = ArchGDAL.gettype(geom)
      @show (name, type, geomidx)
   end
end


function overgeometry(toapply, layer)
   nfeature = ArchGDAL.nfeature(layer)
   result = fill(("", zero(Base.return_types(toapply)[1])), nfeature)
   for featureidx in 1:nfeature
      ArchGDAL.getfeature(layer, featureidx - 1) do feature
         admin = ArchGDAL.getfield(feature, 7 - 1)
         geom = ArchGDAL.getgeom(feature, 0)
         envelope = ArchGDAL.envelope(geom)
         @show admin
         result[featureidx] = (admin, toapply(geom))
      end
   end
   result
end


# You need to include rastersweep.jl
struct RasterRead
   data
   no_data
   geo::Geo
end


function build_countup(raster)
   function count_up(geometry)::Float64
      envelope = ArchGDAL.envelope(geometry)
      # upper-left on the left column. lower-right in the right column.
      rect = [envelope.MinX envelope.MaxX; envelope.MaxY envelope.MinY]
      indices = cover_rect(rect, raster.geo, size(raster.data))
      total = 0.0
      # We use a GDAL point to test for inside polygon.
      geom_point = ArchGDAL.createpoint()
      ArchGDAL.addpoint!(geom_point, 0.0, 0.0)
      within = fill(CartesianIndex(0, 0), length(indices))
      inside_cnt = 0
      for idx in indices
         if raster.data[idx] != raster.no_data
            xy = pixel_center(raster.geo, [idx[1], idx[2]])
            ArchGDAL.setpoint!(geom_point, 0, xy[1], xy[2])
            if ArchGDAL.within(geom_point, geometry)
               value = raster.data[idx]
               @assert value >= 0
               total += value
               inside_cnt += 1
               within[inside_cnt] = idx
            end
         end
      end
      largest = 0.0
      for inside_idx in 1:inside_cnt
         idx = within[inside_idx]
         # All values are 1000 times the actual fraction for numerical accuracy.
         pop_rate = 1000 * raster.data[idx] / total
         @assert pop_rate >= 0
         raster.data[idx] = pop_rate
         largest = max(largest, pop_rate)
      end
      largest
   end
end


"""
Make a dataset with the same size and transform as the given dataset
but it has two bands, one with floats and one with integers.
"""
function create_fraction(dataset, path)
   b = ArchGDAL.getband(dataset, 1)
   w = ArchGDAL.width(b)
   h = ArchGDAL.height(b)
   driver = ArchGDAL.getdriver("GTiff")
   ds = ArchGDAL.create(path; driver = driver, width = w, height = h, nbands = 1, dtype = Int32)
   ArchGDAL.setgeotransform!(ds, ArchGDAL.getgeotransform(dataset))
   ArchGDAL.setproj!(ds, ArchGDAL.getproj(dataset))
   ArchGDAL.getband(ds, 1) do band
      ArchGDAL.setnodatavalue!(band, typemin(Int32))
      ArchGDAL.fillraster!(band, 0)
   end
   ds
end


hrsl_path = "hrsl_ls_adjusted.tif"
admin3_path = expanduser(joinpath("~", "data", "projects", "uganda2020", "outputs",
   "uganda_subcounties_2019_topology_fix", "uganda_subcounties_2019_topology_fix.shp"))
admin1_path = expanduser(joinpath("~", "data", "projects", "uganda2020", "inputs",
   "uganda_districts_2019-wgs84", "uganda_districts_2019-wgs84.shp"))
@show ispath(admin1_path)
ArchGDAL.read(hrsl_path) do hrsl_ds
   geo = band_geo(ArchGDAL.getgeotransform(hrsl_ds))
   hrsl_band = ArchGDAL.getband(hrsl_ds, 1)
   hrsl = ArchGDAL.read(hrsl_band)
   hrsl_no_data = ArchGDAL.getnodatavalue(hrsl_band)
   hrsl_raster = RasterRead(hrsl, hrsl_no_data, geo)
   count_up = build_countup(hrsl_raster)

   ArchGDAL.read(admin1_path) do sc_dataset
      @show ArchGDAL.nlayer(sc_dataset)
      admin = ArchGDAL.getlayer(sc_dataset, 0)
      layerinfo(admin)
      totals = overgeometry(count_up, admin)
      writedlm("admin1_largest_fraction.txt", totals)
      largest = maximum(x[2] for x in totals)
      @show largest

      copy_path = "hrsl_proportion.tif"
      if ispath(copy_path)
         rm(copy_path)
      end
      ArchGDAL.copy(hrsl_ds, filename = copy_path) do hrsl_copy_ds
         ArchGDAL.write!(hrsl_copy_ds, hrsl, 0)
      end
   end
end
