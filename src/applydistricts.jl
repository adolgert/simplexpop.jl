using ArchGDAL
using DelimitedFiles
using DataFrames
import GDAL
# You need to load rastersweep.jl before running this.

"""
Gets the canonical names and ids of each admin unit.
```julia
admins = unique(df[!, ["admin1_name", "admin1_id"]])
```
"""
function uganda_canonical()
   csv_path = expanduser(joinpath(
      "~", "data", "projects", "uganda2020", "outputs", "uga_canonical",
      "uga_canonical_names.csv"))
   csv, header = readdlm(csv_path, ',', '\n'; header = true)
   df = DataFrame(csv, reshape(header, length(header)))
   for scol in ["admin1_name", "admin2_name", "admin3_name"]
      df[scol] = String.(df[scol])
   end
   for icol in ["admin1_id", "admin2_id", "admin3_id"]
      df[icol] = Int.(df[icol])
   end
   df
end


function build_identify_district(rasterlist)
   function identify_district(geometry, featuredict)::UInt16
      admin_cnt = length(rasterlist)
      admin = [
         UInt16(featuredict["admin1id"]),
         UInt16(featuredict["admin2id"]),
         UInt16(featuredict["admin3id"])
         ]
      envelope = ArchGDAL.envelope(geometry)
      # upper-left on the left column. lower-right in the right column.
      rect = [envelope.MinX envelope.MaxX; envelope.MaxY envelope.MinY]
      raster = rasterlist[1]
      indices = cover_rect(rect, raster.geo, size(raster.data))
      total = 0.0
      # We use a GDAL point to test for inside polygon.
      geom_point = ArchGDAL.createpoint()
      ArchGDAL.addpoint!(geom_point, 0.0, 0.0)
      for idx in indices
         xy = pixel_center(raster.geo, [idx[1], idx[2]])
         ArchGDAL.setpoint!(geom_point, 0, xy[1], xy[2])
         if ArchGDAL.within(geom_point, geometry)
            for aidx in 1:admin_cnt
               rasterlist[aidx].data[idx] = UInt16(admin[aidx])
            end
         end
      end
      admin[3]
   end
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


function featureinfo(featuredefn, feature)
   nfield = ArchGDAL.nfield(featuredefn)
   values = Dict{String,Any}()
   for fieldidx in 1:nfield
      field = ArchGDAL.getfielddefn(featuredefn, fieldidx - 1)
      name = ArchGDAL.getname(field)
      value = ArchGDAL.getfield(feature, fieldidx - 1)
      values[name] = value
   end
   values
end


function overgeometry(toapply, layer; limit = typemax(Int))
   featuredefn = ArchGDAL.layerdefn(layer)
   nfeature = min(ArchGDAL.nfeature(layer), limit)
   result = nothing
   for featureidx in 1:nfeature
      ArchGDAL.getfeature(layer, featureidx - 1) do feature
         infodict = featureinfo(featuredefn, feature)
         geom = ArchGDAL.getgeom(feature, 0)
         envelope = ArchGDAL.envelope(geom)
         status = toapply(geom, infodict)
         if result === nothing
            result = fill(("", status), nfeature)
         end
         admin = infodict["admin3name"]
         @show admin
         result[featureidx] = (admin, status)
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
   function count_up(geometry, featuredict)::Float64
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
Use 0 as the nodata value.
"""
function create_fraction(dataset, path, admin_cnt = 3)
   b = ArchGDAL.getband(dataset, 1)
   w = ArchGDAL.width(b)
   h = ArchGDAL.height(b)
   driver = ArchGDAL.getdriver("GTiff")
   options = [
      "TILED=YES",
      "BLOCKXSIZE=256",
      "BLOCKYSIZE=256",
      "COMPRESS=LZW",
   ]
   ds = ArchGDAL.create(
         path;
         driver = driver,
         width = w,
         height = h,
         nbands = admin_cnt,
         dtype = UInt16,
         options = options,
         )
   ArchGDAL.setgeotransform!(ds, ArchGDAL.getgeotransform(dataset))
   ArchGDAL.setproj!(ds, ArchGDAL.getproj(dataset))
   for bandidx in 1:admin_cnt
      ArchGDAL.getband(ds, bandidx) do band
         ArchGDAL.setnodatavalue!(band, 0)
         ArchGDAL.fillraster!(band, 0)
         ArchGDAL.setcolorinterp!(band, GDAL.GCI_GrayIndex)
         access = ArchGDAL.accessflag(band)
         @assert access == GDAL.GA_Update
      end
   end
   ds
end


"""
Make a TIFF with a band for each admin unit, where the values are the index
of that admin in the canonical list of admins.
```R
admin2id = raster::raster("~/.julia/dev/simplexpop/hrsl_admin_id.tif", band = 2)
```
"""
function make_index_tiff(; limit = typemax(Int))
   hrsl_path = "hrsl_ls_adjusted.tif"
   admin3_path = expanduser(joinpath("~", "data", "projects", "uganda2020", "outputs",
      "uga_canonical", "uga_canonical_names.shp"))
   copy_path = abspath("hrsl_admin_id.tif")
   @show copy_path
   admin_cnt = 3
   if ispath(copy_path)
      println("$(copy_path) exists, deleting.")
      rm(copy_path)
   end
   id_raster_ds = ArchGDAL.read(hrsl_path) do hrsl_ds
      create_fraction(hrsl_ds, copy_path, admin_cnt)
   end
   geo = band_geo(ArchGDAL.getgeotransform(id_raster_ds))
   hrsl_no_data = 0
   hrsl_raster = [
      RasterRead(ArchGDAL.read(ArchGDAL.getband(id_raster_ds, i)), hrsl_no_data, geo)
      for i in 1:admin_cnt
      ]
   identify_district = build_identify_district(hrsl_raster)

   ArchGDAL.read(admin3_path) do sc_dataset
      @show ArchGDAL.nlayer(sc_dataset)
      admin = ArchGDAL.getlayer(sc_dataset, 0)
      layerinfo(admin)
      totals = overgeometry(identify_district, admin; limit = limit)
      for aidx in 1:admin_cnt
         ArchGDAL.write!(id_raster_ds, hrsl_raster[aidx].data, aidx)
      end
   end
   ArchGDAL.destroy(id_raster_ds)
   ispath(copy_path)
end


function assign_to_admin_unit()
   hrsl_path = "hrsl_ls_adjusted.tif"
   admin3_path = expanduser(joinpath("~", "data", "projects", "uganda2020", "outputs",
      "uga_canonical", "uga_canonical_names.shp"))
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
end
