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
   result = zeros(Base.return_types(toapply)[1], nfeature)
   for featureidx in 1:nfeature
      ArchGDAL.getfeature(layer, featureidx - 1) do feature
         admin1 = ArchGDAL.getfield(feature, 2 - 1)
         admin3 = ArchGDAL.getfield(feature, 6 - 1)
         geom = ArchGDAL.getgeom(feature, 0)
         envelope = ArchGDAL.envelope(geom)
         @show (admin1, admin3, envelope)
         @show (envelope.MinX, envelope.MaxX, envelope.MinY, envelope.MaxY)
         result[featureidx] = toapply(geom)
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
      for idx in indices
         if raster.data[idx] != raster.no_data
            xy = pixel_center(raster.geo, [idx[1], idx[2]])
            ArchGDAL.setpoint!(geom_point, 0, xy[1], xy[2])
            if ArchGDAL.within(geom_point, geometry)
               total += raster.data[idx]
            end
         end
      end
      @show total
      total
   end
end


hrsl_path = "hrsl_ls_adjusted.tif"
admin3_path = expanduser(joinpath("~", "data", "projects", "uganda2020", "outputs",
   "uganda_subcounties_2019_topology_fix", "uganda_subcounties_2019_topology_fix.shp"))
@show ispath(admin3_path)
ArchGDAL.read(hrsl_path) do hrsl_ds
   geo = band_geo(ArchGDAL.getgeotransform(hrsl_ds))
   hrsl_band = ArchGDAL.getband(hrsl_ds, 1)
   hrsl = ArchGDAL.read(hrsl_band)
   hrsl_no_data = ArchGDAL.getnodatavalue(hrsl_band)
   hrsl_raster = RasterRead(hrsl, hrsl_no_data, geo)
   count_up = build_countup(hrsl_raster)

   ArchGDAL.read(admin3_path) do sc_dataset
      @show ArchGDAL.nlayer(sc_dataset)
      admin3 = ArchGDAL.getlayer(sc_dataset, 0)
      layerinfo(admin3)
      totals = overgeometry(count_up, admin3)
      @show sum(totals)
      writedlm("admin3_totals.txt", totals)
   end
end
