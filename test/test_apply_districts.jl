
function copy_hrsl()
   hrsl_path = "hrsl_ls_adjusted.tif"
   admin3_path = expanduser(joinpath("~", "data", "projects", "uganda2020", "outputs",
      "uga_canonical", "uga_canonical_names.shp"))
   copy_path = "ztesting.tif"
   admin_cnt = 3
   if ispath(copy_path)
      rm(copy_path)
   end
   id_raster_ds = ArchGDAL.read(hrsl_path) do hrsl_ds
      create_fraction(hrsl_ds, copy_path, admin_cnt)
   end
   ArchGDAL.destroy(id_raster_ds)
end

copy_hrsl()
@test ispath("ztesting.tif")
