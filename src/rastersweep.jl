using ArchGDAL
using StaticArrays


struct Geo
    origin::Array{Float64, 1}
    transform::Array{Float64, 2}
end

"""
Find the origin and linear transform from pixel to x-y.
"""
function band_geo(geo)
    Geo([geo[1], geo[4]], [geo[2] geo[3]; geo[5] geo[6]])
end

pixel_corner(geo::Geo, ij) = geo.transform * (ij .- 1) .+ geo.origin
pixel_center(geo::Geo, ij) = geo.transform * (ij .- 0.5) .+ geo.origin
# First column is upper-left (x y). Second column is lower-right (x, y).
raster_rect(geo::Geo, ullr) = pixel_corner(geo, hcat(ullr[:, 1], ullr[:, 2] .+ 1))
partial_pixel_containing(geo::Geo, xy) = (geo.transform \ (xy - geo.origin)) .+ 1
offset_geo(geo::Geo, offset) = Geo(pixel_corner(geo, offset), geo.transform)
pixel_rect(geo::Geo, index) = pixel_corner(geo, [index[1] (index[1] + 1); index[2] (index[2] + 1)])

function cover_rect(rect, geo, bounds)
    # This gives you floating point pixel that contain the rectangle.
    ul_ij = partial_pixel_containing(geo, rect[:, 1])
    lr_ij = partial_pixel_containing(geo, rect[:, 2])
    # Include an upper-left if the center of the pixel is included.
    # 2.2 is 2. 2.7 is 3.
    ul = Int32.(floor.(ul_ij .+ 0.5))
    # For lower-right, 10.2 is 9. 10.7 is 10.
    lr = Int32.(floor.(lr_ij .- 0.5))
    ul = max.(ul, 1)
    lr = min.(lr, bounds)
    if lr[1] < ul[1] || lr[2] < ul[2]
        return CartesianIndices((0, 0))
    else
        return CartesianIndices((ul[1]:lr[1], ul[2]:lr[2]))
    end
end


function over_sweep!(geo, coarse, fine, fine_missing)
    # fine_rect = raster_rect(geo[2], [1 size(fine, 1); 1 size(fine, 2)])
    for cidx in CartesianIndices(coarse)
        rect = pixel_rect(geo[1], cidx)
        fine_indices = cover_rect(rect, geo[2], size(fine))
        if length(fine_indices) == 0
            continue
        end
        nonzero = sum([fine[fidx] != fine_missing for fidx in fine_indices])
        if nonzero > 0
            rate = min(100, coarse[cidx] / nonzero)
            for rmidx in fine_indices
                if fine[rmidx] != fine_missing
                    fine[rmidx] = rate
                    # else leave it missing
                end
            end
        end
    end
    return nothing
end


function whole_band(a::ArchGDAL.Dataset, b::ArchGDAL.IDataset)
    band = [ArchGDAL.getband(a, 1), ArchGDAL.getband(b, 1)]
    raw_geo = ArchGDAL.getgeotransform.([a, b])
    geo = band_geo.(raw_geo)
    coarse = ArchGDAL.read(band[1])
    fine = ArchGDAL.read(band[2])
    fine_no_data = ArchGDAL.getnodatavalue(band[2])
    over_sweep!(geo, coarse, fine, fine_no_data)
    ArchGDAL.write!(band[2], fine)
end


function assignweights(points_path, weight_path, reweighted_path)
    # ds = dataset
    fine_copy_ds = ArchGDAL.read(expanduser(points_path)) do fine_ds
        reweighted_path = expanduser(reweighted_path)
        ArchGDAL.copy(fine_ds, filename = reweighted_path)
    end
    ArchGDAL.read(expanduser(weight_path)) do coarse_ds
        whole_band(coarse_ds, fine_copy_ds)
    end
    ArchGDAL.destroy(fine_copy_ds)
end
