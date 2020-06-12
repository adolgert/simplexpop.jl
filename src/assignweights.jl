using ArchGDAL

"""
    geo_to_xy(geo, i, j)

i and j are 1-based.

In case of north up images, the GT(2) and GT(4) coefficients are zero, and the
GT(1) is pixel width, and GT(5) is pixel height. The (GT(0),GT(3)) position is
the top left corner of the top left pixel of the raster.

Note that the pixel/line coordinates in the above are from (0.0,0.0) at the top
left corner of the top left pixel to (width_in_pixels,height_in_pixels) at the
bottom right corner of the bottom right pixel. The pixel/line location of the
center of the top left pixel would therefore be (0.5,0.5).
"""
geo_to_xy_corner(geo, i, j) = (
        geo[1] + (i - 1) * geo[2] + (j - 1) * geo[3],
        geo[4] + (i - 1) * geo[5] + (j - 1) * geo[6]
        )

geo_to_xy_center(geo, i, j) = (
        geo[1] + (i - .5) * geo[2] + (j - .5) * geo[3],
        geo[4] + (i - .5) * geo[5] + (j - .5) * geo[6]
        )


function xy_bounds(geo, ij_dim)
    ul = geo_to_xy_corner(geo, 1, 1)
    lr = geo_to_xy_corner(geo, ij_dim[1] + 1, ij_dim[2] + 1)
    ul, lr
end


"""
Pixels are numbered from 1.
"""
function pixel_containing(geo, xy)
    @assert geo[3] == 0
    @assert geo[5] == 0

    (1 + Int(floor(1.0 / geo[2] * (xy[1] - geo[1]))),
     1 + Int(floor(1.0 / geo[6] * (xy[2] - geo[4]))))
end


"""
The rect is a tuple of (ul, lr).
"""
function ij_cover_rect(geo, rect)
    ul_ij = pixel_containing(geo, rect[1])
    lr_ij = pixel_containing(geo, rect[2])
    ul_ij, lr_ij
end


function block_cover_ij_rect(ij_rect, blocksize)
    ((1 + (ij_rect[1][1] - 1) ÷ blocksize[1], 1 + (ij_rect[1][2] - 1) ÷ blocksize[2]),
     (1 + (ij_rect[2][1] - 1) ÷ blocksize[1], 1 + (ij_rect[2][2] - 1) ÷ blocksize[2]))
end


"""
Given HRSL on 30m grid and LandScan (LS) on 900m grid,
transfer landscan density to HRSL settlement pixels that are nonzero.

1. Find bounds on points b/c that's Uganda.
2. Find ls rect that covers those points.
3. Find block rect that covers the ls rect.
4. For each block in ls.
   1. Load block.
   2. For each pixel inside ls rect that covers points.
      1. Get points rect within pixel. (some pixels only partially covered)
      2. Load points blocks that cover that rect.
      3. Sum number of nonzero points within rect.
      4. Rescale by pixel value.
      5. Write modified points blocks.
"""
function assignweights(points_path, weight_path, reweighted_path)
    # ds = dataset
    fine_copy_ds = ArchGDAL.read(expanduser(points_path)) do fine_ds
        reweighted_path = expanduser(reweighted_path)
        ArchGDAL.copy(fine_ds, filename = reweighted_path)
    end
    ArchGDAL.read(expanduser(weight_path)) do coarse_ds
        fine_geo = ArchGDAL.getgeotransform(fine_copy_ds)
        fine_band = ArchGDAL.getband(fine_copy_ds, 1)

        coarse_geo = ArchGDAL.getgeotransform(coarse_ds)
        coarse_band = ArchGDAL.getband(coarse_ds, 1)

        w = ArchGDAL.width(fine_band)
        h = ArchGDAL.height(fine_band)
        fine_xy_bounds = xy_bounds(fine_geo, (w, h))
        coarse_ij_bounds = ij_cover_rect(coarse_geo, fine_xy_bounds)
        coarse_block_bounds = block_cover_ij_rect(coarse_ij_bounds, ArchGDAL.blocksize(coarse_band))
        block_indices = CartesianIndices((coarse_block_bounds[1][1]:coarse_block_bounds[2][1],
                coarse_block_bounds[1][2]:coarse_block_bounds[2][2]))
    end
end
