using ArchGDAL
using OffsetArrays

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


function corner_pixel(block_idx::CartesianIndex, blocksize)
    CartesianIndex((block_idx[1] - 1) * blocksize[1] + 1, (block_idx[2] - 1) * blocksize[2] + 1)
end


function pixels_of_block(block_idx::CartesianIndex, blocksize)
    CartesianIndices((
        ((block_idx[1] - 1) * blocksize[1] + 1):(block_idx[1] * blocksize[1]),
        ((block_idx[2] - 1) * blocksize[2] + 1):(block_idx[2] * blocksize[2])
    ))
end

"""
Combines integer bounds on a grid with its geometry.
"""
struct PixelGrid
    ul::Tuple{Int64,Int64}
    lr::Tuple{Int64,Int64}
    geo::Array{Float64,1}
end


function pixelgrid(band::ArchGDAL.AbstractRasterBand)
    ds = ArchGDAL.getdataset(band)
    geo = ArchGDAL.getgeotransform(ds)
    w = ArchGDAL.width(band)
    h = ArchGDAL.height(band)
    PixelGrid((1, 1), (w, h), geo)
end


function xy_bounds(pg::PixelGrid)
    ul = geo_to_xy_corner(pg.geo, pg.ul[1], pg.ul[2])
    lr = geo_to_xy_corner(pg.geo, pg.lr[1] + 1, pg.lr[2] + 1)
    ul, lr
end


function crop_to(large::PixelGrid, small::PixelGrid)
    ul, lr = xy_bounds(small)
    ul_ij, lr_ij = ij_cover_rect(large.geo, (ul, lr))
    uln = (max(large.ul[1], ul_ij[1]), max(large.ul[2], ul_ij[2])
    lrn = (min(large.lr[1], lr_ij[1]), min(large.lr[2], lr_ij[2])
    PixelGrid(uln, lrn, large.geo)
end


"""
Every raster data is read and written in blocks. That's a rule.
The inner pixel grid tells us what the desired rectangle is.
The array has offset indices from the original large grid.
"""
struct DataGrid{T}
    A::AbstractArray{T,2}
    ulb::Tuple{Int32,Int32}  # upper left block
    lrb::Tuple{Int32,Int32}  # lower right block
    pg::PixelGrid
end


function load_pixel_grid(band::ArchGDAL.AbstractRasterBand, pg::PixelGrid)
    blocksize = ArchGDAL.blocksize(band)
    ulb, lrb = block_cover_ij_rect((pg.ul, pg.lr), blocksize)
    dtype = ArchGDAL.pixeltype(band)
    A = zeros(dtype, (lrb[1] - ulb[1] + 1) * blocksize[1], (lrb[2] - ulb[2] + 1) * blocksize[2])
    buffer = zeros(dtype, blocksize...)
    for bidx in CartesianIndices(((lrb[1] - ulb[1] + 1), (lrb[2] - ulb[2] + 1)))
        # readblock! is zero-based.
        ArchGDAL.readblock!(band, bidx[1] - 1, bidx[2] - 1, buffer)
        A[((bidx[1] - 1) * blocksize[1] + 1):(bidx[1] * blocksize[1]),
            ((bidx[2] - 1) * blocksize[2] + 1):(bidx[2] * blocksize[2])] = buffer
    end
    offset = OffsetArray(
            A,
            ((ulb[1] - 1) * blocksize[1] + 1):(lrb[1] * blocksize[1]),
            ((ulb[2] - 1) * blocksize[2] + 1):(lrb[2] * blocksize[2]),
            )
    DataGrid{dtype}(offset, ulb, lrb, pg)
end


function load_single_block(band::ArchGDAL.AbstractRasterBand, pg::PixelGrid, block)
    blocksize = ArchGDAL.blocksize(band)
    dtype = ArchGDAL.pixeltype(band)
    A = zeros(dtype, blocksize...)
    ArchGDAL.readblock!(band, block[1] - 1, block[2] - 1, A)
    ul = ((block[1] - 1) * blocksize[1] + 1, (block[2] - 1) * blocksize[2] + 1)
    lr = ((block[1] * blocksize[1]), (block[2] * blocksize[2]))
    xaxis = (ul[1]:lr[1])
    yaxis = (ul[2]:lr[2])
    offset = OffsetArray(A, xaxis, yaxis)
    npg = PixelGrid(ul, lr, pg.geo)
    DataGrid{dtype}(offset, Tuple(block), Tuple(block), npg)
end


"""
Given HRSL on 30m grid and LandScan (LS) on 900m grid,
transfer landscan density to HRSL settlement pixels that are nonzero.

1. Find bounds on points b/c that's a country's HRSL.
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
        fine_band = ArchGDAL.getband(fine_copy_ds, 1)
        fine_pg = pixelgrid(fine_band)
        coarse_band = ArchGDAL.getband(coarse_ds, 1)
        coarse_pg = pixelgrid(coarse_band)
        coarse_crop_pg = crop_to(coarse_pg, fine_pg)

        coarse_blocksize = ArchGDAL.blocksize(coarse_band)
        coarse_block_bounds = block_cover_ij_rect(coarse_ij_bounds, coarse_blocksize)
        block_indices = CartesianIndices((coarse_block_bounds[1][1]:coarse_block_bounds[2][1],
                coarse_block_bounds[1][2]:coarse_block_bounds[2][2]))

        coarse_buffer = zeros(Int32, coarse_blocksize...)
        for block_idx in block_indices
            ArchGDAL.readblock!(coarse_band, block_idx[1] - 1, block_idx[2] - 1, coarse_buffer)
            corner = corner_pixel(block_idx, coarse_blocksize)
            for pixel in pixels_of_block(block_idx, coarse_blocksize)
                pixel_rect = xy_bounds(coarse_geo, ((pixel[1], pixel[2]), (pixel[1], pixel[2])))
                work_xy_rect = intersect_xy(pixel_rect, fine_xy_bounds)
                fine_ij_rect = ij_cover_rect(fine_geo, work_xy_rect)
                # Load fine data in such a way that we can write it back with writeblock!
            end
        end
    end
end
