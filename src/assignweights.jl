using ArchGDAL: AbstractRasterBand, getdataset, getgeotransform, blocksize, pixeltype,
        readblock!, writeblock!, getband, getgeotransform, width, height, getnodatavalue
using ArchGDAL
using OffsetArrays


"""
Find the origin and linear transform from pixel to x-y.
"""
geo_to_transform(geo) = ([geo[1], geo[4]], [geo[2] geo[3]; geo[5] geo[6]])


"""
    geo_to_xy_corner(origin, transform, i, j)

i and j are 1-based.

In case of north up images, the GT(2) and GT(4) coefficients are zero, and the
GT(1) is pixel width, and GT(5) is pixel height. The (GT(0),GT(3)) position is
the top left corner of the top left pixel of the raster.

Note that the pixel/line coordinates in the above are from (0.0,0.0) at the top
left corner of the top left pixel to (width_in_pixels,height_in_pixels) at the
bottom right corner of the bottom right pixel. The pixel/line location of the
center of the top left pixel would therefore be (0.5,0.5).
"""
geo_to_xy_corner(origin, transform, i, j) = transform * [i - 1, j - 1] + origin


"""
    geo_to_xy_center(origin, transform, i, j)

The center of the pixel. See also: [`geo_to_xy_corner`](@ref).
"""
geo_to_xy_center(origin, transform, i, j) = (
        transform * [i - 0.5, j - 0.5] + origin
        )


"""
Which pixel contains this xy point.
"""
function pixel_containing(origin, transform, xy)
    Int.(floor.(transform \ (xy - origin))) .+ 1  # +1 for one-based indexing.
end


"""
Which pixel is closest to this xy point.
"""
function nearest_pixel(origin, transform, xy)
    Int.(round.(transform / (xy - origin))) .+ 1  # +1 for one-based indexing.
end


"""
The rect is a tuple of (ul, lr).
This tends to pull in the next pixel when
the grids are exactly aligned.
"""
function ij_cover_rect(origin, transform, rect)
    ul_ij = pixel_containing(origin, transform, rect[1])
    lr_ij = pixel_containing(origin, transform, rect[2])
    ul_ij, lr_ij
end


function pointinrect(ij, rect)
    rect[1][1] <= ij[1] < rect[2][1] && rect[1][2] <= ij[2] < rect[2][2]
end


function block_cover_ij_rect(ij_rect, blocksize)
    ((1 + (ij_rect[1][1] - 1) ÷ blocksize[1], 1 + (ij_rect[1][2] - 1) ÷ blocksize[2]),
     (1 + (ij_rect[2][1] - 1) ÷ blocksize[1], 1 + (ij_rect[2][2] - 1) ÷ blocksize[2]))
end


function corner_pixel(block_idx::CartesianIndex, blocksize)
    CartesianIndex((block_idx[1] - 1) * blocksize[1] + 1, (block_idx[2] - 1) * blocksize[2] + 1)
end


### PixelGrid methods.

"""
Represents a grid extent and its geometry, without data.

For GDAL, coordinates are (x-offset, y-offset), and it's stored as
(row, col) in the array. We think of the x-offset, visually, as the column,
but forget that. This world is (x-offset, y-offset).

The grid extent is an upper-lefthand `(i0, j0)` and a lower right-hand
`(i1, j1)`, where `i0 < i1` and `j0 < j1`. This orientation comes from
the geo tradition to flip the y-axis, so increasing y goes negative.
The `origin` and `transform` are from the GDAL geo array and represent
the origin of the grid and how each point offsets that origin.
From the GDAL documentation, the geo is used this way.
```
Xgeo = GT(0) + Xpixel*GT(1) + Yline*GT(2)
Ygeo = GT(3) + Xpixel*GT(4) + Yline*GT(5)
```
"""
struct PixelGrid
    ul::Tuple{Int64,Int64}  # As x, y
    lr::Tuple{Int64,Int64}  # As x, y
    origin::Array{Float64,1}  # x, y, too.
    transform::Array{Float64,2}
end


"""
Make a function to get the geotransform so that we can stub the Band during testing.
"""
ArchGDAL.getgeotransform(band::AbstractRasterBand) = getgeotransform(getdataset(band))

"""
Make a PixelGrid that represents all data available in a Band.
"""
function pixelgrid(band::AbstractRasterBand)
    geo = getgeotransform(band)
    w = width(band)
    h = height(band)
    origin, transform = geo_to_transform(geo)
    PixelGrid((1, 1), (w, h), origin, transform)
end


"""
Compute the area of the xy intersection of two pixelgrids.

If the returned value is negative, they don't intersect.
"""
function intersection_area(a::PixelGrid, b::PixelGrid)
    # The +1 is because this gives the same upper-left corner in both cases.
    ula = geo_to_xy_corner(a.origin, a.transform, a.ul[1], a.ul[2])
    lra = geo_to_xy_corner(a.origin, a.transform, a.lr[1] + 1, a.lr[2] + 1)
    ulb = geo_to_xy_corner(b.origin, b.transform, b.ul[1], b.ul[2])
    lrb = geo_to_xy_corner(b.origin, b.transform, b.lr[1] + 1, b.lr[2] + 1)
    ul = [max(ula[1], ulb[1]), max(ula[2], ulb[2])]
    lr = [min(ula[1], ulb[1]), min(ula[2], ulb[2])]
    (lr[1]-ul[1]) * (lr[2] - ul[2])
end


function xy_bounds(pg::PixelGrid)
    ul = geo_to_xy_corner(pg.origin, pg.transform, pg.ul[1], pg.ul[2])
    lr = geo_to_xy_corner(pg.origin, pg.transform, pg.lr[1] + 1, pg.lr[2] + 1)
    ul, lr
end


function crop_to(large::PixelGrid, small::PixelGrid)
    ul_xy, lr_xy = xy_bounds(small)
    ul_ij, lr_ij = ij_cover_rect(large.origin, large.transform, (ul_xy, lr_xy))
    uln = (max(large.ul[1], ul_ij[1]), max(large.ul[2], ul_ij[2]))
    lrn = (min(large.lr[1], lr_ij[1]), min(large.lr[2], lr_ij[2]))
    PixelGrid(uln, lrn, large.origin, large.transform)
end


invalid(pg::PixelGrid) = (pg.ul[1] > pg.lr[1]) || (pg.ul[2] > pg.lr[2])
valid(pg::PixelGrid) = !invalid(pg)

"""
Given a pixel grid, pull out a pixel grid for a single pixel.
"""
function onepixelgrid(pg::PixelGrid, pij)
    PixelGrid(
            (pij[1], pij[2]),
            (pij[1], pij[2]),
            pg.origin,
            pg.transform
            )
end



"""
Iterate over all blocks that cover a given pixel rectangle.
"""
function iter_block_indices(pg::PixelGrid, blocksize)
    coarse_block_bounds = block_cover_ij_rect((pg.ul, pg.lr), blocksize)
    i, j = 1, 2
    CartesianIndices((coarse_block_bounds[1][i]:coarse_block_bounds[2][i],
            coarse_block_bounds[1][j]:coarse_block_bounds[2][j]))
end


### Stub of GDAL's AbstractRasterBand, for testing.

struct StubBand{T} <: AbstractRasterBand
    A::Array{T, 2}  # data
    bs::Array{Int32,1}  # blocksize
    o::Array{Float64,1}  # origin in xy.
    t::Array{Float64,2}  # linear transform from pixel to xy.
end

ArchGDAL.width(b::StubBand) = size(b.A, 1)
ArchGDAL.height(b::StubBand) = size(b.A, 2)
ArchGDAL.blocksize(b::StubBand) = b.bs
ArchGDAL.pixeltype(b::StubBand{T}) where {T} = T
function ArchGDAL.readblock!(b::StubBand, i::Integer, j::Integer, buffer)
    @assert 0 <= b.bs[1] * i < width(b)  # beginning of block is less than width.
    @assert 0 <= b.bs[2] * j < height(b)
    # The last block will likely extend past the end of the actual data.
    xaxis = (b.bs[1]*i + 1):min(b.bs[1]*(i+1), size(b.A, 1))
    yaxis = (b.bs[2]*j + 1):min(b.bs[2]*(j+1), size(b.A, 2))
    buffer[1:length(xaxis), 1:length(yaxis)] = b.A[xaxis, yaxis]
end
ArchGDAL.getnodatavalue(b::StubBand{T}) where {T <: Real} = floatmin(T)
ArchGDAL.getnodatavalue(b::StubBand{T}) where {T <: Integer} = typemin(T)
function ArchGDAL.writeblock!(b::StubBand, i::Integer, j::Integer, buffer)
    @assert 0 <= b.bs[1] * i < width(b)  # beginning of block is less than width.
    @assert 0 <= b.bs[2] * j < height(b)
    # The last block will likely extend past the end of the actual data.
    xaxis = (b.bs[1]*i + 1):min(b.bs[1]*(i+1), size(b.A, 1))
    yaxis = (b.bs[2]*j + 1):min(b.bs[2]*(j+1), size(b.A, 2))
    b.A[xaxis, yaxis] = buffer[1:length(xaxis), 1:length(yaxis)]
end
ArchGDAL.getgeotransform(b::StubBand) = [b.o[1], b.t[1,1], b.t[1,2], b.o[2], b.t[2,1], b.t[2,2]]


### DataGrid methods.
"""
A DataGrid represents a matrix of raster values.

The data is stored as an array of data blocks, but its extent is usually smaller than
the extent of the blocks. We can always assume that the data is complete blocks.
We keep the data, A, indexed in absolute ids, so the axis range starts with the actual
values, not usually from 1.
"""
struct DataGrid{T}
    A::AbstractArray{T,2}
    nodatavalue::T
    blocksize::Array{Int32,1}
    ulb::Tuple{Int32,Int32}  # upper left block
    lrb::Tuple{Int32,Int32}  # lower right block
    pg::PixelGrid
end


function load_data_grid(band::AbstractRasterBand, pg::PixelGrid)
    bs = blocksize(band)
    ulb, lrb = block_cover_ij_rect((pg.ul, pg.lr), bs)
    dtype = pixeltype(band)
    ul_ij = ((ulb[1] - 1) * bs[1] + 1, (ulb[2] - 1) * bs[2] + 1)
    lr_ij = (lrb[1] * bs[1], lrb[2] * bs[2])
    @show ul_ij
    @show lr_ij
    A = zeros(dtype, lr_ij[1] - ul_ij[1] + 1, lr_ij[2] - ul_ij[2] + 1)
    buffer = zeros(dtype, bs...)
    for bidx in CartesianIndices(((lrb[1] - ulb[1] + 1), (lrb[2] - ulb[2] + 1)))
        # readblock! is zero-based and is col, row b/c it's xaxis, yaxis.
        readblock!(band, bidx[1] - 1, bidx[2] - 1, buffer)
        # Blocksize is still x, y, so it's flipped here.
        A[((bidx[1] - 1) * bs[1] + 1):(bidx[1] * bs[1]),
            ((bidx[2] - 1) * bs[2] + 1):(bidx[2] * bs[2])] = buffer
    end
    offset = OffsetArray(A, ul_ij[1]:lr_ij[1], ul_ij[2]:lr_ij[2])
    nodata = getnodatavalue(band)
    DataGrid{dtype}(offset, nodata, bs, ulb, lrb, pg)
end


"""
Count nonzero pixels that are inside the boundsrect.
"""
function count_nonzero_pixels(dg::DataGrid, boundsrect)
    cnt = 0
    ul = dg.pg.ul
    lr = dg.pg.lr
    # i, j are x, y for GDAL but flipped for the array.
    for j in ul[2]:lr[2]
        for i in ul[1]:lr[1]
            if dg.A[i, j] != dg.nodatavalue
                pt = geo_to_xy_center(dg.pg.origin, dg.pg.transform, i, j)
                if pointinrect(pt, boundsrect)
                    cnt += 1
                end
            end
        end
    end
    percent_filled = cnt / ((lr[1] - ul[1]) * (lr[2] - ul[2]))
    @show percent_filled
    cnt
end


function fill_nonzero_pixels!(dg::DataGrid, value)
    ul = dg.pg.ul
    lr = dg.pg.lr
    for j in ul[2]:lr[2]
        for i in ul[1]:lr[1]
            if dg.A[i, j] != dg.nodatavalue
                dg.A[i, j] = value
            end
        end
    end
end


function write_datagrid!(dg::DataGrid, band)
    bs = dg.blocksize
    for j in dg.ulb[2]:dg.lrb[2]
        for i in dg.ulb[1]:dg.lrb[1]
            xrange = ((i - 1) * bs[1] + 1):(i * bs[1])
            yrange = ((j - 1) * bs[2] + 1):(j * bs[2])
            writeblock!(band, i-1, j-1, dg.A[xrange, yrange])
        end
    end
end


function load_single_block(band::AbstractRasterBand, pg::PixelGrid, block)
    bs = blocksize(band)
    @show bs
    dtype = pixeltype(band)
    A = zeros(dtype, bs[1], bs[2])
    readblock!(band, block[1] - 1, block[2] - 1, A)
    ul = ((block[1] - 1) * bs[1] + 1, (block[2] - 1) * bs[2] + 1)
    lr = (block[1] * bs[1], block[2] * bs[2])
    offset = OffsetArray(A, ul[1]:lr[1], ul[2]:lr[2])
    npg = PixelGrid(ul, lr, pg.origin, pg.transform)
    nodata = getnodatavalue(band)
    DataGrid{dtype}(offset, nodata, bs, Tuple(block), Tuple(block), npg)
end


function tune_band(fine_band, coarse_band, max_pixels)
    fine_pg = pixelgrid(fine_band)
    coarse_pg = pixelgrid(coarse_band)
    coarse_crop_pg = crop_to(coarse_pg, fine_pg)
    coarse_blocksize = blocksize(coarse_band)

    outside_cnt = 0
    coarse_buffer = zeros(Int32, coarse_blocksize...)
    pixel_cnt = 0
    for block_idx in iter_block_indices(coarse_crop_pg, coarse_blocksize)
        block = load_single_block(coarse_band, coarse_pg, block_idx)
        cartesian = CartesianIndices(block.A)
        for linear_index in eachindex(block.A)
            pij = cartesian[linear_index]
            single_coarse_pg = onepixelgrid(coarse_crop_pg, pij)
            fine_cover_pg = crop_to(fine_pg, single_coarse_pg)
            if valid(fine_cover_pg)
                fine_grid = load_data_grid(fine_band, fine_cover_pg)
                pixel_cnt += 1
                area = intersection_area(single_coarse_pg, fine_cover_pg)
                @assert area >= 0
                landscan_value = block.A[
                    single_coarse_pg.ul[1],
                    single_coarse_pg.ul[2]
                ]
                # This should be restricted to pixels under the coarse.
                pixel_cnt = count_nonzero_pixels(fine_grid, xy_bounds(single_coarse_pg))
                @show pixel_cnt
                rate = min(landscan_value / pixel_cnt, 100)
                fill_nonzero_pixels!(fine_grid, rate)
                write_datagrid!(fine_grid, fine_band)
            else
                outside_cnt += 1
            end
            if pixel_cnt > max_pixels
                return nothing
            end
        end
    end
    @show outside_cnt
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
function assignweights(
        points_path, weight_path, reweighted_path;
        max_pixels = typemax(Int64)
        )
    # ds = dataset
    fine_copy_ds = ArchGDAL.read(expanduser(points_path)) do fine_ds
        reweighted_path = expanduser(reweighted_path)
        ArchGDAL.copy(fine_ds, filename = reweighted_path)
    end
    ArchGDAL.read(expanduser(weight_path)) do coarse_ds
        fine_band = getband(fine_copy_ds, 1)
        coarse_band = getband(coarse_ds, 1)
        tune_band(fine_band, coarse_band, max_pixels)
    end
end
