using ArchGDAL
using StaticArrays

struct Geo
    origin::SVector{2,Float64}
    transform::SMatrix{2, 2, Float64}
end

"""
Find the origin and linear transform from pixel to x-y.
"""
function band_geo(band::ArchGDAL.AbstractRasterBand)
    geo = ArchGDAL.getgeotransform(band)
    Geo([geo[1], geo[4]], [geo[2] geo[3]; geo[5] geo[6]])
end

pixel_corner(geo::Geo, ij) = geo.transform * (ij .- 1) + geo.origin
pixel_center(geo::Geo, ij) = geo.transform * (ij .- 0.5) + geo.origin
# First column is upper-left (x y). Second column is lower-right (x, y).
raster_rect(geo::Geo, ullr) = pixel_corner(geo, hcat(ullr[:, 1], ullr[:, 2] .+ 1))
partial_pixel_containing(geo::Geo, xy) = (geo.transform \ (xy - origin)) .+ 1
offset_geo(geo::Geo, offset) = Geo(pixel_corner(geo, offset), geo.transform)

function ij_cover_rect(geo::Geo, rect; eps = 1e-7)
    ul = partial_pixel_containing(geo, rect[:, 1])
    ceil_dx = ceil.[ul] - ul
    ulij = Int32.(floor.(ul))
    for uli in 1:2
        if ceil_dx[uli] < eps
            ulij[uli] = ceil_dx[uli]
        end
    end
    lr = partial_pixel_containing(geo, rect[:, 2])
    floor_dx = floor.[lr] - lr
    lrij = Int32.(ceil.(lr))
    for lri in 1:2
        if floor_dx[lri] < eps
            lrij[lri] = floor_dx[lri]
        end
    end
    hcat(ulij, lrij)
end


function block_cover_rect(ul_lr, bounds, blocksize)
    # available bounds may be less than what covers this rectangle.
    ula = [  # first column is upper-left. second column is lower-right.
        max(ul_lr[1, 1], bounds[1][1])  min(ul_lr[1, 2], bounds[1][end])
        max(ul_lr[2, 1], bounds[2][1])  min(ul_lr[2, 2], bounds[2][end])
    ]
    # return as two UnitRanges, one for each axis.
    toblock(ij, blocksize) = (ij .- 1) ÷ blocksize .+ 1
    ul = toblock(ula[:, 1])
    lr = toblock(ula[:, 2])
    (ul[1]:lr[1], ul[2]:lr[2])
end


struct RasterBlock{T}
    A::Array{T}
    range::Tuple{UnitRange{Int64},UnitRange{Int64}}
    geo::Geo
    nodatavalue::T
end


"""
The processing function works one block at a time and returns whether
    it wrote to a block.
"""
function process_block(blocka::RasterBlock, blockb::RasterBlock)
    for pixela in iterblock(blocka)
        rect = pixel_rect(blocka, pixela)
        cnt = 0
        for pixel in inside_rect(blockb, rect)
            if pixel != blockb.nodatavalue
                cnt += 1
            end
        end
        rate = min(blocka.A[pixela] / cnt, 100.0)
        for pixel in inside_rect(blockb, rect)
            if pixel != blockb.nodatavalue
                blockb[pixel] = rate
            end
        end
    end
    false, true
end


function over_block_pairs(processor::Function, band, geo)
    blocksize = ArchGDAL.blocksize.(band)
    bounds = (bd -> [1:ArchGDAL.width(bd), 1:ArchGDAL.height(bd)]).(band)
    dtype = ArchGDAL.pixeltype.(band)
    # Separate these because they are different types.
    nodata1 = ArchGDAL.getnodatavalue(band1)
    nodata2 = ArchGDAL.getnodatavalue(band2)

    rect2 = raster_rect(geo[2], bounds[2])
    ul_lr1 = ij_cover_rect(geo[1], rect2)
    block_bounds1 = block_cover_rect(ul_lr1, bounds[1], blocksize[1])
    buffer1 = zeros(dtype[1], bounds[1])
    for j in block_bounds1[2]
        for i in block_bounds1[1]
            fill!(buffer1, nodata[1])
            readblock!(band[1], i - 1, j - 1, buffer1)
            block_geo1 = offset_geo(geo[1], blocksize[1] .* ([i, j] .- 1))

            block1_rect = raster_rect(block_geo1, axes(buffer1))
            ul_lr2 = ij_cover_rect(geo[2], block1_rect)
            block_bounds2 = block_cover_rect(ul_lr2, bounds[2], blocksize[2])
            xdim = block_bounds2[1][end] - block_bounds2[1][1]
            ydim = block_bounds2[2][end] - block_bounds2[2][1]
            buffer2 = zeros(dtype[2], xdim, ydim)
            smallbuffer = zeros(dtype[2], blocksize[2]...)
            for j2 in block_bounds2[2]
                for i2 in block_bounds2[1]
                    readblock!(band[2], i2 - 1, j2 - 1, smallbuffer)
                    buffer2[
                        ((i2 - 1) * blocksize[2][1] + 1):(i2 * blocksize[2][1]),
                        ((j2 - 1) * blocksize[2][2] + 1):(j2 * blocksize[2][2]),
                        ] = smallbuffer
                end
            end
            # These are missing the pixel ranges that are allowed.
            rb1 = RasterBlock{dtype[1]}(buffer1, block_geo1)
            rb2 = RasterBlock{dtype[2]}(buffer2, block_geo2)
            towrite = processor(rb1, rb2)
        end
    end
end


function adjust_hrsl_with_landscan(a::ArchGDAL.Dataset, b::ArchGDAL.Dataset)
    band = [ArchGDAL.getband(a, 1), ArchGDAL.getband(b, 1)]
    geo = band_geo.([a, b])
    over_block_pairs(process_block, band, geo)
end
