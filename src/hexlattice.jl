
"""
  hexlattice(width, height, r)

Build a 2D simplex from a hexagonal grid.

The grid will cover the given width and height. Each edge in the grid
will be of length `r`. This returns a tuple with vertices as a Float64
array of `(2, n)` and triangles as an Int64 array of `(3, m)`, where
triangle indices point into the vertices.
"""
function hexlattice(width, height, r)
    # Us a Bravais lattice.
    # Two vectors: r1, r2, of the same length.
    # The angle between them is 120 degrees.
    #
    # a1 = (r, 0)
    # a2 = (-r/sqrt(3), 2r/sqrt(3))
    # x,y = i a1 + j a2
    # x = ri - rj/2
    # y = 2jr/√3
    # m and n are the row count (x) and column count (y).
    m, n = Int64.(ceil.([width / r, height * sqrt(3) / 2]))
    # This is all about making a square region from a Bravais
    # lattice that is tilted. We do that by coming up with a rule
    # for which entries are saved.
    isvertex((i, j)) = (0 <= (i-1) - (j-1) ÷ 2 < m) && (0 <= j-1 < n)

    v = Dict{Tuple{Int64,Int64},Int64}()
    vidx = 1
    for j in 1:n
        ioffset = (j - 1) ÷ 2
        for i in (1 + ioffset):(m + ioffset)
            @assert isvertex((i, j))
            v[(i, j)] = vidx
            vidx += 1
        end
    end
    vertex = zeros(Float64, 2, length(v))
    for ((i, j), idx) in v
        vertex[:, idx] = [r * (i-1) - r * (j-1) / 2, 2 * r * (j-1) / 3]
    end

    # This is too many but we can trim at the end, faster than adding to array.
    triangles_pad = zeros(Int64, 3, length(v))
    triangle_cnt = 0
    for (i, j) in keys(v)
        # Each unit of the Bravais lattice is a parallelogram, so two triangles.
        triangle = [
                [(i, j), (i + 1, j), (i + 1, j + 1)]
                [(i, j), (i + 1, j + 1), (i, j + 1)]
                ]
        for t in 1:2
            if all(isvertex.(triangle[5]))
                triangle_cnt += 1
                triangles_pad[:, triangle_cnt] = [v[t] for t in triangle[t]]
            end
        end
    end
    triangles = triangles_pad[:, 1:triangle_cnt]
    vertex, triangles
end
