
"""
  hexlattice(width, height, r)

Hexagonal grid built from a Bravais lattice.
Two vectors: r1, r2, of the same length.
The angle between them is 120 degrees.

a1 = (r, 0)
a2 = (-r/sqrt(3), 2r/sqrt(3))
x,y = i a1 + j a2
x = ri - rj/2
y = 2jr/√3
"""
function hexlattice(width, height, r)
    # Make rows and columns that fit.
    m, n = Int64.(floor.([width / r, height * sqrt(3) / 2]))
    @show (m, n)
    # This is all about making a square region from a Bravais
    # lattice that is tilted. We do that by coming up with a rule
    # for which entries are saved.
    isvertex((i, j)) = (0 <= (i-1) - (j-1) ÷ 2 < m) && (0 <= j-1 < n)

    v = Dict{Tuple{Int64,Int64},Int64}()
    vidx = 1
    for j in 1:n
        ioffset = (j - 1) ÷ 2
        for i in (1 + ioffset):(m + ioffset)
            if !isvertex((i, j))
                @show (i, j)
            end
            v[(i, j)] = vidx
            vidx += 1
        end
    end
    vertex = zeros(Float64, 2, length(v))
    for ((i, j), idx) in v
        vertex[:, idx] = [r * (i-1) - r * (j-1) / 2, 2 * r * (j-1) / 3]
    end
    # This is too many but we can trim at the end
    triangles_pad = zeros(Int64, 3, length(v))
    triangle_idx = 1
    for (i, j) in keys(v)
        triangle = [
                [(i, j), (i + 1, j), (i + 1, j + 1)]
                [(i, j), (i + 1, j + 1), (i, j + 1)]
                ]
        for t in 1:2
            if all(isvertex.(triangle[5]))
                triangles_pad[:, triangle_idx] = [v[t] for t in triangle[t]]
                triangle_idx += 1
            end
        end
    end
    triangle_cnt = triangle_idx - 1
    triangles = triangles_pad[:, 1:triangle_cnt]
end
