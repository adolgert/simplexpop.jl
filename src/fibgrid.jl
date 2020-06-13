floatmod(x, y) = x - y * floor(x / y)

"""
Fibonacci Grid on a Sphere

There are 2n + 1 points on the sphere.
It returns a 2xn grid of lat-long.

Measurement of areas on a sphere using Fibonacci and latitude–longitude lattices
 ́Alvaro Gonz ́alez.

This makes a grid. Then make an r^2 potential among
neighboring points and use successive over-relaxation
to make final adjustments.
"""
function fibonacci_points(n)
    ϕ = (1 + sqrt(5)) / 2
    points = zeros(Float64, 2, 2n + 1)
    for i in -n:n
        points[1, i + n + 1] = asin(2i / (2n + 1)) * 180 / π
        lon = floatmod(i, ϕ) * 360 / ϕ
        if lon < -180
            points[2, i + n + 1] = lon + 360
        elseif lon > 180
            points[2, i + n + 1] = lon - 360
        else
            points[2, i + n + 1] = lon
        end
    end
    points
end


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
    triangles = zeros(Int64, 3, length(v))
    triangle_idx = 1
    for (i, j) in keys(v)
        triangle = [(i, j), (i + 1, j), (i + 1, j + 1)]
        if all(isvertex.(triangle))
            triangles[:, triangle_idx] = [v[t] for t in triangle]
            triangle_idx += 1
        end
    end
    triangle_cnt = triangle_idx - 1
    vertex, triangles[:, 1:triangle_cnt]
end
