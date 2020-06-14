"""
    floatmod(x, y)

This is the equivalent of `x // y` for floating point.
"""
floatmod(x, y) = x - y * floor(x / y)

"""
    fibonacci_points(n)

Constructs lat-long points on a sphere from a Fibonacci sequence.

There are 2n + 1 points on the sphere. It returns a `(2,2n+1)` grid of lat-long.

These points are a good start on an evenly-spaced set of points on the
surface of a sphere. This algorithm is from González [^1]. The usual way
to use these points is to make them a little more even by adjustment
with a successive over-relaxation algorithm. That algorithm moves them
a little at each step due to a repulsive force among neighbors.
Then you make a Voronoi on the sphere and create triangles from this.
The Voronoi will be almost all hexagons, but there will be a few polygons
with 5 and 7 sides.

[^1]: González, Á. (2010). Measurement of areas on a sphere using Fibonacci
and latitude–longitude lattices. Mathematical Geosciences, 42(1), 49.
"""
function fibonacci_points(n)
    ϕ = (1 + sqrt(5)) / 2
    points = zeros(Float64, 2, 2n + 1)
    for i in -n:n
        points[1, i + n + 1] = asin(2i / (2n + 1)) * 180 / π
        lon = floatmod(i, ϕ) * 360 / ϕ
        if lon < -180
            points[2, i + n + 1] = lon + 360
        elseif lon >= 180
            points[2, i + n + 1] = lon - 360
        else
            points[2, i + n + 1] = lon
        end
    end
    points
end
