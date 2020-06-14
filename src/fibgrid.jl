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
