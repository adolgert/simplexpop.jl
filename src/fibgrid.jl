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
function fibonacci_grid(n)
    invϕ = 2 / (1 + sqrt(5))
    points = zeros(Float64, 2, 2n + 1)
    for i in -n:n
        points[1, i + n + 1] = asin(2i / (2n + 1))
        points[2, i + n + 1] = 2π * i/invϕ
    end
    points
end


"""
Hexagonal grid built from a Bravais lattice.
Two vectors: r1, r2, of the same length.
The angle between them is 120 degrees.

a1 = (r, 0)
a2 = (-r/sqrt(3), 2r/sqrt(3))
x = n1 a1 + n2 a2
"""
