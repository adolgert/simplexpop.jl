# Columns are odd resolution, even resolution (rmin)
# Columns are quadrants in x and y like you'd graph it.
# 1 2
# 0 3
hilbert_variant_encode_table = Function[
    ((x, y, w) -> (x, y))          ((x, y, w) -> (x, y))
    ((x, y, w) -> (y-w, x))        ((x, y, w) -> (y, x-w))
    ((x, y, w) -> (y-w, x-w))      ((x, y, w) -> (y-w, x-w))
    ((x, y, w) -> (2w-x-1, w-y-1)) ((x, y, w) -> (w-x-1, 2w-y-1))
]

"""
Hilbert curve variant.
A new algorithm for encoding
and decoding the Hilbert order
Ningtao Chen∗,†, Nengchao Wang and Baochang Shi
Softw. Pract. Exper. 2007; 37:897–908

0 <= x < 2^n
0 <= y < 2^n
0 <= z < 4^n
"""
function encode_hilbert_zero(x, y)
    z = zero(x)
    if x == 0 && y == 0
        return z
    end
    rmin = convert(typeof(x), floor(log2(max(x, y))) + 1)
    w = 1<<(rmin - 1)
    while rmin > 0
        if rmin&1 == 1  # odd
            quadrant = x < w ? (y < w ? 0 : 1) : (y < w ? 3 : 2)
            parity = 1
        else  # even
            quadrant = x < w ? (y < w ? 0 : 3) : (y < w ? 1 : 2)
            parity = 2
        end
        z = (z<<2) + quadrant
        x, y = hilbert_variant_encode_table[quadrant+1, parity](x, y, w)
        rmin -= 1
        w >>= 1
    end
    z
end


hilbert_variant_decode_table = Function[
    ((x, y, w) -> (x, y))          ((x, y, w) -> (x, y))
    ((x, y, w) -> (y, x+w))        ((x, y, w) -> (y+w, x))
    ((x, y, w) -> (y+w, x+w))      ((x, y, w) -> (y+w, x+w))
    ((x, y, w) -> (2w-x-1, w-y-1)) ((x, y, w) -> (w-x-1, 2w-y-1))
]
function decode_hilbert_zero(z)
    r = z & 3
    x, y = typeof(z).([(0, 0), (0, 1), (1, 1), (1, 0)][r + 1])
    z >>= 2
    rmin = 2
    w = convert(typeof(z), 2)
    while z > 0
        r = z & 3
        parity = 2 - rmin&1
        x, y = hilbert_variant_decode_table[r+1, parity](x, y, w)
        z >>= 2
        rmin += 1
        w <<= 1
    end
    x, y
end


"""
A 1-based Hilbert curve.
"""
encode_hilbert(x, y) = encode_hilbert_zero(x - 1, y - 1) + 1
function decode_hilbert(z)
    x, y = decode_hilbert_zero(z - 1)
    x + 1, y + 1
end
