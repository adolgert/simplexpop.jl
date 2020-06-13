table1 = [
  1 2
  0 3
]
table1 = table1[end:-1:1, :]
table2 = [
  15 12 11 10
  14 13  8  9
   1  2  7  6
   0  3  4  5
]
table2 = table2[end:-1:1, :]

## Test encode
zz = Set(Int64[])
for x in 0:15
  for y in 0:15
    z = encode_hilbert_zero(x, y)
    @show (x, y, z)
    @test z ∉ zz
    push!(zz, z)
  end
end

## Test decode
table2ind = CartesianIndices(table2)
for z in 0:15
  x, y = decode_hilbert_zero(z)
  found = table2ind[table2 .== z][1]
  @show z
  @show found
  @show (x, y)
  @test found[2] - 1 == x
  @test found[1] - 1 == y
end

xy = Set(Tuple{Int64, Int64}[])
for z in 0:255
  x, y = decode_hilbert_zero(z)
  @test (x, y) ∉ xy
  push!(xy, (x, y))
end
for i in 0:15
  for j in 0:15
    @test (i, j) in xy
  end
end


## Test against each other.

for z in 1:400
  x, y = decode_hilbert_zero(z)
  zz = encode_hilbert_zero(x, y)
  @test zz == z
end
