
function xorshift64(a::UInt64)::UInt64
	a ⊻= (a << 13)
	a ⊻= (a >> 7)
	a ⊻= (a << 17)
	a
end


function unique_sequence(n, seed)
	a = UInt64(seed)
	m = zeros(UInt64, n)
	for i in 1:n
		a = xorshift64(a)
		m[i] = a
	end
	m
end


function unique_sequence(::Type{T}, n, seed) where {T <: Int64}
	m = unique_sequence(n, seed)
	for i in 1:n
		m[i] >>= 1
	end
	Int64.(m)
end

function unique_sequence(::Type{T}, n, seed) where {T <: Float64}
	m = unique_sequence(n, seed)
	Float64.(m)
end
