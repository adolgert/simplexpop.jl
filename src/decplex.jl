"""
Simplex in 2D.
PyDEC: Software and Algorithms for Discretization of Exterior Calculus
by Bell and Hirani.
"""
struct Simplex2D
    # vertex locations
    v::Array{Float64,2}
    # vertices that form triangles
    s::Array{Int64,2}
end


function buildsimplex2d(v, s)
    # The simplex has to be normalized.
    # The hilbert order increases locality of the geometry.
    reorder = hilbert_order(v)
    v = v[:, reorder]
    remap = sortperm(reorder)
    for j in 1:size(s, 2)
        for i in 1:size(s, 1)
            s[i, j] = remap[s[i, j]]
        end
    end
    # We will rely on simplex vertices being sorted.
    # circular permuations perserve orientation.
    for shift_idx in 1:size(s, 2)
        verts = s[:, shiftidx]
        s[:, shiftidx] = circshift(verts, findmin(verts))
    end
    s = sort(s, dims = 1)
    Simplex2D(v, s)
end


"""
Assume s is a list of triangles.
Pydec paper page 13. But transposed b/c this is column first ordering.
"""
function boundary_faces(s)
    d, n = size(s)  # dimension (3) and number of simplices
    splus = zeros(Int64, d + 2, n)
    splus[1:d, :] = s
    # need orientation of simplex for edges.
    for orient_idx in 1:d
        orient = iseven(countedinsertionsort!(splus[1:d, orient_idx]))
        splus[d + 1, orient_idx] = orient ? 1 : -1
        splus[d + 2, orient_idx] = orient_idx
    end
    resultant = zeros(Int64, d + 1, d * n)
    for vert in 1:d
        for tri in 1:n
            for copy_idx in 1:vert
                cnt = 1
                if copy_idx != vert
                    resultant[cnt, tri] = splus[copy_idx, tri]
                    cnt += 1
                end
                resultant[d, tri] = splus[d + 1, tri] * (-1)^(vert - 1)
                resultant[d + 1, tri] = splus[d + 2, tri]
            end
        end
    end
    ordered_resultant = sort(resultant, dims = 1)

    faces_pad = zeros(Int64, d - 1, size(resultant, 2))
    face_cnt = 1
    faces_pad[:, 1] = resultant[1:(d - 1), 1]
    for fidx in 2:size(resultant, 2)
        if faces_pad[1:(d - 1), fidx] != faces_pad[:, face_cnt]
            face_cnt += 1
            faces_pad[:, face_cnt] = resultant[1:(d - 1), fidx]
        end
    end
    faces = faces_pad[:, 1:face_cnt]

    edges = SparseMatrixCSC{Int64,Int64}(
        n,  # row count is number of simplices
        size(resultant, 2),  # column count is number of edges
        1:size(resultant, 2),  # one column per edge
        resultant[d + 1, :],  # row indices are the simplex index
        resultant[d, :]  # values are the orientation of the edge
        )
    edges, faces
end
