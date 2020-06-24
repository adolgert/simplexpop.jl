"""
    countedinsertionsort!(a)

Sorts an array and returns the number of swaps required to sort.

This is useful when you want to know the number
of transpositions of vertices in a simplex. This determines the new
orientation of the simplex relative to the old one.
"""
function countedinsertionsort!(a)
    swaps = zero(Int64)
    for i in 2:length(a)
        key = a[i]
        j = i - 1
        while j >= 1 && key < a[j]
            a[j + 1] = a[j]
            j -= 1
            a[j + 1] = key
            swaps += 1
        end
    end
    swaps
end
