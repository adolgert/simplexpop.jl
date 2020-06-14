"""
Does an insertion sort, returning the sorted array and a count of the
number of swaps.
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
