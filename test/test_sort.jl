struct SortTrial
    a
    cnt
end

trials = [
    SortTrial([1,2,3], 0),
    SortTrial([2, 1, 3], 1),
    SortTrial([2, 1], 1),
    SortTrial([2, 1, 4, 3], 2),
    SortTrial([4, 1, 2, 3], 3),
]
for t in trials
    cnt = countedinsertionsort!(t.a)
    @test t.a == sort(t.a)
    @test cnt == t.cnt
end
