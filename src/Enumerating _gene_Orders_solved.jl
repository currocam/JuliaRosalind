using Combinatorics
using Test

solve = x -> collect(permutations(1:x))

@test length(solve(3)) == 6

file = open("outputs/rosalind_perm.txt", "w")
result = solve(7)
println(file, length(result))
for row in result
    println(file, join(row, " "))
end
close(file)
