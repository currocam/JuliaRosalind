using Markdown
using Test
using BenchmarkTools
md"""
Assume that an alphabet A has a predetermined order; that is, we write the alphabet as a permutation A=(a1,a2,…,ak), where a1<a2<⋯<ak. For instance, the English alphabet is organized as (A,B,…,Z).

Given two strings s
and t having the same length n, we say that s precedes t in the lexicographic order (and write s<Lext) if the first symbol s[j] that doesn't match t[j] satisfies sj<tj in A.

Given: A collection of at most 10 symbols defining an ordered alphabet, and a positive integer n
(n≤10).

Return: All strings of length n
that can be formed from the alphabet, ordered lexicographically (use the standard order of symbols in the English alphabet).
"""

lexicographic_permutations(alphabet, n) = Iterators.product(Iterators.repeated(alphabet, n)...) |> collect |> vec |> sort

@test all(lexicographic_permutations(split("A C G T"), 2) .== lexicographic_permutations(split("A G C T"), 2))

# Print the permutations separated by new line into file using println
open("outputs/rosalind_lexf.txt", "w") do io
  @time permutations = lexicographic_permutations(split("A B C D E F"), 3)
  for permutation in permutations
    write(io, join(permutation, ""))
    write(io, "\n")
  end
end

@benchmark permutations = lexicographic_permutations(split("A B C D E F"), 5)

