using Markdown
using Test
using BenchmarkTools
using Combinatorics
using IterTools
md"""
A signed permutation of length n is some ordering of the positive integers {1,2,…,n} in which each integer is then provided with either a positive or negative sign (for the sake of simplicity, we omit the positive sign). For example, π=(5,−3,−2,1,4) is a signed permutation of length 5.

Given: A positive integer n≤6.

Return: The total number of signed permutations of length n, followed by a list of all such permutations (you may list the signed permutations in any order).
"""

function signed_permutations(n)
  is_valid = x -> (abs.(x) |> sort) == (1:n)
  all_permutations = vcat(1:n, -1:-1:-n) |> x -> permutations(x, n)
  Iterators.filter(is_valid, all_permutations)
end
@test signed_permutations(2) |> collect |> length == 8

@time open("outputs/rosalind_sign.txt", "w") do io
  n = 4
  x = signed_permutations(n)
  println(io, 2^n * factorial(n))
  for permutation in x
    println(io, join(permutation, " "))
  end
end
