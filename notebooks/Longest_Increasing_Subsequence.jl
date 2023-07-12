using Markdown
using Test
using BenchmarkTools
md"""
A subsequence of a permutation is a collection of elements of the permutation in the order that they appear. For example, (5, 3, 4) is a subsequence of (5, 1, 3, 4, 2).

A subsequence is increasing if the elements of the subsequence increase, and decreasing if the elements decrease. For example, given the permutation (8, 2, 1, 6, 5, 7, 4, 3, 9), an increasing subsequence is (2, 6, 7, 9), and a decreasing subsequence is (8, 6, 5, 4, 3). You may verify that these two subsequences are as long as possible.

Given: A positive integer n≤10000 followed by a permutation π of length n.

Return: A longest increasing subsequence of π, followed by a longest decreasing subsequence of π.
"""

function fill_dynamic_matrix(A)
  n = length(A)
  LIS, LDS = fill(1, n), fill(1, n)
  for i in eachindex(A)
    for j in 1:i-1
      if A[i] > A[j] && LIS[i] < LIS[j] + 1
        LIS[i] = LIS[j] + 1
      end
      if A[i] < A[j] && LDS[i] < LDS[j] + 1
        LDS[i] = LDS[j] + 1
      end
    end
  end
  LIS, LDS
end

function longest_sorted_subsequence(A)
  LIS, LDS = fill_dynamic_matrix(A)
  longest_increasing, longest_decreasing = max(LIS...), max(LDS...)
  increasing_subsequence, decreasing_subsequence = Int[], Int[]
  for i in reverse(1:length(A))
    if LIS[i] == longest_increasing
      push!(increasing_subsequence, A[i])
      longest_increasing -= 1
    end
    if LDS[i] == longest_decreasing
      push!(decreasing_subsequence, A[i])
      longest_decreasing -= 1
    end
  end
  reverse!(decreasing_subsequence)
  reverse!(increasing_subsequence)
  return (increasing_subsequence, decreasing_subsequence)
end

@test longest_sorted_subsequence([5, 1, 4, 2, 3]) == ([1, 2, 3], [5, 4, 3])

# Read the second line as array of integers of datasets/rosalind_lgis.txt
A = parse.(Int, split(readlines("datasets/rosalind_lgis.txt")[2], " "))
# Get the longest increasing and decreasing subsequence
increasing, decreasing = longest_sorted_subsequence(A)
# Print the longest increasing and decreasing subsequence separated by space into file
open("outputs/rosalind_lgis.txt", "w") do io
  write(io, join(increasing, " "))
  write(io, "\n")
  write(io, join(decreasing, " "))
end

@benchmark longest_sorted_subsequence(A)