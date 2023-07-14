using Markdown, Test, BenchmarkTools
md"""
Given: Three positive integers k, m, and n, representing a population containing k+m+n organisms: k individuals are homozygous dominant for a factor, m are heterozygous, and n are homozygous recessive.
Return: The probability that two randomly selected mating organisms will produce an individual possessing a dominant allele (and thus displaying the dominant phenotype). Assume that any two organisms can mate.
"""

function first_law(k::Int, m::Int, n::Int)
  pop = k + m + n
  (4 * (k * (k - 1) + 2 * k * m + 2 * k * n + m * n) + 3 * m * (m - 1)) / (4 * pop * (pop - 1))
end

@test first_law(2, 2, 2) ≈ 0.78333333333