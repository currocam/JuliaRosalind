using BioSequences
using Markdown
using Test
using BenchmarkTools
md"""
Given: A positive integer N≤100000, a number x between 0 and 1, and a DNA string s of length at most 10 bp.

Return: The probability that if N random DNA strings having the same length as s are constructed with GC-content x (see “Introduction to Random Strings”), then at least one of the strings equals s. We allow for the same random string to be created more than once.
"""


weights = x -> Dict(DNA_A => (1 - x) / 2, DNA_C => x / 2, DNA_G => x / 2, DNA_T => (1 - x) / 2) # A C G T
function probability_motif(N, w, s)
  log_p_sucess = BigFloat(0)
  for char in s
    log_p_sucess += log(w[char])
  end
  1 - exp(log(1 - exp(log_p_sucess)) * N)
end

@test isapprox(probability_motif(90000, weights(0.6), dna"ATAGCCGA"), 0.689, atol=0.001)

N = 88418
cg_content = 0.506549
s = dna"TAACCGCCTG"
p = probability_motif(N, weights(cg_content), s)

@benchmark probability_motif(N, weights(cg_content), s)