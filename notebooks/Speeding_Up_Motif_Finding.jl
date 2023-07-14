using Test
using Markdown
using BioSequences
using FASTX
using BenchmarkTools
md"""
A prefix of a length n string s is a substring s[1:j]; a suffix of s is a substring s[k:n].

The failure array of s is an array P of length n for which P[k] is the length of the longest substring s[j:k] that is equal to some prefix s[1:k−j+1], where j cannot equal 1 (otherwise, P[k] would always equal k). By convention, P[1]=0.

Given: A DNA string s (of length at most 100 kbp) in FASTA format.

Return: The failure array of s.
"""
function failure_array(x::LongSequence)
  n = length(x)
  f = zeros(Int, n)
  for k in 2:n
    j = f[k-1]
    while j > 0 && x[k] != x[j+1]
      j = f[j]
    end
    if x[k] == x[j+1]
      f[k] = j + 1
    end
  end
  return f
end

@test failure_array(dna"CAGCATGGTATCACAGCAGAG") == [0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 1, 2, 1, 2, 3, 4, 5, 3, 0, 0]

# read FASTA
records = [record for record in open(FASTA.Reader, "datasets/rosalind_kmp.txt")]
# use sequence function from FASTX to get the sequence
s = FASTX.sequence(LongDNA{2}, records[1])
@time f = failure_array(s)
# write to file
open("outputs/rosalind_kmp.txt", "w") do io
  write(io, join(f, " "))
end

@benchmark failure_array(s)