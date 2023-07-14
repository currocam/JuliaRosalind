using BioSequences
using BioAlignments
using FASTX
using BenchmarkTools
using Test
using Markdown
md"""
Problem

For two strings s1 and s2 of equal length, the p-distance between them, denoted dp(s1,s2), is the proportion of corresponding symbols that differ between s1 and s2.

For a general distance function d on n taxa s1,s2,…,sn (taxa are often represented by genetic strings), we may encode the distances between pairs of taxa via a distance matrix D in which Di,j=d(si,sj).

Given: A collection of n (n≤10) DNA strings s1,…,sn of equal length (at most 1 kbp). Strings are given in FASTA format.

Return: The matrix D corresponding to the p-distance dp on the given strings. As always, note that your answer is allowed an absolute error of 0.001.
"""
s1, s2 = dna"TTTCCATTTA", dna"GATTCATTTC"

function distance(s1::LongSequence, s2::LongSequence)
  return sum(s1 .!= s2) / length(s1)
end

function distance_matrix(sequences)
  n = length(sequences)
  D = zeros(Float64, n, n)
  # Fill upper matrix
  for i in 1:n
    for j in i+1:n
      if i != j
        D[i, j] = D[j, i] = distance(sequences[i], sequences[j])
      end
    end
  end
  D
end

test_sequences = [dna"TTTCCATTTA", dna"GATTCATTTC", dna"TTTCCATTTT", dna"GTTCCATTTA"]
@test distance_matrix(sequences) ≈ [0.0 0.4 0.1 0.1; 0.4 0.0 0.4 0.3; 0.1 0.4 0.0 0.2; 0.1 0.3 0.2 0.0]

# Read fasta file
sequences = [LongDNA{2}(FASTX.sequence(record)) for record in open(FASTA.Reader, "datasets/rosalind_pdst.txt")]

@time A = distance_matrix(sequences)

# Write to file
open("outputs/rosalind_pdst.txt", "w") do io
  for row in eachrow(A)
    write(io, join(row, " "), "\n")
  end
end

@benchmark distance_matrix(sequences)