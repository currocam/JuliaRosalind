using BioSequences
using Markdown
using Test
using FASTX
using BenchmarkTools
md"""
For a fixed positive integer k, order all possible k-mers taken from an underlying alphabet lexicographically.

Then the k-mer composition of a string s can be represented by an array A for which A[m] denotes the number of times that the mth k-mer (with respect to the lexicographic order) appears in s.

Given: A DNA string s in FASTA format (having length at most 100 kbp).

Return: The 4-mer composition of s.
"""
# Get the integer representation of the sequence
nucleotide2int = Dict{DNA,Int8}(DNA_A => 0, DNA_C => 1, DNA_G => 2, DNA_T => 3)
function lexicographical_position(s::LongDNA, k::Int)
  pos = 0
  for i in 1:k
    pos += nucleotide2int[s[i]] * 4^(k - i)
  end
  pos + 1
end

# create test set for lexicographical_position
@test lexicographical_position(dna"AAAA", 4) == 1
@test lexicographical_position(dna"AAAC", 4) == 2
@test lexicographical_position(dna"AAAG", 4) == 3
@test lexicographical_position(dna"TTTT", 4) == 256

function count_kmers(s::LongDNA, k::Int)
  counts = zeros(Int, 4^k)
  for i in 1:length(s)-k+1
    counts[lexicographical_position(s[i:i+k-1], k)] += 1
  end
  counts
end

@test count_kmers(dna"", 1) == [0, 0, 0, 0]
@test count_kmers(dna"ACG", 1) == [1, 1, 1, 0]

s = open(FASTA.Reader, "datasets/rosalind_kmer.txt") |>
    first |> FASTX.sequence |> LongDNA{2}

@time counts = count_kmers(s, 4)

open("outputs/rosalind_kmer.txt", "w") do io
  for i in counts
    write(io, "$i ")
  end
end

@benchmark count_kmers(s, 4)
