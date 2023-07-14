using BioAlignments, BioSequences, FASTX, Markdown, Test, BenchmarkTools
md"""
Given: Two DNA strings s and t (each having length at most 1 kbp) in FASTA format.
Return: A longest common subsequence of s
"""

function longest_common_subsequence(s1::LongDNA, s2::LongDNA)
  scoremodel = AffineGapScoreModel(DichotomousSubstitutionMatrix(1, -1), gap_open=0, gap_extend=0)
  aln = pairalign(GlobalAlignment(), s1, s2, scoremodel) |> alignment
  # Remove gaps from alignment
  LongDNA{4}([a for (a, b) in aln if a != DNA_Gap && b != DNA_Gap])
end

s1, s2 = dna"AACCTTGG", dna"ACACTGTGA"
@test longest_common_subsequence(s1, s2) |> length == 6

# Read fasta file
sequences = [
  LongDNA{4}(FASTX.sequence(record)) for
  record in open(FASTA.Reader, "datasets/rosalind_lcsq.txt")
]

@time lcss = longest_common_subsequence(sequences[1], sequences[2])

# Write to file
open("outputs/rosalind_lcsq.txt", "w") do io
  write(io, String(lcss))
end

@benchmark longest_common_subsequence(sequences[1], sequences[2])
