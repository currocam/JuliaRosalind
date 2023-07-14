using Markdown
using BioAlignments
using BioSequences
using Test
using FASTX
using BenchmarkTools
md"""
A subsequence of a string is a collection of symbols contained in order (though not necessarily contiguously) in the string (e.g., ACG is a subsequence of TATGCTAAGATC). The indices of a subsequence are the positions in the string at which the symbols of the subsequence appear; thus, the indices of ACG in TATGCTAAGATC can be represented by (2, 5, 9).

As a substring can have multiple locations, a subsequence can have multiple collections of indices, and the same index can be reused in more than one appearance of the subsequence; for example, ACG is a subsequence of AACCGGTT in 8 different ways.

Given: Two DNA strings s
and t

(each of length at most 1 kbp) in FASTA format.

Return: One collection of indices of s
in which the symbols of t appear as a subsequence of s. If multiple solutions exist, you may return any one.
"""

function spliced_motif(s, t)
    scoremodel = AffineGapScoreModel(
        DichotomousSubstitutionMatrix(1, -1),
        gap_open = 0,
        gap_extend = 0,
    )
    aln = pairalign(GlobalAlignment(), s, t, scoremodel) |> alignment
    filter(i -> seq2aln(aln, i)[2] == OP_SEQ_MATCH, eachindex(s))
end

@test spliced_motif(dna"ACGTACGTGACG", dna"GTA") == [3, 4, 5]

# Read fasta sequence one and Two
records = [record for record in open(FASTA.Reader, "datasets/rosalind_sseq.txt")]
# use sequence function from FASTX to get the sequence
s = FASTX.sequence(LongDNA{2}, records[1])
t = FASTX.sequence(LongDNA{2}, records[2])
@benchmark positions = spliced_motif(s, t)

# Print the positions separated by space into file
open("outputs/rosalind_sseq.txt", "w") do io
    write(io, join(positions, " "))
end
