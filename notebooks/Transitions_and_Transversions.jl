using BioSequences
using Test
using FASTX
using Markdown
using BenchmarkTools

md"""
For DNA strings s1 and s2 having the same length, their transition/transversion ratio R(s1,s2)

is the ratio of the total number of transitions to the total number of transversions, where symbol substitutions are inferred from mismatched corresponding symbols as when calculating Hamming distance (see “Counting Point Mutations”).

Given: Two DNA strings s1
and s2

of equal length (at most 1 kbp).

Return: The transition/transversion ratio R(s1,s2)
.
"""

PURINES = [DNA_A, DNA_G]
PYRIMIDINES = [DNA_C, DNA_T]
is_transition(a, b) =
    (a in PURINES && b in PURINES) || (a in PYRIMIDINES && b in PYRIMIDINES)
function transition_transversions_ratio(s1, s2)
    transitions, transversions = 0, 0
    for (i, j) in zip(s1, s2)
        if i == j
            continue
        end
        is_transition(i, j) ? transitions += 1 : transversions += 1
    end
    transitions / transversions
end

@test transition_transversions_ratio(
    dna"GCAACGCACAACGAAAACCCTTAGGGACTGGATTATTTCGTGATCGTTGTAGTTATTGGAAGTACGGGCATCAACCCAGTT",
    dna"TTATCTGACAAAGAAAGCCGTCAACGGCTGGATAATTTCGCGATCGTGCTGGTTACTGGCGGTACGAGTGTTCCTTTGGGT",
) ≈ 1.21428571429

# Read fasta sequence one and Two
records = [record for record in open(FASTA.Reader, "datasets/rosalind_tran.txt")]
# use sequence function from FASTX to get the sequence
s1 = FASTX.sequence(LongDNA{2}, records[1])
s2 = FASTX.sequence(LongDNA{2}, records[2])
ratio = transition_transversions_ratio(s1, s2)
@benchmark transition_transversions_ratio(s1, s2)
