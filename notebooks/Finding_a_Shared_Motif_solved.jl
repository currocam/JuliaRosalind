using Markdown
md"""
Problem
A common substring of a collection of strings is a substring of every member of the collection. We say that a common substring is a longest common substring if there does not exist a longer common substring. For example, "CG" is a common substring of "ACGTACGT" and "AACCGTATA", but it is not as long as possible; in this case, "CGTA" is a longest common substring of "ACGTACGT" and "AACCGTATA".

Note that the longest common substring is not necessarily unique; for a simple example, "AA" and "CC" are both longest common substrings of "AACC" and "CCAA".

Given: A collection of k
 (k≤100
) DNA strings of length at most 1 kbp each in FASTA format.

Return: A longest common substring of the collection. (If multiple solutions exist, you may return any single solution.)
"""
using BioSequences
using Test
using FASTX
using BenchmarkTools
using Distributed


sequences = [dna"ATGTTTT", dna"ATCGTTT", dna"ATGTTCT"]

function generate_substrings(s::LongSequence)
    n = length(s)
    f(length) = (LongDNA{4}(s[i:i+length-1]) for i = 1:n-length+1)
    return Iterators.flatten((f(length) for length = n:-1:1))
end


## Check if a substring is a substring of all sequences
function find_substring(sequences, substring)
    query = ExactSearchQuery(substring)
    for sequence in sequences
        if !occursin(query, sequence)
            return false
        end
    end
    return true
end

function LCSM(sequences)
    ## Find shortst sequence
    shortest = length.(sequences) |> argmin
    sequences_except_shortest = sequences[setdiff(1:end, shortest)]
    for substring in generate_substrings(sequences[shortest])
        ## All sequences except shortest
        if find_substring(sequences_except_shortest, substring)
            return substring
        end
    end
    return nothing
end


@test [dna"GATTACA", dna"TAGACCA", dna"ATACA"] |> LCSM |> length == 2

infile = "datasets/rosalind_lcsm.txt"
## Open fasta file
sequences = [LongDNA{2}(sequence(record)) for record in open(FASTA.Reader, infile)]

longest_substring = LCSM(sequences)
# Print into file
open("outputs/rosalind_lcsm.txt", "w") do io
    write(io, String(longest_substring))
end
