using Markdown
using Test
using BioSequences
using Combinatorics
using FASTX
md"""
For a collection of strings, a larger string containing every one of the smaller strings as a substring is called a superstring.

By the assumption of parsimony, a shortest possible superstring over a collection of reads serves as a candidate chromosome.

Given: At most 50 DNA strings of approximately equal length, not exceeding 1 kbp, in FASTA format (which represent reads deriving from the same strand of a single linear chromosome).

The dataset is guaranteed to satisfy the following condition: there exists a unique way to reconstruct the entire chromosome from these reads by gluing together pairs of reads that overlap by more than half their length.

Return: A shortest superstring containing all the given strings (thus corresponding to a reconstructed chromosome).
"""

# Define overlap function, which returns the length of the maximum overlap between the end of left and the beginning of right
function overlap(left, right)
    for i in eachindex(left)
        if left[i:end] == right[1:min(length(left) - i + 1, end)]
            return length(left) - i + 1
        end
    end
    return 0
end

function overlap_more_than_half(left, right)
    overlap(left, right) > length(left) ÷ 2
    # overlap = length(left) ÷ 2 - 1
    # @inbounds left[overlap:end] == right[1:end-overlap+1]
end


function overlap_graph(sequences)
    graph = Dict{Int,Int}()
    # Iterate over all pairs of strings
    for (i, j) in Iterators.product(eachindex(sequences), eachindex(sequences))
        if i == j
            continue
        end
        if overlap_more_than_half(sequences[i], sequences[j])
            graph[i] = j
        end
    end
    graph
end

function merge(left, right)
    for i in eachindex(left)
        length_left = length(left)
        @inbounds if left[i:end] == right[1:min(length_left - i + 1, end)]
            return left[1:i-1] * right
        end
    end
end

function shortest_superstring(sequences)
    g = overlap_graph(sequences)
    # left sequence is the setdiff of keys and values
    left = setdiff(keys(g), values(g)) |> first
    contig = sequences[left]
    while haskey(g, left)
        right = g[left]
        contig = merge(contig, sequences[right])
        left = right
    end
    contig
end

shortest_superstring([dna"ATTAGACCTG", dna"CCTGCCGGAA", dna"AGACCTGCCG", dna"GCCGGAATAC"])


@testset "Shortest Superstring" begin
    for permutation in
        permutations([dna"ATTAGACCTG", dna"CCTGCCGGAA", dna"AGACCTGCCG", dna"GCCGGAATAC"])
        @test shortest_superstring(permutation) == dna"ATTAGACCTGCCGGAATAC"
    end

    for permutation in permutations([dna"ATTCG", dna"TCGACGA", dna"ACGATGAGA"])
        @test shortest_superstring(permutation) == dna"ATTCGACGATGAGA"
    end
end

infile = "datasets/rosalind_long.txt"
## Open fasta file
@time sequences = [LongDNA{2}(sequence(record)) for record in open(FASTA.Reader, infile)]
@time contig = shortest_superstring(sequences)
open("outputs/rosalind_long.txt", "w") do io
    write(io, String(contig))
end
