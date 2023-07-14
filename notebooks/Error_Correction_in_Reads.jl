using BioAlignments
using BioSequences
using FASTX
using Markdown
using Test
using BenchmarkTools
md"""
Given: A collection of up to 1000 reads of equal length (at most 50 bp) in FASTA format. Some of these reads were generated with a single-nucleotide error. For each read s in the dGiven: A collection of up to 1000 reads of equal length (at most 50 bp) in FASTA format. Some of these reads were generated with a single-nucleotide error. For each read s in the dataset, one of the following applies: 

- s was correctly sequenced and appears in the dataset at least twice (possibly as a reverse complement);

- s is incorrect, it appears in the dataset exactly once, and its Hamming distance is 1 with respect to exactly one correct read in the dataset (or its reverse complement).

Return: A list of all corrections in the form "[old read]->[new read]". (Each correction must be a single symbol substitution, and you may return the corrections in any order.)ataset, one of the following applies: 
"""

function distance(s1::LongDNA, s2::LongDNA)
    min(sum(s1 .!= s2), sum(reverse_complement(s1) .!= s2))
end
function distance_matrix(sequences)
    n = length(sequences)
    D = zeros(Int64, n, n)
    for i = 1:n
        @inbounds for j = i+1:n
            D[i, j] = D[j, i] = distance(sequences[i], sequences[j])
        end
    end
    return D
end



test_sequences = [
    dna"TCATC",
    dna"TTCAT",
    dna"TCATC",
    dna"TGAAA",
    dna"GAGGA",
    dna"TTTCA",
    dna"ATCAA",
    dna"TTGAT",
    dna"TTTCC",
]

function corrected_reads(sequences)
    D = distance_matrix(sequences)
    duplicates = Set(x.I[1] for x in findall(D .== 0) if x[1] > x[2])
    corrections = Dict{LongDNA,LongDNA}()
    for snp in findall(D .== 1)
        if snp.I[1] in duplicates
            seq1 = sequences[snp.I[1]]
            seq2 = sequences[snp.I[2]]
            if count(x -> x[1] != x[2], zip(seq1, seq2)) == 1
                corrections[seq2] = seq1
            else
                corrections[seq2] = reverse_complement(seq1)
            end
        end
    end
    return corrections
end


expected =
    Dict(dna"TTCAT" => dna"TTGAT", dna"GAGGA" => dna"GATGA", dna"TTTCC" => dna"TTTCA")
@test expected == corrected_reads(test_sequences)

sequences = [
    LongDNA{2}(FASTX.sequence(record)) for
    record in open(FASTA.Reader, "datasets/rosalind_corr.txt")
]
@time corrections = corrected_reads(sequences)
open("outputs/rosalind_corr.txt", "w") do io
    for old in keys(corrections)
        write(io, "$(old)->$(corrections[old])\n")
    end
end

@benchmark corrected_reads(sequences)
