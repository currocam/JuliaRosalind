using Markdown, Test, BenchmarkTools
import StatsBase: countmap
md"""
Given: A DNA string s of length at most 1000 nt.
Return: Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in s.
"""

function count_dna(x::String)
    y = countmap(x)
    [y['A'], y['C'], y['G'], y['T']]
end

begin
    x = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
    @test join([string(x) for x in count_dna(x)], " ") == "20 12 17 21"
end

begin
    f = open("datasets/rosalind_dna.txt")
    lines = readlines(f)
    sequence = lines[1]
end

join([string(x) for x in count_dna(sequence)], " ")

@benchmark count_dna(sequence)