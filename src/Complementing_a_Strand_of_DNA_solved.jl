using Markdown, Test, BenchmarkTools
md"""
Given: A DNA string s of length at most 1000 bp.
Return: The reverse complement sc of s.
"""
function complementate(x::String)
    dict = Dict('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A')
    map(y -> dict[y], x) |> reverse
end

@test complementate("AAAACCCGGT") == "ACCGGGTTTT"

f = open("datasets/rosalind_revc.txt")
open("outputs/rosalind_revc.txt", "w") do file
    sequence = readlines(f) |> first
    write(file, complementate(sequence))
end

@benchmark complementate(sequence)