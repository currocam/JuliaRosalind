using Markdown, Test, BenchmarkTools
md"""
Given: A DNA string t having length at most 1000 nt.
Return: The transcribed RNA string of t.
"""

function transcribe(dna::String)
    replace(dna, "T" => "U")
end
@test transcribe("GATGGAACTTGACTACGTAAATT") == "GAUGGAACUUGACUACGUAAAUU"
begin
    f = open("datasets/rosalind_rna.txt")
    lines = readlines(f)
    sequence = lines[1]
end

open("outputs/rosalind_rna.txt", "w") do file
    write(file, transcribe(sequence))
end
@benchmark transcribe(sequence)