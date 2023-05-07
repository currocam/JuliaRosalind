
using Markdown
using InteractiveUtils

using Test

using BioSequences

using FASTX

md"""
After identifying the exons and introns of an RNA string, we only need to delete the introns and concatenate the exons to form a new string ready for translation.

Given: A DNA string s
(of length at most 1 kbp) and a collection of substrings of s

acting as introns. All strings are given in FASTA format.

Return: A protein string resulting from transcribing and translating the exons of s
. (Note: Only one solution will exist for the dataset provided.)
"""

function translate_with_introns(sequence, introns)
	queries = ExactSearchQuery.(introns)
	indexes = [findall(query, sequence) for query in queries] |>
		Iterators.flatten |> collect
	is_intron = x -> any([x in intron for intron in indexes])
	exons = sequence[[i for i in 1:length(sequence) if !is_intron(i)]]
	protein = translate(exons)
	termination = findfirst(AA_Term, protein)
	protein[1:termination-1]
	end

sequence_sample = dna"ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG"

introns = [dna"ATCGGTCGAA", dna"ATCGGTCGAGCGTGT"]

@test translate_with_introns(sequence_sample, introns) == aa"MVYIADKQHVASREAYGHMFKVCA"

open("datasets/rosalind_splc.txt") do io
    global sequences = [LongDNA{4}(sequence(r)) for r in FASTAReader(io)]
end;

result = translate_with_introns(sequences[1], sequences[2:end])
println(result)