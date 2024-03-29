### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ b5830450-0a40-45f4-aed6-984cef4698ab
using BioSequences

# ╔═╡ 73de9f95-be32-42bd-a320-0b0b41daa92a
using FASTX

# ╔═╡ d8db9cd6-cb31-11ed-1939-8bc7c10fb486
md"""
Either strand of a DNA double helix can serve as the coding strand for RNA transcription. Hence, a given DNA string implies six total reading frames, or ways in which the same region of DNA can be translated into amino acids: three reading frames result from reading the string itself, whereas three more result from reading its reverse complement.

An open reading frame (ORF) is one which starts from the start codon and ends by stop codon, without any other stop codons in between. Thus, a candidate protein string is derived by translating an open reading frame into amino acids until a stop codon is reached.

Given: A DNA string s of length at most 1 kbp in FASTA format.

Return: Every distinct candidate protein string that can be translated from ORFs of s
. Strings can be returned in any order.
"""

# ╔═╡ 0b0ed89f-c4e2-4c7d-96eb-9fdb8e19bebd
open(FASTA.Reader, "../datasets/rosalind_orf.txt") do reader
    global sequence_dna = reader |> first |> sequence |> LongDNA{4}
end

# ╔═╡ 3e00dbfd-0d40-4a61-9af1-3ee373e954c4
sequence_dna_sample =
    dna"AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG"

# ╔═╡ 315e122b-c380-4a2a-8c12-70c571ce8d9d
sequence_rna = convert(LongRNA{4}, reverse_complement(sequence_dna))

# ╔═╡ 5a68f8be-c82f-437c-94a0-2ab20a037399
function find_codon(x, sequence)
    first.(findall(ExactSearchQuery(x), sequence))
end

# ╔═╡ 69890660-3c13-4dbb-acc0-c872851cbdf7
function find_all_frames(x)
    starts = find_codon(rna"AUG", x)
    terminations = collect(
        Iterators.flatten(find_codon(term, x) for term in [rna"UAA", rna"UAG", rna"UGA"]),
    )
    frames = []
    for s in starts
        for t in terminations
            if s <= t && (t - s) % 3 == 0
                push!(frames, s:t-1)
                break
            end
        end
    end
    frames
end

# ╔═╡ 4537523f-b60f-46b3-afb5-8fbdfaa43e07
directed = translate.([sequence_rna[f] for f in find_all_frames(sequence_rna)])

# ╔═╡ 98603f1a-f6d3-4fc3-8bb9-70e31de5c124
reversed =
    translate.([reverse_complement(sequence_rna)[f] for f in find_all_frames(sequence_rna)])

# ╔═╡ 611739f7-4a30-4aad-b85d-a28d00123cf9
sequences = unique(vcat(directed, reversed))

# ╔═╡ ef30f75b-79b5-4c1f-b14e-824dcee45374
#cat rosalind_orf.txt | cut -f1 -d"*" > rosalind_orf2.txt
open("../outputs/rosalind_orf.txt", "w") do io
    for sequence in sequences
        println(io, String(sequence))
    end
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BioSequences = "7e6ae17a-c86d-528c-b3b9-7f778a29fe59"
FASTX = "c2308a5c-f048-11e8-3e8a-31650f418d12"

[compat]
BioSequences = "~3.1.3"
FASTX = "~2.0.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0"
manifest_format = "2.0"

[[deps.Automa]]
deps = ["Printf", "ScanByte", "TranscodingStreams"]
git-tree-sha1 = "d50976f217489ce799e366d9561d56a98a30d7fe"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "0.8.2"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BioGenerics]]
deps = ["TranscodingStreams"]
git-tree-sha1 = "0b581906418b93231d391b5dd78831fdc2da0c82"
uuid = "47718e42-2ac5-11e9-14af-e5595289c2ea"
version = "0.1.2"

[[deps.BioSequences]]
deps = ["BioSymbols", "Random", "SnoopPrecompile", "Twiddle"]
git-tree-sha1 = "c96ede1c34ac948b108f11e4d9ae66df13d57454"
uuid = "7e6ae17a-c86d-528c-b3b9-7f778a29fe59"
version = "3.1.3"

[[deps.BioSymbols]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "2052c3ec7c41b69efa0e9ff7e2734aa6658d4c40"
uuid = "3c28c6f8-a34d-59c4-9654-267d177fcfa9"
version = "5.1.2"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.FASTX]]
deps = ["Automa", "BioGenerics", "BioSequences", "ScanByte", "SnoopPrecompile", "StringViews", "TranscodingStreams"]
git-tree-sha1 = "d9b3f4020a5afda65de0db1ee150409c1b69f44c"
uuid = "c2308a5c-f048-11e8-3e8a-31650f418d12"
version = "2.0.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.SIMD]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "8b20084a97b004588125caebf418d8cab9e393d1"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.4"

[[deps.ScanByte]]
deps = ["Libdl", "SIMD"]
git-tree-sha1 = "2436b15f376005e8790e318329560dcc67188e84"
uuid = "7b38b023-a4d7-4c5e-8d43-3f3097f304eb"
version = "0.3.3"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.StringViews]]
git-tree-sha1 = "6f58c457b8adeab3ac97cca6b3434ffefd7bd71e"
uuid = "354b36f9-a18e-4713-926e-db85100087ba"
version = "1.3.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "94f38103c984f89cf77c402f2a68dbd870f8165f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.11"

[[deps.Twiddle]]
git-tree-sha1 = "29509c4862bfb5da9e76eb6937125ab93986270a"
uuid = "7200193e-83a8-5a55-b20d-5d36d44a0795"
version = "1.1.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
"""

# ╔═╡ Cell order:
# ╠═d8db9cd6-cb31-11ed-1939-8bc7c10fb486
# ╠═b5830450-0a40-45f4-aed6-984cef4698ab
# ╠═73de9f95-be32-42bd-a320-0b0b41daa92a
# ╠═0b0ed89f-c4e2-4c7d-96eb-9fdb8e19bebd
# ╠═3e00dbfd-0d40-4a61-9af1-3ee373e954c4
# ╠═315e122b-c380-4a2a-8c12-70c571ce8d9d
# ╠═5a68f8be-c82f-437c-94a0-2ab20a037399
# ╠═69890660-3c13-4dbb-acc0-c872851cbdf7
# ╠═4537523f-b60f-46b3-afb5-8fbdfaa43e07
# ╠═98603f1a-f6d3-4fc3-8bb9-70e31de5c124
# ╠═611739f7-4a30-4aad-b85d-a28d00123cf9
# ╠═ef30f75b-79b5-4c1f-b14e-824dcee45374
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
