### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ a1fa8018-d621-4751-8a7c-442108feced0
using BioSequences

# ╔═╡ 1ff9ac81-4610-4e4a-8375-a3d023475353
using FASTX

# ╔═╡ c505364d-ba7a-4a2b-aaba-3edd5947e667
using Combinatorics

# ╔═╡ c023e2e1-23af-42a8-8e9e-bac435639b27
using Test

# ╔═╡ afe74614-ca56-11ed-0c86-45e6a934728e
md"""
A graph whose nodes have all been labeled can be represented by an adjacency list, in which each row of the list contains the two node labels corresponding to a unique edge.

A directed graph (or digraph) is a graph containing directed edges, each of which has an orientation. That is, a directed edge is represented by an arrow instead of a line segment; the starting and ending nodes of an edge form its tail and head, respectively. The directed edge with tail v
and head w is represented by (v,w) (but not by (w,v)). A directed loop is a directed edge of the form (v,v)

.

For a collection of strings and a positive integer k
, the overlap graph for the strings is a directed graph Ok in which each string is represented by a node, and string s is connected to string t with a directed edge when there is a length k suffix of s that matches a length k prefix of t, as long as s≠t; we demand s≠t

to prevent directed loops in the overlap graph (although directed cycles may be present).

Given: A collection of DNA strings in FASTA format having total length at most 10 kbp.

Return: The adjacency list corresponding to O3
. You may return edges in any order.
"""

# ╔═╡ f999b4e1-567e-463e-8db4-896c5a15c208
function overlap_graphs(records, k)
    graph = Array{Tuple{String,String},1}()
    for seq1 in records, seq2 in records
        if seq1 == seq2
            continue
        end
        if sequence(seq1)[end-k+1:end] == sequence(seq2)[1:k]
            push!(graph, (identifier(seq1), identifier(seq2)))
        end
    end
    graph
end

# ╔═╡ 44da0576-5626-4541-afaa-624e051cb8fe
begin
    headers = [
        "Rosalind_0498",
        "Rosalind_2391",
        "Rosalind_2323",
        "Rosalind_0442",
        "Rosalind_5013",
    ]
    sequences = [dna"AAATAAA", dna"AAATTTT", dna"TTTTCCC", dna"AAATCCC", dna"GGGTGGG"]
    records = [FASTA.Record(h, s) for (h, s) in zip(headers, sequences)]
    expected = [
        ("Rosalind_0498", "Rosalind_2391"),
        ("Rosalind_0498", "Rosalind_0442"),
        ("Rosalind_2391", "Rosalind_2323"),
    ]
    @test overlap_graphs(records, 3) == expected
end

# ╔═╡ 9edde78d-b1b0-4259-922d-0080bd5c0444
begin
    open(FASTA.Reader, "../datasets/rosalind_grph.txt") do reader
        global graph = overlap_graphs([x for x in reader], 3)
    end
end

# ╔═╡ 705aba37-01cf-4d63-9f1b-479609552959
begin
    outfile = open("../outputs/rosalind_grph.txt", "w")
    for (from, dest) in graph
        println(outfile, join([from, dest], " "))
    end
    close(outfile)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BioSequences = "7e6ae17a-c86d-528c-b3b9-7f778a29fe59"
Combinatorics = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
FASTX = "c2308a5c-f048-11e8-3e8a-31650f418d12"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[compat]
BioSequences = "~3.1.3"
Combinatorics = "~1.0.2"
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

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

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
# ╠═afe74614-ca56-11ed-0c86-45e6a934728e
# ╠═a1fa8018-d621-4751-8a7c-442108feced0
# ╠═1ff9ac81-4610-4e4a-8375-a3d023475353
# ╠═c505364d-ba7a-4a2b-aaba-3edd5947e667
# ╠═c023e2e1-23af-42a8-8e9e-bac435639b27
# ╠═f999b4e1-567e-463e-8db4-896c5a15c208
# ╠═44da0576-5626-4541-afaa-624e051cb8fe
# ╠═9edde78d-b1b0-4259-922d-0080bd5c0444
# ╠═705aba37-01cf-4d63-9f1b-479609552959
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
