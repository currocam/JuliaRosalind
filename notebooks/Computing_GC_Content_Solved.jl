### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 1a760eb7-1c2a-4af7-bb6d-25061fbe7808
using Test

# ╔═╡ 54b41d24-9c85-416d-b308-9894022c462f
using FASTX

# ╔═╡ a549e3f3-4b56-4f7c-ab93-4e3261d9adef
using BioSequences

# ╔═╡ 98969194-c5b9-11ed-2261-0fd636e41342
md"""
Problem

The GC-content of a DNA string is given by the percentage of symbols in the string that are 'C' or 'G'. For example, the GC-content of "AGCTATAG" is 37.5%. Note that the reverse complement of any DNA string has the same GC-content.

DNA strings must be labeled when they are consolidated into a database. A commonly used method of string labeling is called FASTA format. In this format, the string is introduced by a line that begins with '>', followed by some labeling information. Subsequent lines contain the string itself; the first line to begin with '>' indicates the label of the next string.

In Rosalind's implementation, a string in FASTA format will be labeled by the ID "Rosalind_xxxx", where "xxxx" denotes a four-digit code between 0000 and 9999.

Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).

Return: The ID of the string having the highest GC-content, followed by the GC-content of that string. Rosalind allows for a default error of 0.001 in all decimal answers unless otherwise stated; please see the note on absolute error below.
"""

# ╔═╡ 87dca189-c323-4bbe-bf90-0d5b56b55aa7
function calculate_gc_content(x)
	count(char->char=='G' || char=='C', x) / length(x)
end

# ╔═╡ 80e3e391-80d0-4c50-836c-f500dccc4714
begin
	ϵ = 0.001
	@test abs(calculate_gc_content("AGCTATAG") - 0.375) < ϵ 
end

# ╔═╡ 5c91a41e-ff9a-42b6-af67-dc381fca5c7d
function find_largest_gc(reader)
		best = ("", -Inf)
	    for record in reader
        	gc = record |> FASTA.sequence|> calculate_gc_content
			if gc > best[2]
				best = (FASTA.identifier(record), gc)
			end
    	end
	best
end

# ╔═╡ 323b808e-b22f-4032-8f20-c9d07dbbaa59
begin
	headers = ["Rosalind_6404", "Rosalind_5959", "Rosalind_0808"]
	sequences = [
		dna"CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG",
		dna"CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC",
		dna"CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT"
	]
end


# ╔═╡ 434652b8-0e2c-4a18-8ec3-060739a8a54d
begin
	records = [FASTA.Record(x[1], x[2]) for x in zip(headers, sequences)]
	@test find_largest_gc(records)[1] == "Rosalind_0808"
	@test abs(find_largest_gc(records)[2] - 0.60919540) < ϵ 
end

# ╔═╡ fd832fd0-6e06-4b8e-9841-9b8e639e8c3f
begin
	infile = "../datasets/rosalind_gc.txt"
	open(FASTA.Reader, infile) |> 
		find_largest_gc |>
		(x-> (x[1], 100x[2]))
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BioSequences = "7e6ae17a-c86d-528c-b3b9-7f778a29fe59"
FASTX = "c2308a5c-f048-11e8-3e8a-31650f418d12"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

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
# ╠═98969194-c5b9-11ed-2261-0fd636e41342
# ╠═1a760eb7-1c2a-4af7-bb6d-25061fbe7808
# ╠═54b41d24-9c85-416d-b308-9894022c462f
# ╠═a549e3f3-4b56-4f7c-ab93-4e3261d9adef
# ╠═87dca189-c323-4bbe-bf90-0d5b56b55aa7
# ╠═80e3e391-80d0-4c50-836c-f500dccc4714
# ╠═5c91a41e-ff9a-42b6-af67-dc381fca5c7d
# ╟─323b808e-b22f-4032-8f20-c9d07dbbaa59
# ╠═434652b8-0e2c-4a18-8ec3-060739a8a54d
# ╠═fd832fd0-6e06-4b8e-9841-9b8e639e8c3f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
