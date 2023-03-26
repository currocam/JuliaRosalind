### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 8d033401-1b57-4839-a560-1229bcc18618
using BioSequences

# ╔═╡ 3a455006-33aa-4743-bf98-b0be83701599
using Test

# ╔═╡ 4eb27171-083e-48ec-908d-62daf1d025a4
using DelimitedFiles

# ╔═╡ 4d2c3f74-dffe-4ed4-9d98-f02ccee49c20
using FASTX

# ╔═╡ 90585e50-cb26-11ed-355a-afecd43116cd
md"""
A matrix is a rectangular table of values divided into rows and columns. An m×n matrix has m rows and n columns. Given a matrix A, we write Ai,j to indicate the value found at the intersection of row i and column j

.

Say that we have a collection of DNA strings, all having the same length n
. Their profile matrix is a 4×n matrix P in which P1,j represents the number of times that 'A' occurs in the jth position of one of the strings, P2,j represents the number of times that C occurs in the j

th position, and so on (see below).

A consensus string c
is a string of length n formed from our collection by taking the most common symbol at each position; the jth symbol of c therefore corresponds to the symbol having the maximum value in the j

-th column of the profile matrix. Of course, there may be more than one most common symbol, leading to multiple possible consensus strings.

	A T C C A G C T
	G G G C A A C T
	A T G G A T C T
DNA Strings	A A G C A A C C
	T T G G A A C T
	A T G C C A T T
	A T G G C A C T
	A   5 1 0 0 5 5 0 0
Profile	C   0 0 1 4 2 0 6 1
	G   1 1 6 3 0 1 0 0
	T   1 5 0 0 0 1 1 6
Consensus	A T G C A A C T

Given: A collection of at most 10 DNA strings of equal length (at most 1 kbp) in FASTA format.

Return: A consensus string and profile matrix for the collection. (If several possible consensus strings exist, then you may return any one of them.)
"""

# ╔═╡ ec53ac85-fe1f-4d76-a508-399ff710e939
records_sample = FASTA.Record.("seq", ["ATCCAGCT", "GGGCAACT", "ATGGATCT", "AAGCAACC","TTGGAACT", "ATGCCATT", "ATGGCACT"])

# ╔═╡ a34d02ab-129e-4352-a3dc-340fe22742b3
function profile_matrix(r)
	records = [LongDNA{4}(sequence(x)) for x in r]
	n = length(records[1])
	profile = zeros(Int, 4, n)
	for x in records
	for (i, char) in enumerate(x)
		if char== DNA_A profile[1, i]+=1 end 
		if char== DNA_C profile[2, i]+=1 end 
		if char== DNA_G profile[3, i]+=1 end 
		if char== DNA_T profile[4, i]+=1 end 
	end end
	profile
end

# ╔═╡ 004997a1-e07a-4b10-8ee3-3f7e87a5e132
function consensus_sequence(profile)
	dict = Dict(1=>DNA_A, 2=>DNA_C, 3=>DNA_G, 4=>DNA_T)
	LongDNA{4}([dict[argmax(col)] for col in eachcol(profile)])
end

# ╔═╡ 85de25cd-b4d2-44d6-b465-b9a272db725b
@test records_sample|> profile_matrix |> consensus_sequence |> String == "ATGCAACT"

# ╔═╡ 3a8790ca-75b5-4351-8152-1a885725a56f
open(FASTA.Reader, "../datasets/rosalind_cons.txt") do reader
global records = [r for r in reader]
end

# ╔═╡ 023ee250-3665-47a2-b545-81eea0a35044
records

# ╔═╡ 593cb794-2435-459c-ab97-48c3dfc8c43b
profile = profile_matrix(records)

# ╔═╡ 69bcfaa9-b3ad-4974-ba52-9b764ef60d20
consensus = profile |> consensus_sequence |> String

# ╔═╡ 07716619-b165-4b9f-950d-ec5c199f261b
open("../outputs/rosalind_cons.txt", "w") do io
	println(io, consensus)
	writedlm(io, profile, ' ')
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BioSequences = "7e6ae17a-c86d-528c-b3b9-7f778a29fe59"
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
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

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

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

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

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
# ╠═90585e50-cb26-11ed-355a-afecd43116cd
# ╠═8d033401-1b57-4839-a560-1229bcc18618
# ╠═3a455006-33aa-4743-bf98-b0be83701599
# ╠═4eb27171-083e-48ec-908d-62daf1d025a4
# ╠═ec53ac85-fe1f-4d76-a508-399ff710e939
# ╠═a34d02ab-129e-4352-a3dc-340fe22742b3
# ╠═004997a1-e07a-4b10-8ee3-3f7e87a5e132
# ╠═85de25cd-b4d2-44d6-b465-b9a272db725b
# ╠═4d2c3f74-dffe-4ed4-9d98-f02ccee49c20
# ╠═3a8790ca-75b5-4351-8152-1a885725a56f
# ╠═023ee250-3665-47a2-b545-81eea0a35044
# ╠═593cb794-2435-459c-ab97-48c3dfc8c43b
# ╠═69bcfaa9-b3ad-4974-ba52-9b764ef60d20
# ╠═07716619-b165-4b9f-950d-ec5c199f261b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
