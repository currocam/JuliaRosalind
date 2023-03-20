### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 3fcb541d-aceb-4648-a943-40b25f988021
using Test

# ╔═╡ 35f10cc8-bb22-439b-a983-1125848464d4
using BioSequences

# ╔═╡ 621524d6-c5d3-11ed-0412-33752871797c
md"""
Given two strings s and t, t is a substring of s if t is contained as a contiguous collection of symbols in s (as a result, t must be no longer than s).

The position of a symbol in a string is the total number of symbols found to its left, including itself (e.g., the positions of all occurrences of 'U' in "AUGCUUCAGAAAGGUCUUACG" are 2, 5, 6, 15, 17, and 18). The symbol at position i
of s is denoted by s[i].

A substring of s
can be represented as s[j:k], where j and k represent the starting and ending positions of the substring in s; for example, if s = "AUGCUUCAGAAAGGUCUUACG", then s[2:5]= "UGCU".

The location of a substring s[j:k] is its beginning position j; note that t will have multiple locations in s if it occurs more than once as a substring of s(see the Sample below).

Given: Two DNA strings s and t (each of length at most 1 kbp).

Return: All locations of t as a substring of s.
"""

# ╔═╡ 4a59a570-a8e0-40f0-b5e1-fbd43020389b
function find_motif(s, t)
	query = ExactSearchQuery(LongDNA{2}(t))
	[motif.start for motif in findall(query, LongDNA{2}(s))]
end

# ╔═╡ d352fc48-af67-4286-9c3e-bf4c8f189e6c
@test find_motif("GATATATGCATATACTT", "ATAT") == [2, 4, 10]

# ╔═╡ 22d99059-3969-4d60-a201-9b73116e0754
begin
	infile = open("../datasets/rosalind_subs.txt");
	lines = readlines(infile)
	s = lines[1]
	t = lines[2]
end


# ╔═╡ fd96e5a5-a100-4512-a213-2e23f4b15192
join(find_motif(s, t), " ")

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BioSequences = "7e6ae17a-c86d-528c-b3b9-7f778a29fe59"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[compat]
BioSequences = "~3.1.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0"
manifest_format = "2.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

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

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

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

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Twiddle]]
git-tree-sha1 = "29509c4862bfb5da9e76eb6937125ab93986270a"
uuid = "7200193e-83a8-5a55-b20d-5d36d44a0795"
version = "1.1.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
"""

# ╔═╡ Cell order:
# ╠═621524d6-c5d3-11ed-0412-33752871797c
# ╠═3fcb541d-aceb-4648-a943-40b25f988021
# ╠═35f10cc8-bb22-439b-a983-1125848464d4
# ╠═4a59a570-a8e0-40f0-b5e1-fbd43020389b
# ╠═d352fc48-af67-4286-9c3e-bf4c8f189e6c
# ╠═22d99059-3969-4d60-a201-9b73116e0754
# ╠═fd96e5a5-a100-4512-a213-2e23f4b15192
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
