### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 776ab84d-bf85-4a01-baab-9226d7306187
using Test

# ╔═╡ a7d315ae-69fd-4fcb-8aa9-e4b5990464d4
using BioSequences

# ╔═╡ 60f2cebe-c5c4-11ed-2c46-dba0a2814dea
md"""
# Translating RNA into Protein
The 20 commonly occurring amino acids are abbreviated by using 20 letters from the English alphabet (all letters except for B, J, O, U, X, and Z). Protein strings are constructed from these 20 symbols. Henceforth, the term genetic string will incorporate protein strings along with DNA strings and RNA strings.

The RNA codon table dictates the details regarding the encoding of specific codons into the amino acid alphabet.

Given: An RNA string s

corresponding to a strand of mRNA (of length at most 10 kbp).

Return: The protein string encoded by s
.
"""

# ╔═╡ 065318bc-3f3f-457a-a749-604415d19d83
function translate2string(x)
    x |> LongRNA{2} |> translate |> String |> y -> split(y, "*")[1]
end

# ╔═╡ dfe6f738-d735-4d89-877c-f31763b75867
@test translate2string("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA") ==
      "MAMAPRTEINSTRING"

# ╔═╡ 5a22c898-df41-47fe-bf0d-eac5925c764c
begin
    infile = open("../datasets/rosalind_prot.txt")
    sequence = readlines(infile)[1]
end

# ╔═╡ a29d775b-fdec-4b0d-9a5e-7d82a322c891
open("../outputs/rosalind_prton.txt", "w") do outfile
    write(outfile, translate2string(sequence))
end

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
# ╠═60f2cebe-c5c4-11ed-2c46-dba0a2814dea
# ╠═776ab84d-bf85-4a01-baab-9226d7306187
# ╠═a7d315ae-69fd-4fcb-8aa9-e4b5990464d4
# ╠═065318bc-3f3f-457a-a749-604415d19d83
# ╠═dfe6f738-d735-4d89-877c-f31763b75867
# ╠═5a22c898-df41-47fe-bf0d-eac5925c764c
# ╠═a29d775b-fdec-4b0d-9a5e-7d82a322c891
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
