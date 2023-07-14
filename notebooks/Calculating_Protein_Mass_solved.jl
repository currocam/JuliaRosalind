### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 3b32bb6d-860e-4d94-83e4-f4568b4ee3c5
using Test

# ╔═╡ 5e5686bd-48c0-4fe0-be60-c691df5586e6
using BioSequences

# ╔═╡ 9df067d4-c9b8-11ed-141c-ef74eeea7dd0
md"""
In a weighted alphabet, every symbol is assigned a positive real number called a weight. A string formed from a weighted alphabet is called a weighted string, and its weight is equal to the sum of the weights of its symbols.

The standard weight assigned to each member of the 20-symbol amino acid alphabet is the monoisotopic mass of the corresponding amino acid.

Given: A protein string P

of length at most 1000 aa.

Return: The total weight of P
. Consult the monoisotopic mass table.
"""

# ╔═╡ 8420ad83-d7ba-4cf8-b9c6-be3e2c70bc32
mass_table = Dict(
    AA_A => 71.03711,
    AA_C => 103.00919,
    AA_D => 115.02694,
    AA_E => 129.04259,
    AA_F => 147.06841,
    AA_G => 57.02146,
    AA_H => 137.05891,
    AA_I => 113.08406,
    AA_K => 128.09496,
    AA_L => 113.08406,
    AA_M => 131.04049,
    AA_N => 114.04293,
    AA_P => 97.05276,
    AA_Q => 128.05858,
    AA_R => 156.10111,
    AA_S => 87.03203,
    AA_T => 101.04768,
    AA_V => 99.06841,
    AA_W => 186.07931,
    AA_Y => 163.06333,
)

# ╔═╡ 527395b3-c26d-4ec2-a02b-dee49d4ee7f4
function calculate_mass(x, table)
    sum(mass_table[char] for char in x)
end

# ╔═╡ 6d0a6eec-382f-40aa-b228-c6b7dedbad0d
begin
    ϵ = 0.0001
    @test calculate_mass(aa"SKADYEK", mass_table) - 821.392 < ϵ
end

# ╔═╡ e4a952fa-600d-4b6c-b3cd-299cd372e966
x = aa"SKADYEK"

# ╔═╡ 2b033f4d-5a83-42fd-905c-28c148c82c77
begin
    f = open("../datasets/rosalind_prtm.txt")
    lines = readlines(f)
    sequence = LongAA(lines[1])
end

# ╔═╡ f548976e-509c-4bbc-bb2b-f785e16aff19
calculate_mass(sequence, mass_table)

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
# ╠═9df067d4-c9b8-11ed-141c-ef74eeea7dd0
# ╠═3b32bb6d-860e-4d94-83e4-f4568b4ee3c5
# ╠═5e5686bd-48c0-4fe0-be60-c691df5586e6
# ╠═8420ad83-d7ba-4cf8-b9c6-be3e2c70bc32
# ╠═527395b3-c26d-4ec2-a02b-dee49d4ee7f4
# ╠═6d0a6eec-382f-40aa-b228-c6b7dedbad0d
# ╠═e4a952fa-600d-4b6c-b3cd-299cd372e966
# ╠═2b033f4d-5a83-42fd-905c-28c148c82c77
# ╠═f548976e-509c-4bbc-bb2b-f785e16aff19
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
