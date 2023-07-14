### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 79506edc-7177-4532-9e0e-a1cfad281fb5
using Test

# ╔═╡ f981820a-c42b-11ed-1bc5-d3e337d80440
md"""
# Transcribing DNA into RNA

An RNA string is a string formed from the alphabet containing 'A', 'C', 'G', and 'U'.

Given a DNA string t
corresponding to a coding strand, its transcribed RNA string u is formed by replacing all occurrences of 'T' in t with 'U' in u

.

Given: A DNA string t

having length at most 1000 nt.

Return: The transcribed RNA string of t
.
"""

# ╔═╡ 06214641-de7c-45c7-bb5b-b2dd8a68caa8
function transcribe(dna::String)
    replace(dna, "T" => "U")
end

# ╔═╡ 23927e5a-e18c-4172-92bf-3d69cf2fd110
@test transcribe("GATGGAACTTGACTACGTAAATT") == "GAUGGAACUUGACUACGUAAAUU"

# ╔═╡ 0aa14549-2bd2-4849-b5c6-f763e99dc7db
begin
    f = open("../datasets/rosalind_rna.txt")
    lines = readlines(f)
    sequence = lines[1]
end

# ╔═╡ 1c81b55c-87e9-440c-acb6-3a12af14654f
transcribe(sequence)

# ╔═╡ 6c33e8ae-d6e1-4dd0-8819-635d52f7b5ff
open("../outputs/rosalind_rna.txt", "w") do file
    write(file, transcribe(sequence))
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0"
manifest_format = "2.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
"""

# ╔═╡ Cell order:
# ╠═f981820a-c42b-11ed-1bc5-d3e337d80440
# ╠═79506edc-7177-4532-9e0e-a1cfad281fb5
# ╠═06214641-de7c-45c7-bb5b-b2dd8a68caa8
# ╠═23927e5a-e18c-4172-92bf-3d69cf2fd110
# ╠═0aa14549-2bd2-4849-b5c6-f763e99dc7db
# ╠═1c81b55c-87e9-440c-acb6-3a12af14654f
# ╠═6c33e8ae-d6e1-4dd0-8819-635d52f7b5ff
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
