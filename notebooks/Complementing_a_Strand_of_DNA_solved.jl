### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ d1b0dcc0-078d-4b70-833f-ece47c443f85
using Test

# ╔═╡ a088ee32-c42f-11ed-1c29-b5d13a884476
md"""
# Complementing a Strand of DNA

In DNA strings, symbols 'A' and 'T' are complements of each other, as are 'C' and 'G'.

The reverse complement of a DNA string s
is the string sc formed by reversing the symbols of s

, then taking the complement of each symbol (e.g., the reverse complement of "GTCA" is "TGAC").

Given: A DNA string s of length at most 1000 bp.

Return: The reverse complement sc
of s.
"""

# ╔═╡ 8b69d3a8-c33e-4411-9822-e14eaf6b5874
function complementate(x::String)
    dict = Dict('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A')
    map(y -> dict[y], x) |> reverse
end

# ╔═╡ 23793c37-4f30-4937-a5e0-004f315a3016
complementate("AAAACCCGGT")

# ╔═╡ 1ebfb42f-9771-4e89-ab73-2c4cbd8ece57
@test complementate("AAAACCCGGT") == "ACCGGGTTTT"

# ╔═╡ 5c62cad4-bfdd-4844-8b3d-7e50d88f00fb
begin
    f = open("../datasets/rosalind_revc.txt")
    lines = readlines(f)
    sequence = lines[1]

end

# ╔═╡ 59646b4e-849a-4ff9-9cd3-f87bd5cf73e1
open("../outputs/rosalind_revc.txt", "w") do file
    write(file, complementate(sequence))
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
# ╠═a088ee32-c42f-11ed-1c29-b5d13a884476
# ╠═d1b0dcc0-078d-4b70-833f-ece47c443f85
# ╠═8b69d3a8-c33e-4411-9822-e14eaf6b5874
# ╠═23793c37-4f30-4937-a5e0-004f315a3016
# ╠═1ebfb42f-9771-4e89-ab73-2c4cbd8ece57
# ╠═5c62cad4-bfdd-4844-8b3d-7e50d88f00fb
# ╠═59646b4e-849a-4ff9-9cd3-f87bd5cf73e1
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
