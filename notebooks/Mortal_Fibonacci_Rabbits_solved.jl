### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 460fb007-c2ab-4857-b4f4-b76ff7c2a256
using Test;

# ╔═╡ d52d1690-c687-11ed-02d3-d5d6b6aeb38e
md"""
Recall the definition of the Fibonacci numbers from “Rabbits and Recurrence Relations”, which followed the recurrence relation Fn=Fn−1+Fn−2

and assumed that each pair of rabbits reaches maturity in one month and produces a single pair of offspring (one male, one female) each subsequent month.

Our aim is to somehow modify this recurrence relation to achieve a dynamic programming solution in the case that all rabbits die out after a fixed number of months. See Figure 4 for a depiction of a rabbit tree in which rabbits live for three months (meaning that they reproduce only twice before dying).

Given: Positive integers n≤100
and m≤20.

Return: The total number of pairs of rabbits that will remain after the n
-th month if all rabbits live for m months.

37395433272364347
"""

# ╔═╡ 663784fd-af31-416f-8d77-bff1a1c53014
function mortal_fibonacci(n, m)
	if n < 3 return 1  end
	if n == 3 return 2  end	
	alive = Vector{Int}(undef,n)
	alive[1]= 1
	alive[2]=1
	for i in range(3, m+1)
		alive[i] = alive[i - 1] + alive[i - 2]
	end
	alive[m+1] -= 1 
	for i in range(m+2, n)
		alive[i] = alive[i - 1] + alive[i - 2] - alive[i - m - 1]
	end		
	alive[n]
end

# ╔═╡ 6ec985f5-320b-4e07-8653-4f68b72e0370
begin
	rabbits = [1, 1, 2, 2, 3, 4]
	for (month, value) in enumerate(rabbits)
	           @test mortal_fibonacci(month, 3) == value
	       end
end

# ╔═╡ da13fb36-a8d2-4ef1-8b8e-c76162f1dcad
@test mortal_fibonacci(82, 19) == 61115936848684058

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
# ╠═d52d1690-c687-11ed-02d3-d5d6b6aeb38e
# ╠═460fb007-c2ab-4857-b4f4-b76ff7c2a256
# ╠═663784fd-af31-416f-8d77-bff1a1c53014
# ╠═6ec985f5-320b-4e07-8653-4f68b72e0370
# ╠═da13fb36-a8d2-4ef1-8b8e-c76162f1dcad
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
