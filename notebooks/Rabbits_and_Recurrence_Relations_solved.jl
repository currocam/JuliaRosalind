### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ fac28afc-723a-4266-bbd2-188dd9ae239a
using Test

# ╔═╡ effac4c2-c5ab-11ed-164b-5348b9c41a6b
md"""
## Rabbits and Recurrence Relations

A sequence is an ordered collection of objects (usually numbers), which are allowed to repeat. Sequences can be finite or infinite. Two examples are the finite sequence (π,−2–√,0,π)
and the infinite sequence of odd numbers (1,3,5,7,9,…). We use the notation an to represent the n -th term of a sequence.

A recurrence relation is a way of defining the terms of a sequence with respect to the values of previous terms. In the case of Fibonacci's rabbits from the introduction, any given month will contain the rabbits that were alive the previous month, plus any new offspring. A key observation is that the number of offspring in any month is equal to the number of rabbits that were alive two months prior. As a result, if Fn represents the number of rabbit pairs alive after the n-th month, then we obtain the Fibonacci sequence having terms Fn that are defined by the recurrence relation Fn=Fn−1+Fn−2 (with F1=F2=1

to initiate the sequence). Although the sequence bears Fibonacci's name, it was known to Indian mathematicians over two millennia ago.

When finding the n
-th term of a sequence defined by a recurrence relation, we can simply use the recurrence relation to generate terms for progressively larger values of n

. This problem introduces us to the computational technique of dynamic programming, which successively builds up solutions by using the answers to smaller cases.

Given: Positive integers n≤40
and k≤5

.

Return: The total number of rabbit pairs that will be present after n
months, if we begin with 1 pair and in each generation, every pair of reproduction-age rabbits produces a litter of k rabbit pairs (instead of only 1 pair).
"""

# ╔═╡ 444a68e0-3f08-49ee-a5dd-35aba303c99d
function fibonacci_rabbit(n::Int, k::Int)
    if n < 3 return 1 end
    acc = acc1 = 1
    while n > 2
        acc, acc1 = acc + k * acc1, acc
        n -= 1
    end
    return acc
end


# ╔═╡ 1e155166-ad64-40b3-bae3-16770d30e6e4
@test fibonacci_rabbit(5, 3) == 19

# ╔═╡ 1f6b1358-596e-446a-b469-2274ca4fc1f9
@test fibonacci_rabbit(1, 3) == 1

# ╔═╡ 1cfe09eb-37c5-4197-b153-bcdfd2d871cd
@test fibonacci_rabbit(2, 3) == 1

# ╔═╡ 5e5aadf4-81ad-456b-9230-f03fc841360e
fibonacci_rabbit(35, 2)

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
# ╠═effac4c2-c5ab-11ed-164b-5348b9c41a6b
# ╠═fac28afc-723a-4266-bbd2-188dd9ae239a
# ╠═444a68e0-3f08-49ee-a5dd-35aba303c99d
# ╠═1e155166-ad64-40b3-bae3-16770d30e6e4
# ╠═1f6b1358-596e-446a-b469-2274ca4fc1f9
# ╠═1cfe09eb-37c5-4197-b153-bcdfd2d871cd
# ╠═5e5aadf4-81ad-456b-9230-f03fc841360e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
