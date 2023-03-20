### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 0cf3123c-509f-4182-a56d-a93de79e2b6c
using Test

# ╔═╡ 03fdd8a6-c629-11ed-1311-4338f52b7e08
md"""
# Calculating Expected Offspring

For a random variable X taking integer values between 1 and n, the expected value of X is E(X)=∑nk=1k×Pr(X=k). The expected value offers us a way of taking the long-term average of a random variable over a large number of trials.

As a motivating example, let X be the number on a six-sided die. Over a large number of rolls, we should expect to obtain an average of 3.5 on the die (even though it's not possible to roll a 3.5). The formula for expected value confirms that E(X)=∑6k=1k×Pr(X=k)=3.5.

More generally, a random variable for which every one of a number of equally spaced outcomes has the same probability is called a uniform random variable (in the die example, this "equal spacing" is equal to 1). We can generalize our die example to find that if X is a uniform random variable with minimum possible value a and maximum possible value b, then E(X)=a+b2. You may also wish to verify that for the dice example, if Y is the random variable associated with the outcome of a second die roll, then E(X+Y)=7.

Given: Six nonnegative integers, each of which does not exceed 20,000. The integers correspond to the number of couples in a population possessing each genotype pairing for a given factor. In order, the six given integers represent the number of couples having the following genotypes:

    AA-AA
    AA-Aa
    AA-aa
    Aa-Aa
    Aa-aa
    aa-aa

Return: The expected number of offspring displaying the dominant phenotype in the next generation, under the assumption that every couple has exactly two offspring.
"""

# ╔═╡ b4ff96c8-6278-4119-8086-0d27d781dad0
function expected_offspring(genotypes)
	probabilities = (1, 1, 1, 0.75, 0.5, 0)
	n_offpring = 2
	sum([n_offpring*k*p for (k, p) in zip(genotypes, probabilities)])
end

# ╔═╡ add6ea75-a24a-44ba-96c4-8ccd14a0fbdb
begin
	ϵ = 0.001
	genotypes = (1,0, 0, 1, 0, 1)
	@test abs(expected_offspring(genotypes) - 3.5) < ϵ
end

# ╔═╡ febfbd8d-dda4-4435-990b-2e993af7624a
genotypes_question = (18941, 18553, 17412, 17338, 18921, 19257)

# ╔═╡ 0d5db693-7740-41e2-9a27-958a764a367c
expected_offspring(genotypes_question)

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
# ╠═03fdd8a6-c629-11ed-1311-4338f52b7e08
# ╠═0cf3123c-509f-4182-a56d-a93de79e2b6c
# ╠═b4ff96c8-6278-4119-8086-0d27d781dad0
# ╠═add6ea75-a24a-44ba-96c4-8ccd14a0fbdb
# ╠═febfbd8d-dda4-4435-990b-2e993af7624a
# ╠═0d5db693-7740-41e2-9a27-958a764a367c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
