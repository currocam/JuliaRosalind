### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 25f6b5f1-50dc-4d18-8e61-8e0ba3ebe9b4
using Test;

# ╔═╡ bf7e9ec5-3208-445a-b044-90ac46fa0de3
using StatsBase

# ╔═╡ 61bce374-55ae-46af-92e3-a48451b022f4
using BioSequences

# ╔═╡ cf9e2946-c8d0-11ed-24a3-7b48dcbfd4a9
md"""
For positive integers a and n, a modulo n (written amodn in shorthand) is the remainder when a is divided by n. For example, 29mod11=7 because 29=11×2+7

.

Modular arithmetic is the study of addition, subtraction, multiplication, and division with respect to the modulo operation. We say that a
and b are congruent modulo n if amodn=bmodn; in this case, we use the notation a≡bmodn

.

Two useful facts in modular arithmetic are that if a≡bmodn
and c≡dmodn, then a+c≡b+dmodn and a×c≡b×dmodn. To check your understanding of these rules, you may wish to verify these relationships for a=29, b=73, c=10, d=32, and n=11

.

As you will see in this exercise, some Rosalind problems will ask for a (very large) integer solution modulo a smaller number to avoid the computational pitfalls that arise with storing such large numbers.

Given: A protein string of length at most 1000 aa.

Return: The total number of different RNA strings from which the protein could have been translated, modulo 1,000,000. (Don't neglect the importance of the stop codon in protein translation.)
"""

# ╔═╡ aa719df0-ad5b-4277-a232-c7d0b3bd51d2
standard_trans = ncbi_trans_table[1] |> values |> countmap

# ╔═╡ 21bddb24-cff4-4220-8d99-078b54d9fe82
keys(standard_trans)

# ╔═╡ 4f2304a1-753d-45f3-9e54-1231387aa0e6
function infer_mRNA(x)
    n = 1000000
    counter = standard_trans[AA_Term]
    for aa in LongAA(x)
        counter = mod(counter * standard_trans[aa], n)
    end
    counter
end

# ╔═╡ b197b6bd-17d0-4542-a22b-b846e783da94
@test infer_mRNA("MA") == 12

# ╔═╡ 79694ae6-ab0c-4504-bf67-4c4459d6111a
x = "MFAWIKAEKVRQFHVAPCCLRLEGLYNFIVHSGKKPSNEIRCQCSTWENNSQKWPCPCFICLMPFQMKVRQQGYEQWQIETQTMYITSGSNVDCLQIANSWIHARTMHDVQDGYCVQERKKKQNRSDHLQWDCMRKHIWICWSENDIHFHHNAIGGVGFQGTDSIIEMHSMRYIADYLHDYLDCTKRIQDFPNYFRWKSQIEQAFNWYDMTRENMVHVRCFWSFTDSGEQMDYATPSFRVKCTEFEKNDLNSAEIWKFYTLGPSVYRCTHSNSQRNKGSVWWFIFRGKTMSMYECMCKIVRAGKVFHEYFYKEWNPFMIHCQGRQHCFWWMQVCDEIKCGCNAHLKFEVRQYLSKQKGLRCQREFKTPNWCCSWQYWIPTWSKQRDKIYCGFHMNMDWNLCPALRIVLTCPVHHDCFRVCDVSSCQMQGNAYDGGCNHPWQGGPSAQRDNKMFLWLISAKPDNWGHYHPGFMCKTPNTHVASCKCGMAMQTETATDSSLHVSFQGHLLKDQMHRKHGVTAKHLHAGFYTMMLVPQMEATRMQNQGWRYMGLGPMSNRIAWRACVAEGKAKPWYFNLGWCSKFLIKNETANTQWQFNNQACYDFPMLCPNWVGPYVEKGQIAYMWQPMMIHYRVIAIGMIRSICVFYVARAAKSDCREESGFEDIIDVWDAMNYKDTLRFRRMENTRDQLQTWMEINFVSPTSELNEYAWSKCNNCAQANDSRMEAAGEDECVPDSSNPRAECWYVFPASCFQQNPLWEQPYRCNLEEMQMKLELSPPIGYGVRQFQPMNAIMYENAYAHNDPNKHAKQTQHDLAIQNAYCSPMIQAHDDWQCCRFDYRWHAGTLANSIPQGKGGMIDQCLGFFDAARWETNSFELNKHWFQIPHPPMSITYESLTQVIENCMYRTRMVELGLWWQKQRILLLRFCYMRRVEGIWEILAVHGAPAPLWETAVHVFVWLWMHKPMSGIVAALECQHFRHHQAHWMTPQRATNCWSAWAKY"

# ╔═╡ 9b9fbe8b-5926-4b6d-9651-de048d30c930
infer_mRNA(x)

# ╔═╡ 6d23929b-0de0-44c4-bc32-4c3b5ee01750
md"""
Real answer using BigInt
"""

# ╔═╡ f8e104fb-d77a-437f-9359-b1d16dacf69e
begin
    n = 1000000
    counter = BigInt(standard_trans[AA_Term])
    for aa in LongAA(x)
        counter = counter * standard_trans[aa]
    end
    counter
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BioSequences = "7e6ae17a-c86d-528c-b3b9-7f778a29fe59"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[compat]
BioSequences = "~3.1.3"
StatsBase = "~0.33.21"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0"
manifest_format = "2.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

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

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "485193efd2176b88e6622a39a246f8c5b600e74e"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.6"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

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

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

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

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
"""

# ╔═╡ Cell order:
# ╠═cf9e2946-c8d0-11ed-24a3-7b48dcbfd4a9
# ╠═25f6b5f1-50dc-4d18-8e61-8e0ba3ebe9b4
# ╠═bf7e9ec5-3208-445a-b044-90ac46fa0de3
# ╠═aa719df0-ad5b-4277-a232-c7d0b3bd51d2
# ╠═21bddb24-cff4-4220-8d99-078b54d9fe82
# ╠═61bce374-55ae-46af-92e3-a48451b022f4
# ╠═4f2304a1-753d-45f3-9e54-1231387aa0e6
# ╠═b197b6bd-17d0-4542-a22b-b846e783da94
# ╠═79694ae6-ab0c-4504-bf67-4c4459d6111a
# ╠═9b9fbe8b-5926-4b6d-9651-de048d30c930
# ╠═6d23929b-0de0-44c4-bc32-4c3b5ee01750
# ╠═f8e104fb-d77a-437f-9359-b1d16dacf69e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
