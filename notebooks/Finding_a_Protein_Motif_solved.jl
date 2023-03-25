### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ f7b84209-07e8-4a2d-84a7-94f8999a2924
using FASTX

# ╔═╡ 4d271bdd-5fec-4b62-bf65-872fca3dfa2f
using HTTP

# ╔═╡ 9d3a62a8-c209-4cf3-8770-abe8d25c7536
using BioSequences

# ╔═╡ b9cc5dae-d52c-4f6e-9c64-fdb389bd3788
using Test

# ╔═╡ 4170c2f7-13e6-4570-9ab0-72592a877386
function get_url(id::String)
	x = split(id, "_")[1]
	"http://www.uniprot.org/uniprot/$x.fasta"
end

# ╔═╡ aa2e1b9d-6351-4673-9951-943ba7e50c1f
begin
	@test get_url("uniprot_id") == "http://www.uniprot.org/uniprot/uniprot.fasta"
	@test get_url("uniprotid") == "http://www.uniprot.org/uniprot/uniprotid.fasta"
end

# ╔═╡ 364f5e5b-e421-4305-abff-96dc8de076ce
function find_motif(motif, sequence)
	[x.captured[1] for x in eachmatch(motif, sequence)]
end

# ╔═╡ 2fe6015e-155e-489f-9d80-d905e7d225fb
@test find_motif(prosite"N{P}[ST]{P}", aa"NSSSPPNSSS") == [1, 7]

# ╔═╡ a16d7edf-eb13-4b4a-9cb1-fc6ee1103ec5
samples = ["A2Z669","B5ZC00","P07204_TRBM_HUMAN","P20840_SAG1_YEAST"]

# ╔═╡ 946b87b7-9371-4552-bbba-fb767c5cb9b9
infiles = map(x -> get_url(x) |> download |> open, samples)

# ╔═╡ 695cfac9-8686-455e-8375-f93920fb5f3c
records = FASTAReader.(first, infiles)

# ╔═╡ a9a64b2c-770f-434b-89a0-3a398f30c173
begin
	observed = [find_motif(prosite"N{P}[ST]{P}",  LongAA(sequence(x))) for x in records]
	expected = [[], [85,118,142, 306, 395], [47,115, 116, 382, 409], [79 ,109 ,135, 248, 306,348 ,364 ,402 ,485 ,501 ,614] ]
	@test observed == expected
end

# ╔═╡ 4f7a331f-09bd-4660-88af-d3f44cd3401b
md"""
# Solution
"""

# ╔═╡ ca7c2512-adc1-4734-9b66-269efb7ef919
samples2 = [
	"A5F5B4","P46096_SYT1_MOUSE","Q9CE42","P0A4Y7","P04141_CSF2_HUMAN",
	"Q7S432","Q4FZD7","P04921_GLPC_HUMAN","Q1JLH6","P12763_A2HS_BOVIN",
	"Q1JHI2","P07204_TRBM_HUMAN","P47002"]

# ╔═╡ 2c1201bb-9ab1-46a1-9c51-c59361fa65d3
infiles2 = map(x -> get_url(x) |> download |> open, samples2)

# ╔═╡ d4759236-6f9c-430e-a0d9-b42c2ff81a1a
records2 = FASTAReader.(first, infiles2)

# ╔═╡ 8b836c7e-786c-40f0-892a-f6d366b8cbb0
begin
	open("../outputs/mprt.txt","w") do io
	for (x, y) in zip(records2, samples2)
		pos = find_motif(prosite"N{P}[ST]{P}",  LongAA(sequence(x)))
		if isempty(pos) continue end
		println(io,y)
		println(io,join(pos, " "))
	end
	end
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BioSequences = "7e6ae17a-c86d-528c-b3b9-7f778a29fe59"
FASTX = "c2308a5c-f048-11e8-3e8a-31650f418d12"
HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[compat]
BioSequences = "~3.1.3"
FASTX = "~2.0.1"
HTTP = "~1.7.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0"
manifest_format = "2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

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

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FASTX]]
deps = ["Automa", "BioGenerics", "BioSequences", "ScanByte", "SnoopPrecompile", "StringViews", "TranscodingStreams"]
git-tree-sha1 = "d9b3f4020a5afda65de0db1ee150409c1b69f44c"
uuid = "c2308a5c-f048-11e8-3e8a-31650f418d12"
version = "2.0.1"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "37e4657cd56b11abe3d10cd4a1ec5fbdb4180263"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.7.4"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "6503b77492fd7fcb9379bf73cd31035670e3c509"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9ff31d101d987eb9d66bd8b176ac7c277beccd09"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.20+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

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

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.StringViews]]
git-tree-sha1 = "6f58c457b8adeab3ac97cca6b3434ffefd7bd71e"
uuid = "354b36f9-a18e-4713-926e-db85100087ba"
version = "1.3.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

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

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═f7b84209-07e8-4a2d-84a7-94f8999a2924
# ╠═4d271bdd-5fec-4b62-bf65-872fca3dfa2f
# ╠═9d3a62a8-c209-4cf3-8770-abe8d25c7536
# ╠═b9cc5dae-d52c-4f6e-9c64-fdb389bd3788
# ╠═4170c2f7-13e6-4570-9ab0-72592a877386
# ╠═aa2e1b9d-6351-4673-9951-943ba7e50c1f
# ╠═364f5e5b-e421-4305-abff-96dc8de076ce
# ╠═2fe6015e-155e-489f-9d80-d905e7d225fb
# ╠═a16d7edf-eb13-4b4a-9cb1-fc6ee1103ec5
# ╠═946b87b7-9371-4552-bbba-fb767c5cb9b9
# ╠═695cfac9-8686-455e-8375-f93920fb5f3c
# ╠═a9a64b2c-770f-434b-89a0-3a398f30c173
# ╠═4f7a331f-09bd-4660-88af-d3f44cd3401b
# ╠═ca7c2512-adc1-4734-9b66-269efb7ef919
# ╠═2c1201bb-9ab1-46a1-9c51-c59361fa65d3
# ╠═d4759236-6f9c-430e-a0d9-b42c2ff81a1a
# ╠═8b836c7e-786c-40f0-892a-f6d366b8cbb0
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
