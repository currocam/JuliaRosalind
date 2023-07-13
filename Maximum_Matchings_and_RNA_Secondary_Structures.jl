using Markdown
using Test
using BenchmarkTools
using BioSequences
md"""
The graph theoretical analogue of the quandary stated in the introduction above is that if we have an RNA string s that does not have the same number of occurrences of 'C' as 'G' and the same number of occurrences of 'A' as 'U', then the bonding graph of s cannot possibly possess a perfect matching among its basepair edges. For example, see Figure 1; in fact, most bonding graphs will not contain a perfect matching.
  
In light of this fact, we define a maximum matching in a graph as a matching containing as many edges as possible. See Figure 2 for three maximum matchings in graphs.
  
A maximum matching of basepair edges will correspond to a way of forming as many base pairs as possible in an RNA string, as shown in Figure 3.
  
Given: An RNA string s of length at most 100.

Return: The total possible number of maximum matchings of basepair edges in the bonding graph of s.
"""
import StatsBase: countmap

function n_maximum_matching(x::LongRNA)
  counts = countmap(x)
  n = BigInt(1)
  n *= max(counts[RNA_A], counts[RNA_U]) |> big |> factorial
  n ÷= abs(counts[RNA_A] - counts[RNA_U]) |> big |> factorial
  n *= max(counts[RNA_C], counts[RNA_G]) |> big |> factorial
  n ÷= abs(counts[RNA_C] - counts[RNA_G]) |> big |> factorial
  return n
end

@test n_maximum_matching(rna"AUGCUUC") == 6

x = rna"UGCAAAUCAAUACGGCGCCGUCACGGAUGCCCAGUUGUCAUUGCUCUAAGCCCUCGAAUAGGACAAGCGAAAUGUCUAUCCAGCACCG"
@time n_maximum_matching(x)

@benchmark n_maximum_matching(x)