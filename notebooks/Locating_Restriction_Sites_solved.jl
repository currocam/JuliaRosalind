using Markdown
md"""
Problem

Figure 2. Palindromic recognition site
A DNA string is a reverse palindrome if it is equal to its reverse complement. For instance, GCATGC is a reverse palindrome because its reverse complement is GCATGC. See Figure 2.

Given: A DNA string of length at most 1 kbp in FASTA format.

Return: The position and length of every reverse palindrome in the string having length between 4 and 12. You may return these pairs in any order.
"""

using BioSequences
using Test
using FASTX

is_palindrome = x -> x == reverse_complement(x)
function find_palindromes(x::LongDNA, min_length=4, max_length=12)
    n = length(x)
    # Create empty array of tuples
    result = Array{Tuple{Int, Int}, 1}()

    for length in max_length:-1:min_length
        for i in 1:n-length+1
            substring = x[i:i+length-1]
            if is_palindrome(substring)
                push!(result, (i, length))
            end
        end
    end
    return result
end

@test find_palindromes(dna"TCAATGCATGCGGGTCTATATGCAT") |> length == 8

infile = "datasets/rosalind_revp.txt"
## Open fasta file
sequences = [LongDNA{2}(sequence(record)) for record in open(FASTA.Reader, infile)]

palindromes = find_palindromes(sequences[1])
# Print into file
open("outputs/rosalind_revp.txt", "w") do io
    for (i, length) in palindromes
        println(io, "$i $length")
    end
end