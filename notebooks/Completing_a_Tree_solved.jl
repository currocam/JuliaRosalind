using Markdown
md"""
Problem

Figure 2. A labeled tree with 6 vertices and 5 edges.
An undirected graph is connected if there is a path connecting any two nodes. A tree is a connected (undirected) graph containing no cycles; this definition forces the tree to have a branching structure organized around a central core of nodes, just like its living counterpart. See Figure 2.

We have already grown familiar with trees in “Mendel's First Law”, where we introduced the probability tree diagram to visualize the outcomes of a random variable.

In the creation of a phylogeny, taxa are encoded by the tree's leaves, or nodes having degree 1. A node of a tree having degree larger than 1 is called an internal node.

Given: A positive integer n
 (n≤1000
) and an adjacency list corresponding to a graph on n
 nodes that contains no cycles.

Return: The minimum number of edges that can be added to the graph to produce a tree.
Sample Dataset
10
1 2
2 8
4 10
5 9
6 10
7 9
Sample Output
3
"""
using Test
struct Tree
    n::Int
    adjacency_list::Array{Tuple{UInt64,UInt64}}
end


# Read data from file to struct
function read_tree(filename)
    open(filename) do file
        n = parse(UInt64, readline(file))
        adjacency_list = []
        for line in eachline(file)
            edge = parse.(UInt64, split(line))
            push!(adjacency_list, (edge[1], edge[2]))

        end
        return Tree(n, adjacency_list)
    end
end

count_missing_edges = x -> x.n - length(x.adjacency_list) - 1

example_tree = Tree(10, [(1, 2), (2, 8), (4, 10), (5, 9), (6, 10), (7, 9)])


@test count_missing_edges(example_tree) == 3

@time "datasets/rosalind_tree.txt" |> read_tree |> count_missing_edges
