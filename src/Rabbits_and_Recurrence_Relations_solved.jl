using Markdown, Test, BenchmarkTools
md"""
Given: Positive integers n≤40 and k≤5.
Return: The total number of rabbit pairs that will be present after n months, if we begin with 1 pair and in each generation, every pair of reproduction-age rabbits produces a litter of k rabbit pairs (instead of only 1 pair).
"""
function fibonacci_rabbit(n::Int, k::Int)
    if n < 3
        return 1
    end
    acc = acc1 = 1
    while n > 2
        acc, acc1 = acc + k * acc1, acc
        n -= 1
    end
    return acc
end
@test fibonacci_rabbit(5, 3) == 19
@test fibonacci_rabbit(1, 3) == 1
@test fibonacci_rabbit(2, 3) == 1
@benchmark fibonacci_rabbit(35, 2)