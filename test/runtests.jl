using PhyNEST
using Test

tests=[
    "test_check.jl"
    #"test_Quartets.jl"
    #"test_Probabilities.jl"
    #"test_Optimization.jl"
]

for t in tests
    include(t)
end