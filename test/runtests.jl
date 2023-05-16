using ExtendRandomVariables
using Test

@testset "ExtendRandomVariables.jl" begin
    X1 = RV(Normal())
    X2 = RV(Exponential(3))
    X3 = RV(Poisson(3))
    X4 = RV(Poisson(4))

    @test 𝔼(X1)
    @test 𝔼(X3)
    @test 𝔼(X1|(X1<0))
    @test 𝔼(X3|(X3≤5))
    @test Z1 = X1 + X2
    @test Z2 = X3 + X4
    @test W1 = X1*X2
    @test W2 = X3*X4
    @test P(Z1 < 3)
    @test P(Z2 < 3)
    @test P(W1 < 3)
    @test P(W2 < 3)
    @test plotCDF(Z1)
    @test plotCDF!(Z1)
    @test plotPDF(Z1)
    @test plotPDF!(Z1)
end