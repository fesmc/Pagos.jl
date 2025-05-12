@testset "exponential_vertical_layers" begin
    nz = 20
    sigma = exponential_vertical_layers(nz)
    for l in 1:nz-1
        if l == 1
            dsigma = sigma[1] 
        else
            dsigma = sigma[l] - sigma[l - 1]
        end
        dsigma_next = sigma[l + 1] - sigma[l]
        @test dsigma_next > dsigma
    end
end