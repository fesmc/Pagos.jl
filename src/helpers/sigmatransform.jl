sigma(z, b, H) = (z - b) / H

function exponential_vertical_layers(n)
    dz = 1 / n
    return (exp.(dz:dz:1) .- exp(0)) ./ (exp(1) - exp(0))
end