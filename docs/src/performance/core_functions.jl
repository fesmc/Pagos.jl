using Pagos, CairoMakie
arrhenius_rate_factor = ArrheniusRateFactor()
T_relative_kelvin = range(-50, stop = 10, step = 0.1) .+ 273.15
A = similar(T_relative_kelvin)
@b rate_factor($T_relative_kelvin, $arrhenius_rate_factor)
@b rate_factor!($A, $T_relative_kelvin, $arrhenius_rate_factor)