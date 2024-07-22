# Module to solve the 1D enthalpy equation for an ice column

module IceEnthalpy

mutable struct icecolumn_params{T<:AbstractFloat}
    nz::Integer
    omega_max::T    # [-] Maximum allowed water fraction inside ice, typically omega_max=0.02 
    T0::T           # [K or degreesCelcius] Reference melting temperature  
    rho_ice::T 
    rho_w::T
    L_ice::T
    sec_year::T
end

mutable struct icecolumn{T<:AbstractFloat}
    ζ::Vector{T}
    ζac::Vector{T}
    H::Vector{T}
    cp::Vector{T}
    kt::Vector{T}
    Qxy::Vector{T}
    uz::Vector{T}
    T_pmp::Vector{T}
    Q_strn::Vector{T}
    Q_b::T
    Q_rock::T
    bmb_grnd::T
    H_ice::T
    H_w::T
    f_grnd::T
    T_srf::T
    T_shlf::T

    dζa::Vector{T}
    dζb::Vector{T}

    # Diagnostic output variables
    T::Vector{T}
    ω::Vector{T}
    Q_ice_b::T
    H_cts::T
end

end