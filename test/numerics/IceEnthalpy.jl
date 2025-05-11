# Module to solve the 1D enthalpy equation for an ice column

module IceEnthalpy

export IceColumnParams
export IceColumn

@kwdef struct IceColumnParams{T<:AbstractFloat}
    nz::Integer
    omega_max::T    # [-] Maximum allowed water fraction inside ice, typically omega_max=0.02 
    T0::T           # [K or degreesCelcius] Reference melting temperature  
    rho_ice::T 
    rho_w::T
    L_ice::T
    sec_year::T
end

mutable struct IceColumn{T<:AbstractFloat}
    ζ::Vector{T}
    ζac::Vector{T}
    dζa::Vector{T}
    dζb::Vector{T}

    # H::Vector{T}
    # cp::Vector{T}
    # kt::Vector{T}
    # Qxy::Vector{T}
    # uz::Vector{T}
    # T_pmp::Vector{T}
    # Q_strn::Vector{T}
    # Q_b::T
    # Q_rock::T
    # bmb_grnd::T
    # H_ice::T
    # H_w::T
    # f_grnd::T
    # T_srf::T
    # T_shlf::T

    
    # # Diagnostic output variables
    # T::Vector{T}
    # ω::Vector{T}
    # Q_ice_b::T
    # H_cts::T
end

function IceColumn(T,nz::Int64)
    ζ   = zeros(T, nz)
    ζac = zeros(T, nz)
    dζa = zeros(T, nz)
    dζb = zeros(T, nz)
    return IceColumn(ζ,ζac,dζa,dζb)
end

# call calc_enth_column(enth(i,j,:),T_ice(i,j,:),omega(i,j,:),bmb_grnd(i,j),Q_ice_b(i,j), &
#                             H_cts(i,j),
#                             T_pmp(i,j,:),cp(i,j,:),kt(i,j,:),advecxy(i,j,:),
#                             uz(i,j,:),Q_strn(i,j,:), &
#                             Q_b(i,j),Q_rock(i,j),T_srf(i,j),T_shlf,H_ice_now,H_w(i,j),f_grnd(i,j),zeta_aa, &
#                             zeta_ac,dzeta_a,dzeta_b,cr,omega_max,T0,rho_ice,rho_w,L_ice,sec_year,dt)

# function enth_update!(icc::IceColumn, ict::IceColumn, icr::IceColumn,
#     ζ::Vector{T},T::AbstractFloat)


#     return
# end

end