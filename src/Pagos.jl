"""
    Pagos

An ice-sheet model written in pure Julia, designed to be accessible, modular and performant.

Pagos provides a hierarchy of momentum balances (from the shallow-ice approximation up to
full Stokes, see [`AbstractMomentumBalance`](@ref)) together with the material laws, basal
friction models and time integrators needed to advance an [`IceSheet`](@ref). Kernels are
written with `KernelAbstractions`, so the same code runs on CPU and GPU backends.

!!! warning
    Pagos is work in progress: it is not yet fully functional and its API is subject to
    major changes.
"""
module Pagos

using Adapt
using DocStringExtensions
using Downloads
using KernelAbstractions
using LinearAlgebra
using LinearSolve           # TODO: externalize
using LoopVectorization
using Random
using StatsBase
using SparseArrays
using Tullio

###########################################################
# Structs
###########################################################

include("api/domain.jl")
export AbstractGrid, RegularGrid, CommonGrid

include("api/state.jl")
export TopographicState, MechanicState, ThermodynamicState, MaterialState

include("api/constants.jl")
export Constants

include("api/icesheet.jl")
export IceSheet, Topography, Tracers, Mechanics, Thermodynamics, Material

include("api/simulation.jl")
export Simulation

include("api/step.jl")

###########################################################
# Utils
###########################################################

include("utils/indices.jl")
export AbstractIndexing, StrictIndexing, FlatIndexing
export ReflectiveIndexing, PeriodicIndexing
export index, stencil_fd, stencil

include("utils/math.jl")
include("utils/debug.jl")

include("utils/integrators.jl")
export AbstractIntegrationMethod, Euler, RungeKutta4, BogackiShampine32, Tsitouras54, RKL2
export SSPRK33, SSPRK43, ssp_coefficient
export Integrator, step!, substep!, estimate_spectral_radius!

include("utils/oceananigans.jl")
export StaggeredGrids

###########################################################
# Topography
###########################################################

include("topography/sigmatransform.jl")
export AbstractSigmaTransform, PowerSigmaTransform
export LinearSigmaTransform, QuadraticSigmaTransform
export VerticalLayering, CorrectedVerticalLayering
export ζ_aa, ζ_ac, sigma
export get_ζ_aa, get_ζ_ac

include("topography/calving.jl")
export AbstractCalving, AbstractGroundedCalving, AbstractFloatingCalving
export PrescribedCalving, RelaxedCalving, ThicknessCalving
export LipscombCalving, LevermannCalving, CrawfordCalving, BassisCalving
export EigenCalving, FlotationCalving, PollardDeContoCalving
export calving_rate, calving_rate!

# TODO: update to new API (State → TopographicState/MechanicState, Domain → RegularGrid)
# include("topography/advection.jl")

###########################################################
# Dynamics
###########################################################

include("mechanics/momentum.jl")
export AbstractMomentumBalance
export SIAMomentumBalance, SSAMomentumBalance, SIASSAMomentumBalance
export DIVAMomentumBalance, BlatterPattynMomentumBalance, StokesMomentumBalance
export InertialDIVAMomentumBalance, InertialSIASSAMomentumBalance

include("mechanics/solvers.jl")
export AbstractMomentumSolver, LinearMomentumSolver2D
export populate_vectors!, velocity!

include("legacy/performance/velocities.jl")
export LegacyLinearMomentumSolver2D

include("mechanics/effective_pressure.jl")
export AbstractEffectivePressure, PrescribedEffectivePressure
export OverburdenEffectivePressure
export LeguyEffectivePressure, TillEffectivePressure
export effective_pressure, effective_pressure!

include("mechanics/basal_friction.jl")
export AbstractFriction, BasalFriction
export AbstractBasalBeta, PrescribedBasalBeta
export PseudoPlasticPowerBasalBeta, CoulombBasalBeta
export basal_shear_stress, basal_shear_stress!

include("mechanics/velocities3D.jl")
export aggregate_viscosity_integral!, aggregated_viscosity_integral!
export velocities3D!, surface_velocity!, depthaveraged_velocity!
export basal_velocity_from_surface_velocity!, depthavg_velocity!

include("mechanics/pseudotransient.jl")
export PseudoTransientSolver
export pseudo_transient!, pseudo_rate!, pseudo_dt, pseudo_vel!, dotvel!

include("mechanics/inertial.jl")
export inertial_velocity!

# include("mechanics/friction/plastic.jl")
# include("mechanics/friction/stagger.jl")

include("mechanics/stress.jl")
export shearstress!, basalstress!, drivingstress!
export deviatoric_stress!

include("mechanics/strainrate.jl")
export strainrate!, scaledstrainrate!, velocitygradients!
export raw_strainrate!, raw_strainrate_effective!, FullColumnMomentumBalance

###########################################################
# Material
###########################################################

include("material/anisotropy.jl")
export AbstractAnisotropy, EnhancementAnisotropy, CAFFEAnisotropy
export anisotropy!, enhancement_factor
export deformability, square_tangential_invariant

include("material/rate_factor.jl")
export AbstractRateFactor, PrescribedRateFactor, ArrheniusRateFactor
export SmithMorlandRateFactor, HookeRateFactor, LliboutryDuvalRateFactor
export FanLowStrainGSIRateFactor, FanLowStrainGSS1RateFactor, FanLowStrainGSS2RateFactor
export FanHighStrainGSIRateFactor
export rate_factor, rate_factor!

include("material/creep.jl")
export AbstractCreep, GlenNyeCreep, SmithMorlandCreep
export PettitWaddingtonCreep, GoldsbyKohlstedtCreep
export FanLowStrainCreep
export creep, creep!

include("material/flow_law.jl")
export AbstractFlowLaw, PrescribedViscosityFlowLaw, RateCreepFlowLaw
export GlenNyeFlowLaw, SmithMorlandFlowLaw
export FanLowStrainFlowLaw, FanHighStrainFlowLaw
export viscosity, viscosity!

include("material/pressure_melting_point.jl")
export AbstractPressureMeltingPoint, LinearPressureMeltingPoint
export PrescribedPressureMeltingPoint, LinearSalinityPressureMeltingPoint
export pressure_melting_point, pressure_melting_point!
export relative_temperature, relative_temperature!
export thermal_forcing, thermal_forcing!

###########################################################
# Numerics
###########################################################

include("numerics/differences.jl")
include("numerics/picard.jl")
export ∂x₁, ∂x₂, ∂x₃, ∂x!, ∂y!, ∂x₃!
export ∂x₁₂, ∂x₁₂!

###########################################################
# Extensions
###########################################################

include("ext/plots.jl")

end