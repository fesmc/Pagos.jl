#=
# Basal friction

Basal friction laws describe the relationship between basal shear stress and basal sliding velocity. Based on laboratory experiments, the most common model to describe this relationship is the (regularized) Coulomb friction law. In Pagos.jl, this is implemented as `RegularizedCoulombBasalFriction`:
=#
using Pagos, CairoMakie
v_basal = [100, 0]
c_basal = 1
friction_coulomb = RegularizedCoulombBasalFriction()
τ_basal_coulomb = basal_shear_stress(v_basal, c_basal, friction_coulomb)

#=
Another commonly used basal friction law is the pseudo-plastic power law, implemented as `PseudoPlasticPowerBasalFriction`. We propose to compare both friction laws in the following figure.
=#
v_basal_norm = range(0, stop = 200, step = 0.1)
v_basal = [[v, 0] for v in v_basal_norm]
c_basal = ones(length(v_basal))
τ_basal_coulomb = basal_shear_stress(v_basal, c_basal, friction_coulomb)

friction_powerlaw = PseudoPlasticPowerBasalFriction()
τ_basal_powerlaw = basal_shear_stress(v_basal, c_basal, friction_powerlaw)

fig = plot_basal_shear_stress(v_basal, [τ_basal_coulomb, τ_basal_powerlaw], ["Coulomb", "Power law"])