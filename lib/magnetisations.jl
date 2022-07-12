###########################
#### magnetisations.jl ####
###########################

sz(ρ) = tr(ρ*σz)

### Statics ###
szGibbs(T) = sz(ρGibbs(T))
szMFGS(prm::Lorentzian, ang::CouplingAngles, n::Levels, T) = sz(ρMFGS(prm, ang, n, T))
szGround(prm::Lorentzian, ang::CouplingAngles, n::Levels, T) = sz(ρGround(prm, ang, n))
szAnalytical3D(T) = -tanh(1/T)

### Dynamics ###
szDyn(ρ, n) = tr(ρ*kronecker(σz, 𝕀(Int(hspace_size(n)/2))))