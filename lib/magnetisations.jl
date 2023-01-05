###########################
#### magnetisations.jl ####
###########################

### Expectation Values ###
sx(ρ) = tr(ρ*σx)
sy(ρ) = tr(ρ*σy)
sz(ρ) = tr(ρ*σz)

### Statics ###
sxGibbs(T) = sx(ρGibbs(T))
sxMFGS(prm::Lorentzian, ang::CouplingAngles, n::Levels, T) = sx(ρMFGS(prm, ang, n, T))
sxGround(prm::Lorentzian, ang::CouplingAngles, n::Levels, T) = sx(ρGround(prm, ang, n))

syGibbs(T) = sy(ρGibbs(T))
syMFGS(prm::Lorentzian, ang::CouplingAngles, n::Levels, T) = sy(ρMFGS(prm, ang, n, T))
syGround(prm::Lorentzian, ang::CouplingAngles, n::Levels, T) = sy(ρGround(prm, ang, n))

szGibbs(T) = sz(ρGibbs(T))
szMFGS(prm::Lorentzian, ang::CouplingAngles, n::Levels, T) = sz(ρMFGS(prm, ang, n, T))
szGround(prm::Lorentzian, ang::CouplingAngles, n::Levels, T) = sz(ρGround(prm, ang, n))
szAnalytical3D(T) = -tanh(1/T)

### Dynamics ###
sxDyn(ρ, n) = tr(ρ*kronecker(σx, 𝕀(Int(hspace_size(n)/2))))
syDyn(ρ, n) = tr(ρ*kronecker(σt, 𝕀(Int(hspace_size(n)/2))))
szDyn(ρ, n) = tr(ρ*kronecker(σz, 𝕀(Int(hspace_size(n)/2))))