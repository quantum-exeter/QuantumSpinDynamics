###########################
#### magnetisations.jl ####
###########################

sz(ρ, T) = tr(ρ*σz)

szGibbs(T) = sz(ρGibbs(T), T)
szMFGS(prm::Lorentzian, ang::CouplingAngles, n::Levels, T) = sz(ρMFGS(prm, ang, n, T), T)
szGround(prm::Lorentzian, ang::CouplingAngles, n::Levels, T) = sz(ρGround(prm, ang, n), T)
szAnalytical3D(T) = -tanh(1/T)