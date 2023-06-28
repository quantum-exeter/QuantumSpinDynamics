###########################
#### magnetisations.jl ####
###########################

### Expectation values ###
sx(ρ) = tr(ρ*σx)  # expectation value sx in state ρ
sy(ρ) = tr(ρ*σy) # expectation value sy in state ρ
sz(ρ) = tr(ρ*σz) # expectation value sz in state ρ

### Statics ###
sxGibbs(T) = sx(ρGibbs(T)) # expectation value sx in Gibbs state at temperature T
sxMFGS(prm::Lorentzian, ang::CouplingAngles, n::Levels, T) = sx(ρMFGS(prm, ang, n, T)) # expectation value sx in MF state at temperature T
sxGround(prm::Lorentzian, ang::CouplingAngles, n::Levels, T) = sx(ρGround(prm, ang, n)) # expectation value sx in ground state at temperature T

syGibbs(T) = sy(ρGibbs(T)) # expectation value sy in Gibbs state at temperature T
syMFGS(prm::Lorentzian, ang::CouplingAngles, n::Levels, T) = sy(ρMFGS(prm, ang, n, T)) # expectation value sy in MF state at temperature T
syGround(prm::Lorentzian, ang::CouplingAngles, n::Levels, T) = sy(ρGround(prm, ang, n)) # expectation value sy in ground state at temperature T

szGibbs(T) = sz(ρGibbs(T)) # expectation value sz in Gibbs state at temperature T
szMFGS(prm::Lorentzian, ang::CouplingAngles, n::Levels, T) = sz(ρMFGS(prm, ang, n, T)) # expectation value sz in MF state at temperature T
szGround(prm::Lorentzian, ang::CouplingAngles, n::Levels, T) = sz(ρGround(prm, ang, n)) # expectation value sz in ground state at temperature T

### Dynamics ###
sxDyn(ρ, n) = tr(ρ*kronecker(σx, 𝕀(Int(hspace_size(n)/2)))) # expectation value sx in total dynamical state ρ, tracing out space of dimension n
syDyn(ρ, n) = tr(ρ*kronecker(σy, 𝕀(Int(hspace_size(n)/2)))) # expectation value sy in total dynamical state ρ, tracing out space of dimension n
szDyn(ρ, n) = tr(ρ*kronecker(σz, 𝕀(Int(hspace_size(n)/2)))) # expectation value sz in total dynamical state ρ, tracing out space of dimension n