######################
#### variables.jl ####
######################

γ = 1 # Gyromagnetic ratio for an electron
Λ = 10^10 # Cutoff frequency for spectral density
s0 = 1/2 # Spin size

### Lorentzian Parameters ###
abstract type Lorentzian end

struct LorPrm1D{T <: Real} <: Lorentzian
    ω01::T
    Γ1::T
    α1::T
end

struct LorPrm2D{T <: Real} <: Lorentzian
    ω01::T
    Γ1::T
    α1::T
    ω02::T
    Γ2::T
    α2::T
end

struct LorPrm3D{T <: Real} <: Lorentzian
    ω01::T
    Γ1::T
    α1::T
    ω02::T
    Γ2::T
    α2::T
    ω03::T
    Γ3::T
    α3::T
end

### Truncation of Harmonic Oscillator ###
abstract type Levels end

struct Lev1D <: Levels
    n1::Int
end

struct Lev2D <: Levels
    n1::Int
    n2::Int
end

struct Lev3D <: Levels
    n1::Int
    n2::Int
    n3::Int
end

### Coupling Angles ###
abstract type CouplingAngles end

struct CouplAng1D{T<:Real} <: CouplingAngles
    θ1::T
    ϕ1::T
end

struct CouplAng2D{T<:Real} <: CouplingAngles
    θ1::T
    ϕ1::T
    θ2::T
    ϕ2::T
end

struct CouplAng3D{T<:Real} <: CouplingAngles
    θ1::T
    ϕ1::T
    θ2::T
    ϕ2::T
    θ3::T
    ϕ3::T
end

hspace_size(n::Lev1D) = 2*n.n1
hspace_size(n::Lev2D) = 2*n.n1*n.n2
hspace_size(n::Lev3D) = 2*n.n1*n.n2*n.n3

dim(n::Lev1D) = 1
dim(n::Lev2D) = 2
dim(n::Lev3D) = 3