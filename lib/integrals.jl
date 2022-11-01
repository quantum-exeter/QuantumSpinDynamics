######################
#### integrals.jl ####
######################

### Import Packages ###
using LinearAlgebra
using Kronecker
using QuadGK
using ForwardDiff

### Inclusions ###
include("variables.jl")
include("hamiltonians.jl")
include("magnetisations.jl")
include("maths.jl")
include("operators.jl")
include("states.jl")
include("spectralDensity.jl")


xcoth(x) = iszero(x) ? one(x) : x*coth(x)

function dblquadgk(f, a::AbstractArray{T}, b::AbstractArray{T};
    rtol=sqrt(eps(T)), atol=zero(T), maxevals=10^7, order=7) where T<:AbstractFloat
J(x) = quadgk(y -> f(x,y), a[2], b[2], atol=atol, maxevals=maxevals, order=order)[1]
K = quadgk(x -> J(x), a[1], b[1], atol=atol, maxevals=maxevals, order=order)[1]
return K
end

function quadgk_cauchy(f, a, c, b)
    fc = f(c)
    g(x) = (f(x) - fc)/(x - c)
    I = quadgk(g, a, c, b)
    C = fc*log(abs((b - c)/(a - c))) 
    return (I[1] + C, I[2], C)
  end
  
  function quadgk_hadamard(f, a, c, b)
    df(x) = ForwardDiff.derivative(f,x)
    I = quadgk_cauchy(df, a, c, b)
    C = f(a)/(a-c) - f(b)/(b-c)
    return (C + I[1], I[2], C)
  end

function Σ(prm::LorPrm1D)
    I(ω) = spectral_density(ω, prm)*ω/(ω + 1)
    Il = quadgk_cauchy(I, 0.0, 1.0, 2.0)[1]
    Ir = quadgk(ω -> I(ω)/(ω - 1), 2.0, Inf)[1]
    return Il + Ir
  end
  
  function Σ′(prm::LorPrm1D)
    I(ω) = spectral_density(ω, prm)*ω/(ω + 1)^2
    Il = quadgk_hadamard(I, 0.0, 1.0, 2.0)[1]
    Ir = quadgk(ω -> I(ω)/(ω - 1)^2, 2.0, Inf)[1]
    return 2*(Il + Ir)
  end
  
  function Δ(prm::LorPrm1D, β)
    I(ω) = spectral_density_inv(ω, prm)*xcoth(β*ω)/(β)/(ω + 1)
    Il = quadgk_cauchy(I, 0.0, 1.0, 2.0)[1]
    Ir = quadgk(ω -> I(ω)/(ω - 1), 2.0, Inf)[1]
    return Il + Ir
  end
  
  function Δ′(prm::LorPrm1D, β)
    I(ω) = 0.5*spectral_density_inv(ω, prm)*(ω^2 + 1)*xcoth(β*ω)/(β)/(ω + 1)^2
    Il = quadgk_hadamard(I, 0.0, 1.0, 2.0)[1]
    Ir = quadgk(ω -> I(ω)/(ω - 1)^2, 2.0, Inf)[1]
    return Il + Ir
  end

  prm = LorPrm1D(7.0, 5.0, 1.0)
  
Δ′(prm, 1.0)