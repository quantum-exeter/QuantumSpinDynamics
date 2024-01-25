using MKL
using NPZ
using DataFrames
using ProgressMeter

include("../lib/Statics.jl")
using .Statics

### Parameters ###

## Lorentzian parameters ##
prm1D = LorPrm1D(2., 0.6, 10.) 

## Coupling angles ##
ang3D =  CouplAng3D(π/2, 0.0, π/2, π/2, 0.0, 0.0)

## RC levels ##
n1D = Lev1D(10)
n3D = Lev3D(4, 4, 4)

θ = range(atan(sqrt(2)/2), atan((3*sqrt(2))/2), 50)

T = big(0.01)

ℰ_list1D = zeros(length(θ))
ℰ_list3D = zeros(length(θ))

@showprogress for i in eachindex(θ)
    ang1D = CouplAng1D(θ[i], π/4)
    λx = sqrt(3)*sin(θ[i])*cos(π/4)
    λy = sqrt(3)*sin(θ[i])*sin(π/4)
    λz = sqrt(3)*cos(θ[i])
    prm3D = LorPrm3D(2., 0.6, 10*λx^2, 2., 0.6, 10*λy^2, 2., 0.6, 10*λz^2)
    ℰ_list1D[i] = realIfClose(entropy(ComplexF64.(ρMFGS(prm1D, ang1D, n1D, T))))
    ℰ_list3D[i] = realIfClose(entropy(ComplexF64.(ρMFGS(prm3D, ang3D, n3D, T))))
end

### Export ###
npzwrite("anisotropy/anisotropy.npz",
    Dict("theta" => θ, 
         "entropy1D" => ℰ_list1D,
         "entropy3D" => ℰ_list3D))