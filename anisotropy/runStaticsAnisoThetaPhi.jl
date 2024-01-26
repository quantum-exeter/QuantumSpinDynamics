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
n1D = Lev1D(100)
n3D = Lev3D(5, 5, 5)

θ = range(atan(sqrt(2)/2), atan((3*sqrt(2))/2), 20)
ϕ = range(atan(0.5), atan(1.5), 20)

T = big(0.01)

ℰ_list1D = Float64[]
ℰ_list3D = Float64[]

for i in eachindex(θ)
    @showprogress for j in eachindex(ϕ)
        ang1D = CouplAng1D(θ[i], ϕ[j])
        λx = sqrt(3)*sin(θ[i])*cos(ϕ[j])
        λy = sqrt(3)*sin(θ[i])*sin(ϕ[j])
        λz = sqrt(3)*cos(θ[i])
        prm3D = LorPrm3D(2., 0.6, 10*λx^2, 2., 0.6, 10*λy^2, 2., 0.6, 10*λz^2)
        push!(ℰ_list1D, realIfClose(entropy(ComplexF64.(ρMFGS(prm1D, ang1D, n1D, T)))))
        push!(ℰ_list3D, realIfClose(entropy(ComplexF64.(ρMFGS(prm3D, ang3D, n3D, T)))))
    end
end

### Export ###
npzwrite("anisotropy/anisotropy_theta_phi.npz",
    Dict("theta" => θ,
         "phi" => ϕ,
         "entropy1D" => ℰ_list1D,
         "entropy3D" => ℰ_list3D)) 