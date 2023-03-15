using MKL
using CSV
using DataFrames
using ProgressMeter

include("../lib/Statics.jl")
using .Statics

### Parameters ###

## Lorentzian Spectral Density ##
# {ω_0, Γ, α}
#prma = 2., 0.001, 10.
#prmb = 2., 0.001, 1. 
#prmc = 2., 0.001, 0.1 

#prmd = 2., 0.6, 10.
#prme = 2., 0.6, 1.
#prmf = 2., 0.6, 0.1

#prmg = 2., 0.6, 20.
#prmh = 2., 0.6, 60.
#prmi = 2., 0.6, 80.
#prmj = 2., 0.6, 500.
#prmk = 2., 0.6, 1000.

prm = LorPrm1D(2., 0.6, 10.) 
# prm = LorPrm2D(2., 0.6, 1000., 2., 0.6, 1000.) 
# prm = LorPrm3D(2., 0.6, 0.1, 2., 0.6, 0.1, 2., 0.6, 0.1)

## Coupling Angles ##
ang =  CouplAng1D(atan(sqrt(2)), π/4)
# ang =  CouplAng2D(π/2, 0.0, π/2, π/2)
# ang =  CouplAng3D(π/2, 0.0, π/2, π/2, 0.0, 0.0)

## RC Levels ##
n = Lev1D(100) # Number of RC levels
# n = Lev2D(10, 10) # Number of RC levels
# n = Lev3D(6, 6, 6) # Number of RC levels

## Temperature Range ##
T = exp10.(range(-3, 2, length=100))

sxGS_list = zeros(length(T))
syGS_list = zeros(length(T))
szGS_list = zeros(length(T))

@showprogress for i in eachindex(T)
    sxGS_list[i] = real(sxGibbs(big(T[i])))
    syGS_list[i] = real(syGibbs(big(T[i])))
    szGS_list[i] = real(szGibbs(big(T[i])))
end

sxMFGS_list = zeros(length(T))
syMFGS_list = zeros(length(T))
szMFGS_list = zeros(length(T))

@showprogress for i in eachindex(T)
    if T[i] < 0.01
        Ti = big(T[i])
    else
        Ti = T[i]
    end
    sxMFGS_list[i] = real(sxMFGS(prm, ang, n, Ti))
    syMFGS_list[i] = real(syMFGS(prm, ang, n, Ti))
    szMFGS_list[i] = real(szMFGS(prm, ang, n, Ti))
end

### Store Values ###
dfGibbs = DataFrame(hcat(T, sxGS_list, syGS_list, szGS_list), :auto)
dfMFGS = DataFrame(hcat(T, sxMFGS_list, syMFGS_list, szMFGS_list), :auto)

### Export for Mac ###
CSV.write("paper_data/qu_Gibbs.csv",  dfGibbs, header = ["T", "sxGS", "syGS", "szGS"])
CSV.write("paper_data/qu_MFGS_1D_prmd_100.csv",  dfMFGS, header = ["T", "sxMFGS", "syMFGS", "szMFGS"])

### Export for Windows ###
# CSV.write(".//paper_data/qu_Gibbs.csv",  dfGibbs, header = ["T", "sxGS", "syGS", "szGS"])
# CSV.write(".//paper_data/qu_MFGS_1D_prmd_200.csv",  dfGibbs, header = ["T", "sxMFGS", "syMFGS", "szMFGS"])