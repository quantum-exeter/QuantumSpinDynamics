using MKL
using CSV
using DataFrames
using ProgressMeter

include("../lib/Dynamics.jl")
using .Dynamics

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
# prm = LorPrm2D(2., 0.6, 1., 2., 0.6, 1.) 
# prm = LorPrm3D(2., 0.6, 1., 2., 0.6, 1., 2., 0.6, 1.)

## Coupling Angles ##
ang =  CouplAng1D(atan(1/sqrt(2)), π/4)
# ang =  CouplAng2D(π/2, 0.0, π/2, π/2)
# ang =  CouplAng3D(π/2, 0.0, π/2, π/2, 0.0, 0.0)

## RC Levels ##
n = Lev1D(10) # Number of RC levels
# n = Lev2D(5, 5) # Number of RC levels
# n = Lev3D(3, 3, 3) # Number of RC levels

## Temperature ##
T = 10

### Time Range ###
ti, tf, dt = [0 100 1];
tspan = (ti, tf);
t = ti:dt:tf;

sz_list = complex(zeros(length(t)))

ρ = dsolve(prm, ang, n, T, tspan)

sz_list = [realIfClose(szDyn(ρ(i), n)) for i in t]

### Store Values ###
df = DataFrame(hcat(t, sz_list), :auto)

### Export for Mac ###
CSV.write("/Users/charliehogg/1D_prmd_10_10.csv",  df, header = ["t", "sz"])

### Export for Windows ###
# CSV.write("C:/Users/crh222/filename.csv",  dfMFGS, header = ["t", "sz"])