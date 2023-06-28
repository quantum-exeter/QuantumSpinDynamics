using MKL
using CSV
using DataFrames
using ProgressMeter

include("../lib/Dynamics.jl")
using .Dynamics

### Parameters ###

## Lorentzian spectral density ##
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
# prm = LorPrm2D(2., 0.6, 10., 2., 0.6, 10.) 
# prm = LorPrm3D(2., 0.6, 10., 2., 0.6, 10., 2., 0.6, 10.)

## Coupling angles ##
ang =  CouplAng1D(atan(sqrt(2)), π/4)
# ang =  CouplAng2D(π/2, 0.0, π/2, π/2)
# ang =  CouplAng3D(π/2, 0.0, π/2, π/2, 0.0, 0.0)

## RC levels ##
n = Lev1D(10)
# n = Lev2D(5, 5)
# n = Lev3D(3, 3, 3)

## Temperature ##
T = 0.01

### Time range ###
ti, tf, dt = [0 100 1];
tspan = (ti, tf);
t = ti:dt:tf;

### Dynamic state expectation values ###

## Initialise lists ##
sz_list = complex(zeros(length(t)))
ρ = dsolve(prm, ang, n, T, tspan)
sz_list = [realIfClose(szDyn(ρ(i), n)) for i in t]

### Store values ###
df = DataFrame(hcat(t, sz_list), :auto)

### Export ###
CSV.write("filename.csv",  dfGibbs, header = ["T", "sz"])