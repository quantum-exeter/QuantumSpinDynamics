using CSV
using DataFrames

include("../lib/Dynamics.jl")
using .Dynamics

### Low Gamma ###
#prma = 2.0, 0.001, 10.0
#prmb = 2.0, 0.001, 1.0

### Weak ###
#prmc = 2.0, 0.001, 0.1

### High Gamma ###
#prmd = 2.0, 0.6, 10.0
#prme = 2.0, 0.6, 1.0

### Parameters ###
prm = LorPrm1D(2., 0.6, 1.) 
# prm = LorPrm2D(2., 0.001, 10., 2., 0.001, 10.) 
# prm = LorPrm3D(2., 0.6, 1., 2., 0.6, 1., 2., 0.6, 1.) # Lorentzian parameters
ang =  CouplAng1D(π/4, 0.0) # Coupling angles
# ang =  CouplAng2D(π/2, 0.0, π/2, π/2) # Coupling angles
# ang =  CouplAng3D(π/2, 0.0, π/2, π/2, 0.0, 0.0) # Coupling angles
n = Lev1D(10) # Number of RC levels
# n = Lev2D(10, 10) # Number of RC levels
# n = Lev3D(3, 3, 3) # Number of RC levels

T = 10 # Temperature

### Time Range ###
ti, tf, dt = [0 50 1];
tspan = (ti, tf);
t = ti:dt:tf;

ρ = dsolve(prm, ang, n, T, tspan)
sz_list = [realIfClose(szDyn(ρ(i), n)) for i in t]

### Store Values ###
df = DataFrame(hcat(t, sz_list), :auto)
CSV.write("C://Users//crh222//Dropbox//PhD//1. 3D Project//Dynamics Data//1D_pi4_prmh_10_10.csv",  df, header = ["T", "sz"]);