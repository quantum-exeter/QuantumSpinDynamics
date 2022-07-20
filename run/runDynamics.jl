using CSV
using DataFrames

include("../lib/Dynamics.jl")
using .Dynamics

### Parameters ###
prm = LorPrm3D(2.0, 0.001, 10.0, 2.0, 0.001, 10.0, 2.0, 0.001, 10.0) # Lorentzian parameters
ang =  CouplAng3D(π/2, 0.0, π/2, π/2, 0.0, 0.0) # Coupling angles
n = Lev3D(3, 3, 3) # Number of RC levels
T = 10 # Temperature

### Time Range ###
ti, tf, dt = [8000 11000 100];
tspan = (ti, tf);
t = ti:dt:tf;

ρ = dsolve(prm, ang, n, T, tspan)
sz_list = [realIfClose(szDyn(ρ(i), n)) for i in t]

### Store Values ###
df = DataFrame(hcat(t, sz_list), :auto)
CSV.write("C://Users//crh222//Dropbox//PhD//1. RC Mapping//Dynamics Data//3d_prma_10_3.csv",  df, header = ["T", "sz"]);