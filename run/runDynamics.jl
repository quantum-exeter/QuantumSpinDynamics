using CSV
using DataFrames

include("../lib/Dynamics.jl")
using .Dynamics

### Parameters ###
prm = LorPrm3D(2.0, 0.001, 10.0, 2.0, 0.001, 10.0, 2.0, 0.001, 10.0) # Lorentzian parameters
ang =  CouplAng2D(π/2, 0.0, π/2, π/2) # Coupling angles
n = Lev2D(3, 3) # Number of RC levels
T = 0.01 # Temperature

### Time Range ###
ti, tf, dt = [7000 10000 100];
tspan = (ti, tf);
t = ti:dt:tf;

szDyn(ρ(1), n)

ρ = dsolve(prm, ang, n, T, tspan)
sz_list = [realIfClose(szDyn(ρ(i), n)) for i in t]

### Store Values ###
df = DataFrame(hcat(t, sz_list), :auto)
CSV.write("C://Users//crh222//Dropbox//PhD//1. RC Mapping//Dynamics Data//test.csv",  df, header = ["T", "sz"]);