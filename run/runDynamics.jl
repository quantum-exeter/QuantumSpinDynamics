using CSV
using DataFrames

include("../lib/Dynamics.jl")
using .Dynamics

### Parameters ###
prm = LorPrm1D(2.0, 0.001, 10.0) # Lorentzian parameters
ang =  CouplAng1D(0.0, 0.0) # Coupling angles
n = Lev1D(2) # Number of RC levels
T = 10 # Temperature

### Time Range ###
ti, tf, dt = [0 100 2];
tspan = (ti, tf);
t = ti:dt:tf;

sz_list = [szDyn(prm, ang, n, T, tspan, i) for i in tspan]

### Store Values ###
df = DataFrame(hcat(T, sz_list), :auto)
CSV.write("C://Users//crh222//Data//Classical//test.csv",  df, header = ["T", "sz"]);