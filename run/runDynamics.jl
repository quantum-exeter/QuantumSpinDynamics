using CSV
using DataFrames

include("../lib/Dynamics.jl")
using .Dynamics

prm = LorPrm1D(2.0, 0.001, 10.0)
ang =  CouplAng1D(0.0, 0.0)
n = Lev1D(2)
T = 10

### Time Range ###
ti, tf, dt = [0 100 2];
tspan = (ti, tf);
t = ti:dt:tf;

using LinearAlgebra
### z-Magnetisations ###
ρ = dsolve(prm, ang, n, T, tspan)
tr(ρ(0))

### Store Values ###
df = DataFrame(hcat(T, magz_list), :auto)
CSV.write("C://Users//crh222//Data//Classical//test.csv",  df, header = ["T", "sz"]);