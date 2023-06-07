using MKL
using CSV
using DataFrames
using ProgressMeter
using Plots
using LinearAlgebra

include("../lib/Statics.jl")
using .Statics

prm = LorPrm2D(2., 0.6, 100., 2., 0.6, 100.)
ang =  CouplAng2D(π/2, 0.0, 0.0, 0.0)
n = Lev2D(20, 20) # Number of RC levels

x = osc_positions(prm, ang, n)[1]

vline!([-x_max, x_max])

α = (exp10.(range(-1, 5, length=100)))
szMFGS_list = zeros(length(α))
szMFGS_list_transformed = zeros(length(α))
occupancy_list = zeros(length(α))
occupancy_list_transformed = zeros(length(α))
occupancy_list_transformed_sq = zeros(length(α))

@showprogress for i in eachindex(α)
    prm = LorPrm2D(2., 0.6, α[i], 2., 0.6, α[i])
    szMFGS_list[i] = real(szTzero(prm, ang, n))
    szMFGS_list_transformed[i] = real(szTzero′(prm, ang, n))
    occupancy_list[i] = real(nTzero(prm, ang, n))
    occupancy_list_transformed[i] = real(nTzero′(prm, ang, n))
end

df = DataFrame(hcat(α, szMFGS_list, szMFGS_list_transformed, occupancy_list, occupancy_list_transformed), :auto)
CSV.write("other_data/sz_alpha_2D_sx_sz_30.csv",  df, header = ["α", "szMFGS", "szMFGS′", "occupancy", "occupancy′"])