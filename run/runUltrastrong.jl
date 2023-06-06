using MKL
using CSV
using DataFrames
using ProgressMeter

include("../lib/Statics.jl")
using .Statics

ang =  CouplAng2D(π/2, 0.0, 0.0, 0.0)
n = Lev2D(20, 20) # Number of RC levels

nTzero′(prm, ang, n)

α = (exp10.(range(-1, 5, length=100)))
szMFGS_list = zeros(length(α))
szMFGS_list_transformed = zeros(length(α))
occupancy_list = zeros(length(α))

@showprogress for i in eachindex(α)
    prm = LorPrm2D(2., 0.6, α[i], 2., 0.6, α[i])
    szMFGS_list[i] = real(szTzero(prm, ang, n))
    szMFGS_list_transformed[i] = real(szTzero′(prm, ang, n))
    occupancy_list[i] = real(nTzero′(prm, ang, n))
end

df = DataFrame(hcat(α, szMFGS_list, szMFGS_list_transformed), :auto)
CSV.write("other_data/sz_alpha_2D_sx_sz_20.csv",  df, header = ["α", "szMFGS", "szMFGS′"])