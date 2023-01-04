module Variational

    ### Import Packages ###
    using LinearAlgebra
    using Kronecker
    using QuadGK
    using ForwardDiff

    ### Inclusions ###
    include("variables.jl")
    include("hamiltonians.jl")
    include("magnetisations.jl")
    include("maths.jl")
    include("operators.jl")
    include("SDintegrals.jl")
    include("spectralDensity.jl")
    include("states.jl")
    include("wkCoupling.jl")

    ####################################
    ####################################
    ####################################

    ### Exports ###
    export LorPrm1D, LorPrm2D, LorPrm3D,
           CouplAng1D, CouplAng2D, CouplAng3D,
           Lev1D, Lev2D, Lev3D,
           sxGibbs, sxMFGS, sxGround, 
           syGibbs, syMFGS, syGround, 
           szGibbs, szMFGS, szGround, szAnalytical3D, 
           szWK, szWK2, szWKZT,
           realIfClose, ρMFGS,
           ptrace, ψGround, ϕp, ϕm
           
end