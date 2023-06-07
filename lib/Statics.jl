module Statics

    ### Import Packages ###
    using LinearAlgebra
    using Kronecker
    using QuadGK
    using ForwardDiff
    using QuantumOptics

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
    include("entanglement.jl")
    include("ultrastrong.jl")

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
           szWK, szWKZT,
           realIfClose, ρMFGS,
           entropy,
           szTzero, szTzero′,
           nTzero, nTzero′,
           osc_positions

end