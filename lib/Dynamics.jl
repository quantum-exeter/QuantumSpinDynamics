module Dynamics

    ### Import Packages ###
    using LinearAlgebra
    using Kronecker
    using DifferentialEquations
    using SparseArrays

    ### Inclusions ###
    include("variables.jl")
    include("diffEqSolver.jl")
    include("hamiltonians.jl")
    include("magnetisations.jl")
    include("maths.jl")
    include("operators.jl")
    include("spectralDensity.jl")
    include("states.jl")
    include("superoperators.jl")
    include("transitions.jl")

    ####################################
    ####################################
    ####################################

    ### Exports ###
    export LorPrm1D, LorPrm2D, LorPrm3D,
           CouplAng1D, CouplAng2D, CouplAng3D,
           Lev1D, Lev2D, Lev3D,
           dsolve, szDyn, realIfClose

end