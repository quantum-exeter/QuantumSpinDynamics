module Statics

    ### Import Packages ###
    using LinearAlgebra
    using Kronecker

    ### Inclusions ###
    include("hamiltonians.jl")
    include("magnetisations.jl")
    include("maths.jl")
    include("operators.jl")
    include("states.jl")
    include("variables.jl")

    ####################################
    ####################################
    ####################################

    ### Exports ###
    export szGibbs, szMFGS, szGround, szAnalytical3D, realIfClose

end