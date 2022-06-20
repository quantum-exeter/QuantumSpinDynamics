module Statics

    ### Import Packages ###
    using LinearAlgebra
    using Kronecker

    ### Inclusions ###
    include("constants.jl")
    include("variables.jl")
    include("maths.jl")
    include("initialStates.jl")
    include("hamiltonians.jl")

    ####################################
    ####################################
    ####################################

    ### Exports ###
    export ρG, ρMFGS, sG, sMFGS, groundState, szAnalytical, realIfClose, T
    
    ### Statics-Specific Functions ###

    ## Gibbs State ##
    ρG(T) = gibbs(HG(), T)

    ## Mean-Force Gibbs State (MFGS) ##
    function ρMFGS(dim, T)
        H = HS(dim) # Hamiltonian for the given coupling dimension
        proj = eigen(H).vectors # Projector onto the basis of H
        HTr = adjoint(proj)*H*proj # Transformed H
        stateTot = proj*gibbs(HTr, T)*adjoint(proj)
        n = Int(hspace_dimension(dim)/2)
        return ptrace(stateTot, n)
    end

    # Spin Expectations #
    exps(ρ) = [tr(ρ*σx) tr(ρ*σy) tr(ρ*σz)]
    sG(T) = exps(ρG(T))
    sMFGS(dim, T) = exps(ρMFGS(dim, T))

    ## Ground State ##
    function groundState(dim)
        H = HS(dim)
        n = Int(hspace_dimension(dim)/2)
        state = eigen(H).vectors[:,1]
        sx = adjoint(state)*kronecker(σx, 𝕀(n))*state
        sy = adjoint(state)*kronecker(σy, 𝕀(n))*state
        sz = adjoint(state)*kronecker(σz, 𝕀(n))*state
        return [sx sy sz]
    end

    ## Sz Coupling Analytical Expression ##

    szAnalytical(T) = -tanh(1/T)

end