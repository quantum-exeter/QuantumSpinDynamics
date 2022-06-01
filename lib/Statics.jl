module Statics

    ### Import Packages ###
    using LinearAlgebra
    using Kronecker

    ### Inclusions ###
    include("variables.jl")
    include("maths.jl")
    include("initialStates.jl")
    include("hamiltonians.jl")
    include("expectationValues.jl")

    ####################################
    ####################################
    ####################################

    ### Exports ###
    export sG, s1D, s2D, s3D,
           gs1D, gs2D, gs3D
    
    ### Statics-Specific Functions ###

    function MFGS(H, T, n)
        proj = eigen(H).vectors
        hT = adjoint(proj)*H*proj
        stateTot = proj*gibbs(hT, T)*adjoint(proj)
        return ptrace(stateTot, n)
    end

    ρG(T) = gibbs(HG(), T)
    ρMFGS1D(T) = MFGS(HRC1D(n1, λ1, Ω1), T, n1)
    ρMFGS2D(T) = MFGS(HRC2D(n1, n2, λ1, λ2, Ω1, Ω2), T, n1*n2)
    ρMFGS3D(T) = MFGS(HRC3D(n1, n2, n3, λ1, λ2, λ3, Ω1, Ω2, Ω3), T, n1*n2*n3)

    sG(T) = exps(ρG(T))
    s1D(T) = exps(ρMFGS1D(T))
    s2D(T) = exps(ρMFGS2D(T))
    s3D(T) = exps(ρMFGS3D(T))

    ## Ground State ##
    function ground_state(H, n)
        state = eigen(H).vectors[:,1]
        sx = -adjoint(state)*kronecker(sx0, 𝕀(n))*state
        sy = -adjoint(state)*kronecker(sy0, 𝕀(n))*state
        sz = -adjoint(state)*kronecker(sz0, 𝕀(n))*state
        return [sx sy sz]
    end

    gs1D = ground_state(HRC1D(n1, λ1, Ω1), n1)
    gs2D = ground_state(HRC2D(n1, n2, λ1, λ2, Ω1, Ω2), n1*n2)
    gs3D = ground_state(HRC3D(n1, n2, n3, λ1, λ2, λ3, Ω1, Ω2, Ω3), n1*n2*n3)

end