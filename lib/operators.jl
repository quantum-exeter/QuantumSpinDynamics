######################
#### operators.jl ####
######################

### Pauli matrices ###
σx = [[0 1];[1 0]]
σy = [[0 -im];[im 0]]
σz = [[1 0];[0 -1]]

### Spin coupling operators ###
sc(θ, ϕ) = σx*(sin(θ)*cos(ϕ)) + σy*(sin(θ)*sin(ϕ)) + σz*cos(θ)