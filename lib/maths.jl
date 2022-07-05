##################
#### maths.jl ####
##################

### Commutators/Anticommutators ###
comm(A,B) = A*B - B*A
acomm(A,B) = A*B + B*A

### Square ###
square(n) = n*n

### Identity Matrices ###
𝕀(n) = Matrix(I, n, n) 
𝕀s = 𝕀(2) # Spin

### Partial Trace ###
function ptrace(ρ, n)
    nR = Int(size(ρ, 1)/n) # This is the remaining dimension
    lhs(i) = kronecker(𝕀(nR), (𝕀(n)[[i],:]))
    rhs(i) = kronecker(𝕀(nR), (𝕀(n)[:,i]))
    return sum(lhs(i)*ρ*rhs(i) for i=1:n)
end

### Uhlmann Fidelity ###
fidelity(ρ1, ρ2) = square(tr(sqrt(sqrt(ρ1)*ρ2*sqrt(ρ1))))

### Check for Choppable Components ###
realIfClose(c) = isnan(imag(c)) || imag(c) < 1e-14 ? real(c) : c;
realIfClose(c::AbstractArray) = realIfClose.(c);