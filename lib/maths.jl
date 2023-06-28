##################
#### maths.jl ####
##################

### Commutators/anticommutators ###
comm(A, B) = A*B - B*A # commutator of A and B
acomm(A, B) = A*B + B*A # anticommutator of A and B

### Identity matrices ###
𝕀(n) = Matrix(I, n, n) # identity matrix of dimension n
𝕀s = 𝕀(2) # identity matrix of dimension 2

### Partial trace (tracing out bath) ###
function ptrace(ρ, n)
    nR = Int(size(ρ, 1)/n) # remaining dimension
    lhs(i) = kronecker(𝕀(nR), (𝕀(n)[[i],:]))
    rhs(i) = kronecker(𝕀(nR), (𝕀(n)[:,i]))
    return sum(lhs(i)*ρ*rhs(i) for i=1:n)
end

### Partial trace (tracing out spin) ###
function ptraceSp(ρ)
  n = 2
  nR = Int(size(ρ, 1)/n) # remaining dimension
  lhs(i) = kronecker((𝕀(n)[[i],:]), 𝕀(nR))
  rhs(i) = kronecker((𝕀(n)[:,i]), 𝕀(nR))
  return sum(lhs(i)*ρ*rhs(i) for i=1:n)
end

### Check for choppable components ###
realIfClose(c) = isnan(imag(c)) || imag(c) < 1e-14 ? real(c) : c; # returns the real part of c if the imaginary part is small or NaN
realIfClose(c::AbstractArray) = realIfClose.(c); # returns the real part of each element of c if their imaginary part is small or NaN
chopReal(x) = real(x) < 1e-12 ? imag(x)*1im : x; # chop the real components
chopImag(x) = imag(x) < 1e-12 ? real(x) : x; # chop the imaginary components
chopBoth(x) = chopReal(chopImag(x)); # chop both the real and imaginary components

### Integration ###
xcoth(x) = iszero(x) ? one(x) : x*coth(x)

## Double integration ##
function dblquadgk(f, a::AbstractArray{T}, b::AbstractArray{T};
    rtol=sqrt(eps(T)), atol=zero(T), maxevals=10^7, order=7) where T<:AbstractFloat
J(x) = quadgk(y -> f(x,y), a[2], b[2], atol=atol, maxevals=maxevals, order=order)[1]
K = quadgk(x -> J(x), a[1], b[1], atol=atol, maxevals=maxevals, order=order)[1]
return K
end

## Principal value integrals with quadgk ##
function quadgk_cauchy(f, a, c, b)
  fc = f(c)
  g(x) = (f(x) - fc)/(x - c)
  I = quadgk(g, a, c, b, order=12)
  C = fc*log(abs((b - c)/(a - c))) 
  return (I[1] + C, I[2], C)
end
  
function quadgk_hadamard(f, a, c, b)
  df(x) = ForwardDiff.derivative(f,x)
  I = quadgk_cauchy(df, a, c, b)
  C = f(a)/(a-c) - f(b)/(b-c)
  return (C + I[1], I[2], C)
end