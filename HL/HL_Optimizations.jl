################################################################################
############# OPTIMIZATIONS.JL - solve value function optimization #############
################################################################################

using Distributions, FastGaussQuadrature, Interpolations, Optim

################################################################################

# Working life with habits
function bellOpt{T<:AbstractFloat}(x::T, h::T, a::T, b::T, z::T, wmin::T,
  v_int::Interpolations.GriddedInterpolation, yn::Normal, k::Array, ρ::T,
  r::T, δ::T, λ::T)

  function EVprime(w′::Float64, h=h, a=a, b=b, z=z, yn=yn, k=k, v_int=v_int,ρ=ρ,
    λ=λ, x=x)

    h′ = (1-λ)*h + λ*(x - w′)

    function EVp(y::Float64, w′=w′, h′=h′, v_int=v_int, yn=yn, a=a, b=b, z=z,
                                                                       ρ=ρ, k=k)
      dy = y - mean(yn)
      v_int[w′+exp(y), h′, a+k[1]*dy, b+k[2]*dy, ρ*z+k[3]*dy]
    end

    (n, wgt) = gausshermite(15)
    π^(-0.5)*sum( [wgt[i]*EVp(sqrt(2)*std(yn)*n[i] + mean(yn))
                                                          for i = 1:length(n)] )
  end

  Blmn(w′::Float64, x=x, r=r, δ=δ) = -( u(x - w′, h) + δ*EVprime(r*w′) )

  optimum = optimize(Blmn, wmin/r, x + abs(wmin/r) + 0.001)
  w′ = optimum.minimum
  vopt = -(optimum.f_minimum)

  return w′, vopt
end

################################################################################

# Transition period with habits
function bellOpt_TRANS{T<:AbstractFloat}(x::T, h::T, pension::T, wmin::T,
  v_int::Interpolations.GriddedInterpolation, r::T, δ::T, λ::T)

  Blmn(w′::Float64) = -(u(x-w′, h) + δ*v_int[r*w′, (1-λ)*h+λ*(x-w′), pension])

  optimum = optimize(Blmn, wmin/r, x + abs(wmin/r) + 0.001)
  w′ = optimum.minimum
  vopt = -(optimum.f_minimum)

  return w′, vopt
end

################################################################################

# Retirement period with habits
function bellOpt_R{T<:AbstractFloat}(w::T, h::T, y::T, wmin::T,
  v_int::Interpolations.GriddedInterpolation, r::T, δ::T, λ::T)

  x = w + y

  Blmn(w′) = -( u(x - w′, h) + δ*v_int[r*w′, (1-λ)*h + λ*(x-w′), y] )

  optimum = optimize(Blmn, wmin/r, x + abs(wmin/r) + 0.001)
  w′ = optimum.minimum
  vopt = -(optimum.f_minimum)

  return w′, vopt
end
