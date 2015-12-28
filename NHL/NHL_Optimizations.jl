################################################################################
########## NHL_OPTIMIZATIONS.JL - solve value function optimization ############
################################################################################

using Distributions, FastGaussQuadrature, Interpolations, Optim

################################################################################
# Working life, no habits
function bellOpt{T<:AbstractFloat}(x::T, a::T, b::T, z::T, wmin::T,
   v_int::Interpolations.GriddedInterpolation, yn::Normal, k::Array, ρ::T,
   r::T, δ::T, σ::T)

  function EVprime(w′::Float64, a=a, b=b, z=z, yn=yn, k=k, v_int=v_int, ρ=ρ)

    function EVp(y::Float64, w′=w′, v_int=v_int, yn=yn, a=a, b=b, z=z, ρ=ρ)
      dy = y - mean(yn)
      v_int[w′+exp(y), a+k[1]*dy, b+k[2]*dy, ρ*z+k[3]*dy]
    end

    (n, wgt) = gausshermite(10)
    π^(-0.5)*sum( [wgt[i]*EVp(sqrt(2)*std(yn)*n[i] + mean(yn))
                                                          for i = 1:length(n)] )
  end

  Blmn(w′::Float64, x=x, r=r, δ=δ) = -( u(x-w′, σ) + δ*EVprime(r*w′) )

  optimum = Optim.optimize(Blmn, wmin/r, x + abs(wmin/r) + 1.)
  w′ = optimum.minimum
  vopt = -(optimum.f_minimum)

  return w′, vopt
end

################################################################################

# Transition period, no habits
function bellOpt_TRANS{T<:AbstractFloat}(x::T, pension::T, wmin::T,
                  v_int::Interpolations.GriddedInterpolation, r::T, δ::T, σ::T)

  Blmn(w′) = -( u(x-w′, σ) + δ*v_int[r*w′, pension] )

  optimum = Optim.optimize(Blmn, wmin/r, x)
  w′ = optimum.minimum
  vopt = -(optimum.f_minimum)

  return w′, vopt
end
