################################################################################
############# OPTIMIZATIONS.JL - solve value function optimization #############
################################################################################

using ApproXD, Distributions, Grid, Optim, QuantEcon

################################################################################

# Working life with habits
function bellOpt{T<:AbstractFloat}(x::T, h::T, a::T, b::T, z::T, wmin::T,
  v_int::CoordInterpGrid, yln::LogNormal, k::Array, ρ::T, r::T, δ::T, λ::T)

  function EVprime(w′::Float64, h=h, a=a, b=b, z=z, yln=yln, k=k, v_int=v_int,
                   λ=λ, ρ=ρ, x=x)

    h′ = (1-λ)*h + λ*(x - w′)

    function EVp(y::Array{Float64,1}, w′=w′, h′=h′, v_int = v_int, yln = yln,
                 k=k, a=a, b=b, z=z, ρ=ρ)
      ey = meanlogx(yln)
      result = similar(y)
      @inbounds for i = 1:size(y, 1)
         result[i,:] =
          v_int[w′ + y[i], h′,
                a + k[1]*(y[i]- ey),
                b + k[2]*(y[i]- ey),
                ρ*z + k[3]*(y[i]- ey)]*pdf(yln, y[i])
      end
      result
    end

    y_l = meanlogx(yln) - 3*stdlogx(yln)
    y_h = meanlogx(yln) + 3*stdlogx(yln)

    quadrect(EVp, 9, exp(y_l), exp(y_h))
  end

  Blmn(w′::Float64, x=x, r=r, δ=δ) = -( u_h(x - w′, h) + δ*EVprime(r*w′) )

  optimum = optimize(Blmn, wmin/r, x + abs(wmin/r) + 0.001)
  w′ = optimum.minimum
  vopt = -(optimum.f_minimum)

  return w′, vopt
end

################################################################################

# Transition period with habits
function bellOpt_TRANS{T<:AbstractFloat}(x::T, h::T, pension::T, wmin::T,
  v_int::Lininterp, r::T, δ::T, λ::T)

  Blmn(w′::Float64) = -(u_h(x - w′, h) +
                      δ*getValue(v_int, [r*w′, (1-λ)*h + λ*(x-w′), pension])[1])

  optimum = optimize(Blmn, wmin/r, x)
  w′ = optimum.minimum
  vopt = -(optimum.f_minimum)

  return w′, vopt
end

################################################################################

# Retirement period with habits
function bellOpt_R{T<:AbstractFloat}(w::T, h::T, y::T, wmin::T,
  v_int::Lininterp, r::T, δ::T, λ::T)

  x = w + y

  Blmn(w′) = -(u_h(x - w′, h)
             + δ*getValue(v_int, [r*w′, (1-λ)*h + λ*(x-w′), y])[1])

  optimum = optimize(Blmn, wmin/r, x)
  w′ = optimum.minimum
  vopt = -(optimum.f_minimum)

  return w′, vopt
end

################################################################################

# Working life, no habits
function bellOpt{T<:AbstractFloat}(x::T, a::T, b::T, z::T, wmin::T,
   v_int::Lininterp, yln::LogNormal, k::Array, ρ::T, r::T, δ::T)

  function EVprime(w′::Float64, a=a, b=b, z=z, yln=yln, k=k, v_int=v_int, ρ=ρ)

    function EVp(y::Array{Float64,1}, w′=w′, v_int=v_int, yln=yln,
                 a=a, b=b, z=z, ρ=ρ)
      ey = meanlogx(yln)
      result = similar(y)
      for i = 1:size(y, 1)
        @inbounds result[i, :] =
          (getValue(v_int,[w′ + y[i],
            a + k[1]*(log(y[i]) - ey),
            b + k[2]*(log(y[i]) - ey),
            ρ*z + k[3]*(log(y[i]) - ey)])[1])*pdf(yln, y[i])
      end
      return result
    end

    y_l = meanlogx(yln) - 3*stdlogx(yln)
    y_h = meanlogx(yln) + 3*stdlogx(yln)

    quadrect(EVp, 9, exp(y_l), exp(y_h))
  end

  Blmn(w′::Float64, x=x, r=r, δ=δ) = -( u(x-w′) + δ*EVprime(r*w′) )

  optimum = optimize(Blmn, wmin/r, x + abs(wmin/r) + .001)
  w′ = optimum.minimum
  vopt = -(optimum.f_minimum)

  return w′, vopt
end

################################################################################

# Transition period, no habits
function bellOpt_TRANS{T<:AbstractFloat}(x::T, pension::T, wmin::T,
                                         v_int::Lininterp, r::T, δ::T)

  Blmn(w′) = -( u(x-w′) + δ*(getValue(v_int, [r*w′ + pension, pension])[1]) )

  optimum = optimize(Blmn, wmin/r, x)
  w′ = optimum.minimum
  vopt = -(optimum.f_minimum)

  return w′, vopt
end

################################################################################
