################################################################################
############# OPTIMIZATIONS.JL - solve value function optimization #############
################################################################################

# Contains:
# bellOpt(w, h, a, b, z, wmin, v_int, yln, r, λ, δ, t):
#     Solve value function optimization for state (x,h,a,b,z) using
#     Gauss-Legendre Quadrature
#
# bellOpt_TRANS(w, h, y, pension, wmin, v_int, R, dbeta, λ, t):
#     Solve value function optimization for last period of working life
#     with retirement income given by pension
#
# bellOpt_R(w, h, y, wmin, v_int, r, δ, λ, t)
#     Solve retirement value function optimization with constant income
#
# bellOpt(w, a, b, z, wmin, v_int, yln, R, δ, t):
#     Solve optimization problem without habits
#
# bellOpt_TRANS(w, y, pension, wmin, v_int, R, δ, t)
#     Solve transition period optimizatoin problem without habits
#
# bellOpt_R(w, y, wmin, v_int, R, δ, t)
#     Solve retirement optimization problem without habits

################################################################################

using ApproXD, Distributions, Grid, Optim, QuantEcon

################################################################################

function bellOpt(w::Float64, h::Float64, y::Float64, a::Float64, b::Float64,
                 z::Float64, wmin::Float64, v_int::CoordInterpGrid,
                 yln::LogNormal, k::Array, r::Float64, λ::Float64, δ::Float64)

   x = w + y

  function EVprime(w′::Float64, h=h, a=a, b=b, z=z, yln=yln, k=k, v_int=v_int,
                   λ=λ)

    h′ = (1-λ)*h + λ*(x - w′)

    function EVp(y::Array{Float64,1}, w′=w′, h′=h′, v_int = v_int, yln = yln,
                 k=k, a=a, b=b, z=z)
      ey = meanlogx(yln)
      result = similar(y)
      @inbounds for i = 1:size(y, 1)
         result[i,:] =
          v_int[w′ + y[i], h′,
                a + k[1]*(y[i]- ey),
                b + k[2]*(y[i]- ey),
                z + k[3]*(y[i]- ey)]*pdf(yln, y[i])
      end
      result
    end

    y_l = meanlogx(yln) - 3*stdlogx(yln)
    y_h = meanlogx(yln) + 3*stdlogx(yln)

    quadrect(EVp, 9, exp(y_l), exp(y_h))
  end

  Blmn(w′::Float64, x=x, r=r, δ=δ) = -( u(x - w′, h) + δ*EVprime(r*w′) )

  optimum = optimize(Blmn, wmin/r, x)
  w′ = optimum.minimum
  vopt = -(optimum.f_minimum)

  return w′, vopt
end

################################################################################

function bellOpt_TRANS(w::Float64, h::Float64, y::Float64, pension::Float64,
                       wmin::Float64, v_int::CoordInterpGrid,
                       r::Float64, δ::Float64, λ::Float64)

  x = w + y

  Blmn(w′::Float64) = -(u(x - w′, h) +
                          δ*v_int[r*w′, (1-λ)*h + λ*(x-w′), pension])

  optimum = optimize(Blmn, wmin/r, x)
  w′ = optimum.minimum
  vopt = -(optimum.f_minimum)

  return w′, vopt
end

################################################################################

function bellOpt_R(w::Float64, h::Float64, y::Float64, wmin::Float64,
                   v_int::CoordInterpGrid, r::Float64, δ::Float64, λ::Float64)

  x = w + y

  Blmn(w′) = -(u(x - w′, h) + δ*v_int[r*w′, (1-λ)*h + λ*(x-w′), y])

  optimum = optimize(Blmn, wmin/r, x)
  w′ = optimum.minimum
  vopt = -(optimum.f_minimum)

  return w′, vopt
end

################################################################################

function bellOpt(w::Float64, y::Float64, a::Float64, b::Float64, z::Float64,
                 wmin::Float64, v_int::lininterp, yln::LogNormal, k::Array,
                 ρ::Float64, r::Float64, δ::Float64)

  x = w + y

  function EVprime(w′::Float64, a=a, b=b, z=z, yln=yln, k=k, v_int=v_int)

    function EVp(y::Array{Float64,1}, w′=w′, v_int=v_int, yln=yln, a=a, b=b, z=z)
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

  optimum = optimize(Blmn, wmin/r, x)
  w′ = optimum.minimum
  vopt = -(optimum.f_minimum)

  return w′, vopt
end

################################################################################

function bellOpt_TRANS(w::Float64, y::Float64, pension::Float64, wmin::Float64,
                       v_int::lininterp, r::Float64, δ::Float64)

  x = w + y

  Blmn(w′) = -( u(x-w′) + δ*(getValue(v_int, [r*w′ + pension, pension])[1]) )

  optimum = optimize(Blmn, wmin/r, x)
  w′ = optimum.minimum
  vopt = -(optimum.f_minimum)

  return w′, vopt
end

################################################################################


function bellOpt_R(w::Float64, y::Float64, wmin::Float64, v_int::CoordInterpGrid,
                   r::Float64, δ::Float64)

  x = w + y

  Blmn(w′) = -(u(x - w′) + δ*v_int[r*w′, y])

  optimum = optimize(Blmn, wmin/r, x)
  w′ = optimum.minimum
  vopt = -(optimum.f_minimum)

  return w′, vopt
end

################################################################################
