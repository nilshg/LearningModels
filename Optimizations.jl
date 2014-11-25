################################################################################
############# OPTIMIZATIONS.JL - solve value function optimization #############
################################################################################

# Contains:
# BellOpt(w, h, a, b, z, wmin, V_int, Yln, u, R, lambda, \delta, t):
#     Solve value function optimization for state (x,h,a,b,z) using
#     Gauss-Legendre Quadrature
#
# BellOpt_TRANS(w, h, y, yfixed, wmin, V_int, u, R, dbeta, lambda, t):
#     Solve value function optimization for last period of working life
#     with retirement income given by yfixed
#
# BellOpt_R(w, h, y, wmin, V_int, u, R, dbeta, lambda, t)
#     Solve retirement value function optimization with constant income
#
# BellOpt(w, a, b, z, wmin, V_int, Yln, u, R, dbeta, t):
#     Solve optimization problem without habits
#
# BellOpt_TRANS(w, y, yfixed, wmin, V_interpol, u, R, dbeta, t)
#     Solve transition period optimizatoin problem without habits
#
# BellOpt_R(w, y, wmin, V_int, u, R, dbeta, t)
#     Solve retirement optimization problem without habits

################################################################################

function BellOpt{T<:Float64}(w::T, h::T, y::T, a::T, b::T, z::T, wmin::T,
                             V_int::CoordInterpGrid, Yln::LogNormal, u::Function,
                             R::T, λ::T, δ::T, t::Real)

   x = w + y

  function EVprime_leg(w′::Real, x=x, h=h, a=a, b=b, z=z, t=t, Yln=Yln, V_int=V_int)

    y_l = Yln.meanlog - 3*Yln.sdlog
    y_h = Yln.meanlog + 3*Yln.sdlog

    h′ = (1-λ)*h + λ*(x - w′/R)

    function EVp(y::Array{Float64,1}, w′=w′, h′=h′,
                 V_int = V_int, Yln = Yln, a=a, b=b, z=z)
      result = zeros(size(y, 1), 1)
      for i = 1:size(y, 1)
        result[i, :] = V_int[w′ + y[i], h′, a, b, z]*pdf(Yln, y[i])
      end
      return result
    end

    function EVp(y::Float64, w′=w′, h′=h′,
                 V_int = V_int, Yln =Yln, a = a, b = b, z = z)
      return V_int[w′ + y, h′, a, b, z]*pdf(Yln, y)
    end

    quadrect(EVp, 9, exp(y_l), exp(y_h))
  end

  Blmn(w′::Float64) = -(u(x-w′/R, (1-λ)*h + λ*(x-w′/R)) + δ*EVprime_leg(w′))

  Optimum = optimize(Blmn, wmin, x - 0.01)
  w′ = Optimum.minimum
  Vopt = -(Optimum.f_minimum)

  return w′  , Vopt
end

################################################################################

function BellOpt_TRANS{T<:Float64}(w::T, h::T, y::T, yfixed::T, wmin::T,
                                   V_int::CoordInterpGrid, u::Function,
                                   R::T, δ::T, λ::T, t::Int64)

  x = w + y

  Blmn(w′::Float64) = -(u(x - w′/R, h) + δ*V_int[w′, (1-λ)*h + λ*(x-w′/R), yfixed])

  Optimum = optimize(Blmn, wmin, x - 0.01)
  w′ = Optimum.minimum
  Vopt = -(Optimum.f_minimum)

  return w′, Vopt
end

################################################################################

function BellOpt_R{T<:Float64}(w::T, h::T, y::T, wmin::T,
                   V_int::CoordInterpGrid, u::Function,
                   R::T, δ::T, λ::T, t::Int64)

  x = w + y

  Blmn(w′) = -(u(x - w′/R, h) + δ*V_int[w′, (1-λ)*h + λ*(x-w′/R), y])

  Optimum = optimize(Blmn, wmin, x - 0.01)
  w′ = Optimum.minimum
  Vopt = -(Optimum.f_minimum)

  return w′, Vopt
end

#######################################################################################

function BellOpt{T<:Float64}(wt::T, yt::T, a::T, b::T, z::T, wmin::T,
                             V_int::CoordInterpGrid, Yln::LogNormal, u::Function,
                             R::T, δ::T, t::Int64)

  x = wt + yt

  function EVprime_leg(w′::Float64, a=a, b=b, z=z, t=t, Yln=Yln, V_int=V_int)
    y_l = Yln.meanlog - 3*Yln.sdlog
    y_h = Yln.meanlog + 3*Yln.sdlog

    function EVp(y::Array{Float64,1}, w′=w′, V_int=V_int, Yln=Yln, a=a, b=b, z=z)
      result = zeros(size(y, 1), 1)
      for i = 1:size(y, 1)
        result[i, :] = V_int[w′ + y[i], a, b, z]*pdf(Yln, y[i])
      end
      return result
    end

    function EVp(y::Float64, w′=w′, V_int=V_int, Yln=Yln, a=a, b=b, z=z)
      return V_int[w′ + y, a, b, z]*pdf(Yln, y)
    end

    quadrect(EVp, 9, exp(y_l), exp(y_h))
  end

  Blmn(w′) = -( u(x - w′/R) + δ * EVprime_leg(w′) )

  Optimum = optimize(Blmn, wmin, t - 0.01)
  w′ = Optimum.minimum
  Vopt = -(Optimum.f_minimum)

  return w′, Vopt
end

################################################################################

function BellOpt_TRANS{T<:Float64}(w::T, y::T, yfixed::T, wmin::T,
                   V_int::CoordInterpGrid, u::Function,
                   R::T, δ::T, t::Int64)

  x = w + y

  Blmn(w′) = -( u(x - w′/R) + δ*V_int[w′, yfixed])

  Optimum = optimize(Blmn, wmin, x - 0.01)
  w′ = Optimum.minimum
  Vopt = -(Optimum.f_minimum)

  return w′, Vopt
end

################################################################################

function BellOpt_R{T<:Float64}(w::T, y::T, wmin::T,
                   V_int::CoordInterpGrid, u::Function,
                   R::T, δ::T, t::Int64)

  x = w + y

  Blmn(w′) = -(u(x-w′/R) + δ*V_int[w′, y])

  Optimum = optimize(Blmn, wmin, x - 0.01)
  w′ = Optimum.minimum
  Vopt = -(Optimum.f_minimum)

  return w′, Vopt
end

################################################################################

function BellOpt_R{T<:Float64}(w::T, y::T, wmin::T,
                   V_int::Spline2D, u::Function,
                   R::T, δ::T, t::Int64)

  x = w + y

  Blmn(w′) = -(u(x-w′/R) + δ*evaluate(V_int, w′, y))

  Optimum = optimize(Blmn, wmin, x - 0.01)
  w′ = Optimum.minimum
  Vopt = -(Optimum.f_minimum)

  return w′, Vopt
end
