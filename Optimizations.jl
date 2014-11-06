################################################################################
############# OPTIMIZATIONS.JL - solve value function optimization #############
################################################################################

# Contains:
# BellOpt(w, h, a, b, z, wmin, V_int, Yln, u, R, lambda, dbeta, t):
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
                             R::T, lamdba::T, dbeta::T, t::Real)

   x = w + y

  function EVprime_leg(wp::Real, x = x, h = h, a=a, b=b, z=z, t=t,
                       Yln = Yln, V_int = V_int)
    y_l = y - 3*Yln.sdlog
    y_h = y + 3*Yln.sdlog

    hp = (1-lambda)*h + lambda*(x - wp/R)

    function EVp(y::Array{Float64,1}, wp = wp, hp = hp,
                 V_int = V_int, Yln = Yln, a=a, b=b, z=z)
      result = zeros(size(y, 1), 1)
      for i = 1:size(y, 1)
        result[i, :] = V_int[wp + y[i], hp, a, b, z]*pdf(Yln, y[i])
      end
      return result
    end

    function EVp(y::Float64, wp = wp, hp=hp,
                 V_int = V_int, Yln =Yln, a = a, b = b, z = z)
      return V_int[wp + y, hp, a, b, z]*pdf(Yln, y)
    end

    quadrect(EVp, 9, y_l, y_h)
  end

  Blmn(wp::Float64) = -( u(x - wp/R, (1 - lambda)*h + lambda*(x - wp/R))
                         + dbeta * EVprime_leg(wp) )

  if x - 0.01 <= wmin
    wprime = wmin
    Vopt = Blmn(wmin)
  else
    Optimum = optimize(Blmn, wmin, x - 0.01)
    wprime = Optimum.minimum
    Vopt = -(Optimum.f_minimum)
  end
  return wprime, Vopt
end

################################################################################

function BellOpt_TRANS{T<:Float64}(w::T, h::T, y::T, yfixed::T, wmin::T,
                                   V_int::CoordInterpGrid, u::Function,
                                   R::T, dbeta::T, lambda::T, t::Int64)

  x = w + y

  Blmn(wp) = -( u(x - wp/R, h)
                + dbeta * V_int[wp, (1-lambda)*h + lambda*(x-wp/R), yfixed])

  if x - 0.01 <= wmin
    wprime = wmin
    Vopt = -(Blmn(wmin))
    if h == 1
      @printf "\t At [w,y,q] = [%d,%d,%d] consumer is in a corner\n" w y q
    end
  else
    Optimum = optimize(Blmn, wmin, x - 0.01)
    wprime = Optimum.minimum
    Vopt = -(Optimum.f_minimum)
  end
  return wprime, Vopt
end

################################################################################

function BellOpt_R{T<:Float64}(w::T, h::T, y::T, wmin::T,
                   V_int::CoordInterpGrid, u::Function,
                   R::T, dbeta::T, lambda::T, t::Int64)

  x = w + y

  Blmn(wp) = -( u(x - wp/R, h)
                + dbeta * V_int[wp, (1-lambda)*h + lambda*(x-wp/R), y])

  Optimum = optimize(Blmn, wmin, x - 0.01)

  wprime = Optimum.minimum
  Vopt = -(Optimum.f_minimum)

  return wprime, Vopt
end

#######################################################################################

function BellOpt{T<:Float64}(wt::T, yt::T, at::T, bt::T, zt::T, wmin::T,
                             V_int::CoordInterpGrid, Yln::LogNormal, u::Function,
                             R::T, dbeta::T, t::Int64)

   xt = wt + yt

  function EVprime_leg(wp::Float64, x = xt, a=at, b=bt, z=zt, t=t,
                       Yln = Yln, V_int = V_int)
    y_l = yt - 3*Yln.sdlog
    y_h = yt + 3*Yln.sdlog

    function EVp(y::Array{Float64,1}, wp = wp, V_int = V_int, Yln =Yln, a = a, b = b, z = z)
      result = zeros(size(y, 1), 1)
      for i = 1:size(y, 1)
        result[i, :] = V_int[wp + y[i], a, b, z]*pdf(Yln, y[i])
      end
      return result
    end

    function EVp(y::Float64, wp = wp, V_int = V_int, Yln =Yln, a = a, b = b, z = z)
      return V_int[wp + y, a, b, z]*pdf(Yln, y)
    end

    quadrect(EVp, 9, y_l, y_h)
  end

  Blmn(wp) = -( u(xt - wp/R) + dbeta * EVprime_leg(wp) )

  if xt - 0.01 <= wmin
    wprime = wmin
    Vopt = Blmn(wmin)
  else
    Optimum = optimize(Blmn, wmin, xt - 0.01)
    wprime = Optimum.minimum
    Vopt = -(Optimum.f_minimum)
  end
  return wprime, Vopt
end

################################################################################

function BellOpt_TRANS{T<:Float64}(w::T, y::T, yfixed::T, wmin::T,
                   V_interpol::CoordInterpGrid, u::Function,
                   R::T, dbeta::T, t::Int64)

  x = w + y

  Blmn(wp::Float64) = -( u(x - wp/R) + dbeta * V_interpol[wp, yfixed])

  if x - 0.01 <= wmin
    wprime = wmin
    Vopt = -(Blmn(wmin))
  else
    Optimum = optimize(Blmn, wmin, x - 0.01)
    wprime = Optimum.minimum
    Vopt = -(Optimum.f_minimum)
  end
  return wprime, Vopt
end

################################################################################

function BellOpt_R{T<:Float64}(w::T, y::T, wmin::T,
                   V_int::CoordInterpGrid, u::Function,
                   R::T, dbeta::T, t::Int64)

  x = w + y

  Blmn(wp) = -( u(x - wp/R) + dbeta * V_int[wp, y])

  Optimum = optimize(Blmn, wmin, x - 0.01)

  wprime = Optimum.minimum
  Vopt = -(Optimum.f_minimum)

  return wprime, Vopt
end
