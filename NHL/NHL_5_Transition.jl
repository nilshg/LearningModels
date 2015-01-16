###################################################################################
#########################    TRANSITION PROBLEM      ##############################
###################################################################################
#
# CONTAINS:
#
# solveTransition(v_R, wgrid_R, ygrid_R, wgrid, agrid, bgrid, ymedian, r, δ, tW)
# - solves the transition problem as in Guvenen (2007) with z=0 and expected 
#   retirement value as E[V']
#
# solveTransition(wgrid, agrid, bgrid, zgrid, r, δ, tW)
# - solves the transition problem with a simple V_{T+1}=0 approach, no retirement
#
###################################################################################

using Dierckx

function solveTransition(v_R::Array, wgrid_R::Array, ygrid_R::Array,
                         wgrid::Array, agrid::Array, bgrid::Array, ymedian::Float64,
                         r::Float64, δ::Float64, tW::Int64)

  @printf "5. Solving the problem for the last period of work\n"
  tic()
  wp = Array(Float64,
             (size(wgrid,1), size(agrid,1), size(bgrid,1), size(zgrid,1), tW))
  v = similar(wp)

  # INTERPOLATION
  valueRETIRE = interpolatev(v_R, wgrid_R, ygrid_R, 1)

  # MAXIMIZATION
  wmin = wgrid_R[1, 1]

  for a = 1:length(agrid)
    for b = 1:length(bgrid)
      for z = 1:length(zgrid)
        size(a,2)==1 ? at = agrid[a] : at = agrid[a,t]
        size(b,2)==1 ? bt = bgrid[b] : bt = bgrid[b,t]
        zt = 0.0

        yt = exp(at + bt*tW + zt)

        if (yt<0.3*ymedian)
          yfixed=0.9*yt
        elseif (yt<=2.0*ymedian)
          yfixed=0.27*ymedian+0.32*(yt-0.3*ymedian)
        elseif (yt<4.1*ymedian)
          yfixed = 0.81*ymedian+0.15*(yt-2.0*ymedian)
        else
          yfixed = 1.1*ymedian
        end
        pension = 0.715*yfixed

        for w = 1:wpoints
          wt = wgrid[w, tW]

          (wpopt, vopt) =
            bellOpt_TRANS(wt, yt, pension, wmin, valueRETIRE, r, δ)

          wpopt < wt + yt || @printf "NC @ w=%d,a=%d,b=%d,z=%d" w a b z

          abs(wpopt - wmin/r) > 0.0001 || (constrained[tW] += 1)

          v[w, a, b, z, tW] = vopt
          wp[w, a, b, z, tW] = wpopt
        end
      end
    end
  end
  (dim1, dim2, dim3, dim4) = checkmonotonicity(v, tW)
  dim1 + dim2 + dim3 + dim4 == 0 || @printf "\tMonot. violated, t=%d\n" t
  @printf "\tTransition period problem solved in %.1f seconds.\n" toq()

  return v, wp, constrained
end

###################################################################################

function solveTransition(wgrid::Array, agrid::Array, bgrid::Array, zgrid::Array,
                         r::Float64, δ::Float64, tW::Int64)

  xpoints = 200
  xmax = wgrid[end, 40]
  xmin = wgrid[1, 40]
  power = 3.0

  xgrid_irr = Array(Float64, xpoints)
  xgrid_exp = similar(xgrid_irr)
  xdistexp = (xmax - xmin)^(1/power)
  xinc = xdistexp/(xpoints-1)
  for i = 1:xpoints
    xgrid_exp[i] = (i-1)*xinc
  end

  xgrid_irr = xgrid_exp.^power + xmin

  v_T = Array(Float64, xpoints)

  for x = 1:xpoints
    v_T[x] = u(xgrid_irr[x] - xgrid_irr[1] + 0.01)
  end

  v_int = Spline1D(xgrid_irr, v_T)

  function bellOpt(w::Float64, y::Float64, a::Float64, b::Float64,
                   z::Float64, v_int::Spline1D, yln::LogNormal,
                   r::Float64, δ::Float64, wmin::Float64)

    x = w + y

    function EVprime(w′::Float64, a=a, b=b, z=z, yln=yln, v_int=v_int)

      function EVp(y::Array{Float64,1}, w′=w′, v_int=v_int, yln=yln, a=a, b=b, z=z)
        result = Array(Float64, size(y, 1))
        for i = 1:size(y, 1)
          result[i, :] = evaluate(v_int, r*w′ + y[i])*pdf(yln, y[i])
        end
        return result
      end

      y_l = meanlogx(yln) - 3*varlogx(yln)
      y_h = meanlogx(yln) + 3*varlogx(yln)

      quadrect(EVp, 9, exp(y_l), exp(y_h))
    end

    Blmn(w′) = -( u(x - w′) + δ * EVprime(r*w′) )

    optimum = optimize(Blmn, wmin/r, x)
    w′ = optimum.minimum
    vopt = -(optimum.f_minimum)

    return w′, vopt
  end

  xp_1 = Array(Float64, 
               (size(wgrid,1), size(agrid,1), size(bgrid,1), size(zgrid,1), tW))
  c_1 = similar(xp_1)
  v_1 = similar(xp_1)
  c_over_x = similar(xp_1)

  for w = 1:size(wgrid , 1)
      for a = 1:length(agrid)
          for b = 1:length(bgrid)
              for z = 1:length(zgrid)
          
                  size(a,2)==1 ? at = agrid[a] : at = agrid[a,t]
                  size(b,2)==1 ? bt = bgrid[b] : bt = bgrid[b,t]
                  size(z,2)==1 ? zt = zgrid[z] : zt = zgrid[z,t]

                  y = exp(agrid[a] + bgrid[b]*40 + zgrid[z])
                  yln = LogNormal(at + bt*40 + zt, 0.22)

                  (xp_1[w, a, b, z, tW], v_1[w, a, b, z, tW]) =
                      bellOpt(wgrid[w, 39], y, at, bt, zt, v_int, yln, r, δ, 
                              xgrid_irr[1])

                  c_1[w, a, b, z, tW] = wgrid[w, 39] + y - xp_1[w, a, b, z, tW]

                  c_over_x[w, a, b, z, tW] = 
                    c_1[w, a, b, z, tW]/(wgrid[w, 39] - xgrid_irr[1]/r + y)

                  c_1[w, a, b, z, tW] > 0 || @printf "NC @ x=%d,a=%d,b=%d,z=%d\n" w a b z

              end
          end
      end
  end

  return xp_1, c_1, v_1, c_over_x
end
