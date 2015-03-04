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

function solveTransition(v_R::Array{Float64, 3}, wgrid_R::Array, ygrid_R::Array,
                         wgrid::Array{Float64, 2}, agrid::Array, bgrid::Array,
                         zgrid::Array, yit::Array, r::Float64, δ::Float64)

  @printf "5. Solving the problem for the last period of work\n"
  tic()
  tW = size(wgrid,2)
  wp = Array(Float64,
             (size(wgrid,1), size(agrid,1), size(bgrid,1), size(zgrid,1), tW))
  v = similar(wp)
  c_over_x = similar(wp)

  # INTERPOLATION
  valueRETIRE = interpolatev_A(v_R[:, :, 1], wgrid_R[:, 1], ygrid_R)

  # MAXIMIZATION
  wmin = wgrid_R[1, 1]

  for a = 1:size(agrid,1)
    for b = 1:size(bgrid,1)
      for z = 1:size(zgrid,1)
        size(a,2)==1 ? at = agrid[a] : at = agrid[a,end]
        size(b,2)==1 ? bt = bgrid[b] : bt = bgrid[b,end]
        zt = 0.0

        yt = exp(at + bt*tW + zt)

        ybari = mean(yit, 2)[:]
        (k_0, k_1) = linreg(yit[:, 40], ybari)
        avgy = mean(yit)

        function get_pension(y::Float64, k_0::Float64, k_1::Float64, avgy::Float64)
            ytilde = (k_0 + k_1*y)/avgy
            rratio = 0.0

            if ytilde < 0.3
                rratio = 0.9*ytilde
            elseif ytilde <= 2.0
                rratio = 0.27 + 0.32*(ytilde - 0.3)
            elseif ytilde <= 4.1
                rratio = 0.814 + 0.15*(ytilde - 2.0)
            else
                rratio = 1.129
            end
            return rratio*avgy
        end
        pension = get_pension(yt, k_0, k_1, avgy)

        for w = 1:size(wgrid,1)
          wt = wgrid[w, tW]

          (wpopt, vopt) =
            bellOpt_TRANS(wt, yt, pension, wmin, valueRETIRE, r, δ)

          wpopt < wt + yt || @printf "NC @ w=%d,a=%d,b=%d,z=%d" w a b z

          v[w, a, b, z, tW] = vopt
          wp[w, a, b, z, tW] = wpopt

          c = wgrid[w, tW] + yt - wp[w, a, b, z, tW]

          c_over_x[w, a, b, z, tW] =
            c/(wgrid[w, tW] - wmin/r + yt)
        end
      end
    end
  end
  (dim1, dim2, dim3, dim4) = checkmonotonicity(v, tW)
  dim1 + dim2 + dim3 + dim4 == 0 || @printf "\tMonotonicity violated in transition"
  @printf "\tTransition period problem solved in %.1f seconds.\n" toq()

  return v, wp, c_over_x
end

###################################################################################

function solveTransition(wgrid::Array{Float64, 2}, agrid::Array, bgrid::Array,
                         zgrid::Array, r::Float64, δ::Float64)
  tw = size(wgrid,2)

  xpoints = 200
  xmax = wgrid[end, 40]
  xmin = wgrid[1, 40]
  power = 2.0

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

  xg = Array{Float64, 1}[]
  push!(xg, xgrid_irr)
  v_int = lininterp(v_T, xg)

  function bellOpt(w::Float64, y::Float64, v_int::lininterp,
                   yln::LogNormal, r::Float64, δ::Float64, wmin::Float64)

    x = w + y

    function EVprime(w′::Float64, yln=yln, v_int=v_int)

      function EVp(y::Array{Float64,1}, w′=w′, v_int=v_int, yln=yln)
        result = similar(y)
        for i = 1:size(y, 1)
          result[i, :] = getValue(v_int, [r*w′ + y[i]])*pdf(yln, y[i])
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
               (size(wgrid,1), size(agrid,1), size(bgrid,1), size(zgrid,1), 40))
  v_1 = similar(xp_1)
  c_over_x = similar(xp_1)

  for w = 1:size(wgrid , 1)
    for a = 1:length(agrid)
      for b = 1:length(bgrid)
        for z = 1:length(zgrid)

          size(a,2)==1 ? at = agrid[a] : at = agrid[a,end]
          size(b,2)==1 ? bt = bgrid[b] : bt = bgrid[b,end]
          size(z,2)==1 ? zt = zgrid[z] : zt = zgrid[z,end]

          y = exp(agrid[a] + bgrid[b]*40 + zgrid[z])
          yln = LogNormal(at + bt*40 + zt, 0.22)

          (xp_1[w, a, b, z, end], v_1[w, a, b, z, end]) =
            bellOpt(wgrid[w, end], y, v_int, yln, r, δ, xgrid_irr[1])

          c = wgrid[w, end] + y - xp_1[w, a, b, z, end]

          c_over_x[w, a, b, z, end] =
            c/(wgrid[w, end] - xgrid_irr[1]/r + y)

          c > 0 || @printf "NC @ x=%d,a=%d,b=%d,z=%d\n" w a b z
        end
      end
    end
  end

  return xp_1, v_1, c_over_x
end
