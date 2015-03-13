################################################################################
############################## BELLMAN RECURSION ###############################
################################################################################

function solveWorkingLife(v::Array{Float64, 5}, wp::Array{Float64, 5},
                          c_over_x::Array{Float64, 5},
                          wgrid::Array{Float64, 2}, agrid::Vector{Float64},
                          bgrid::Vector{Float64}, zgrid::Vector{Float64},
                          stdy::Vector{Float64}, k::Array, r::Float64,
                          δ::Float64)

  @printf "6. Recursively solve for optimal decision rules\n"
  @printf "\tSolving on %d points..." length(v[:, :, :, :, end])

  for t = (size(wgrid,2)-1):-1:1
    mod(t,10) != 0 || @printf "%d/%d, " t size(wgrid,2)

    # INTERPOLATION
    v_interpol =
      interpolatev_A(v[:, :, :, :, t+1], wgrid[:, t+1], agrid, bgrid, zgrid)

    # MAXIMIZATION
    wmin = wgrid[1, t+1]
    @inbounds for w = 1:size(wgrid,1)
      for a = 1:size(agrid,1)
        for b = 1:size(bgrid,1)
          for z = 1:size(zgrid,1)
            at = agrid[a]
            bt = bgrid[b]
            zt = zgrid[z]

            wt = wgrid[w, t]
            yt = exp(at + bt*t + zt)
            yln = LogNormal(at + bt*(t+1) + ρ*zt, stdy[t])

            (wt + yt - wmin/r > 0.01) || error("Cash on hand is $wt + $yt, BC is $wmin/r")

            (wpopt, vopt) =
              bellOpt(wt, yt, at, bt, zt, wmin, v_interpol, yln, k[:, t], ρ, r, δ)

            wpopt < wt + yt || @printf "\tw'>xt, w=%d,a=%d,b=%d,z=%d,t=%d\n" w a b z t

            c = wgrid[w, t] + yt - wpopt
            c_over_x[w, a, b, z, t] = c/(wgrid[w, t] + yt - wmin/r)

            wp[w, a, b, z, t] = wpopt
            v[w, a, b, z, t] = vopt
          end
        end
      end
    end
  end

  return v, wp, c_over_x
end

################################################################################

function solveWorkingLife(v::Array{Float64, 5}, wp::Array{Float64, 5},
                          c_over_x::Array{Float64, 5},
                          wgrid::Array{Float64, 2}, agrid::Vector{Float64},
                          bgrid::Vector{Float64}, zgrid::Vector{Float64},
                          stdy::Vector{Float64}, r::Float64,
                          δ::Float64)

  @printf "6. Recursively solve for optimal decision rules\n"
  @printf "\tSolving on %d points..." length(v[:, :, :, :, end])

  for t = (size(wgrid,2)-1):-1:1
    mod(t,10) != 0 || @printf "%d/%d, " t size(wgrid,2)

    # INTERPOLATION
    v_interpol =
      interpolatev_A(v[:, :, :, :, t+1], wgrid[:, t+1], agrid, bgrid, zgrid)

    # MAXIMIZATION
    wmin = wgrid[1, t+1]
    @inbounds for w = 1:size(wgrid,1)
      for a = 1:size(agrid,1)
        for b = 1:size(bgrid,1)
          for z = 1:size(zgrid,1)
            at = agrid[a]
            bt = bgrid[b]
            zt = zgrid[z]
            zgrid
            wt = wgrid[w, t]
            yt = exp(at + bt*t + zt)
            yln = LogNormal(at + bt*(t+1) + ρ*zt, stdy[t])

            (wt + yt - wmin/r > 0.01) || error("Cash on hand is $wt + $yt, BC is $wmin/r")

            (wpopt, vopt) =
              bellOpt(wt, yt, at, bt, zt, wmin, v_interpol, yln, r, δ)

            wpopt < wt + yt || @printf "\tw'>xt, w=%d,a=%d,b=%d,z=%d,t=%d\n" w a b z t

            c = wgrid[w, t] + yt - wpopt
            c_over_x[w, a, b, z, t] = c/(wgrid[w, t] + yt - wmin/r)

            wp[w, a, b, z, t] = wpopt
            v[w, a, b, z, t] = vopt
          end
        end
      end
    end
  end

  return v, wp, c_over_x
end

################################################################################

function solveWorkingLife(v::Array{Float64, 5}, wp::Array{Float64, 5},
                          c_over_x::Array{Float64, 5},
                          wgrid::Array{Float64, 2}, agrid::Array{Float64, 2},
                          bgrid::Array{Float64, 2}, zgrid::Array{Float64, 2},
                          stdy::Vector{Float64}, r::Float64, δ::Float64)

  @printf "6. Recursively solve for optimal decision rules\n"
  @printf "\tSolving on %d points..." length(v[:, :, :, :, end])

  for t = (size(wgrid,2)-1):-1:1
    mod(t,10) != 0 || @printf "%d/%d, " t size(wgrid,2)

    # INTERPOLATION
    v_interpol =
      interpolatev_A(v[:, :, :, :, t+1], wgrid[:, t+1], agrid[:, t+1],
                     bgrid[:, t+1], zgrid[:, t+1])

    # MAXIMIZATION
    wmin = wgrid[1, t+1]
    @inbounds for w = 1:size(wgrid,1)
      for a = 1:size(agrid,1)
        for b = 1:size(bgrid,1)
          for z = 1:size(zgrid,1)
            at = agrid[a, t]
            bt = bgrid[b, t]
            zt = zgrid[z, t]
            wt = wgrid[w, t]
            yt = exp(at + bt*t + zt)
            yln = LogNormal(at + bt*(t+1) + ρ*zt, stdy[t])

            (wt + yt - wmin/r > 0.01) || error("Cash on hand is $wt + $yt, BC is $wmin/r")

            (wpopt, vopt) =
              bellOpt(wt, yt, at, bt, zt, wmin, v_interpol, yln, k, r, δ)

            wpopt < wt + yt || @printf "\tw'>xt, w=%d,a=%d,b=%d,z=%d,t=%d\n" w a b z t

            c = wgrid[w, t] + yt - wpopt
            c_over_x[w, a, b, z, t] = c/(wgrid[w, t] + yt - wmin/r)

            wp[w, a, b, z, t] = wpopt
            v[w, a, b, z, t] = vopt
          end
        end
      end
    end
  end

  return v, wp, c_over_x
end
