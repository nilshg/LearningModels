###################################################################################
#########################    RETIREMENT PROBLEM      ##############################
###################################################################################

function SolveRetirement(wgrid_R::Array, hgrid_R::Array, ygrid_R::Array,
                         wgridpoints_R::Int64, hgridpoints_R::Int64, ygridpoints_R::Int64,
                         u::Function, R::Float64, δ::Float64, λ::Float64, TR::Int64,
                         interpolateV::Function, BellOpt_R::Function)

  @printf "4. Solving the retirement problem on %d gridpoints\n" wgridpoints_R*hgridpoints_R*ygridpoints_R
  tic()

  V_R = zeros(wgridpoints_R, hgridpoints_R, ygridpoints_R, TR)
  wp_R = zeros(wgridpoints_R, hgridpoints_R, ygridpoints_R, TR)

  # Value of last period of retirement
  for w = 1:wgridpoints_R
    for h = 1:hgridpoints_R
      for y = 1:ygridpoints_R
        wt = wgrid_R[w, TR]
        ht = hgrid_R[h, TR]
        yt = ygrid_R[y]
        xt = wt + yt
        V_R[w, h, y, TR] = u(xt, ht)
      end
    end
  end

  for t = (TR-1):-1:1
    # INTERPOLATION
    V_R_interpol = interpolateV(V_R, wgrid_R, hgrid_R, ygrid_R, t+1)

    # MAXIMIZATION
    wmin = wgrid_R[1, t+1]
    for w = 1:wgridpoints_R
      for h = 1:hgridpoints_R
        for y = 1:ygridpoints_R
          wt = wgrid_R[w, t]
          ht = hgrid_R[h, t]
          yt = ygrid_R[y]

          (wpopt, Vopt) = BellOpt_R(wt, ht, yt, wmin, V_R_interpol, u, R, δ, λ, t)

          V_R[w, h, y, t] = Vopt
          wp_R[w, h, y, t] = wpopt
        end
      end
    end
  end
  @printf "\tSolving retirement problem took %.1f seconds.\n" toq()

  return V_R, wp_R
end
