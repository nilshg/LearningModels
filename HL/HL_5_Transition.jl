###################################################################################
#########################    TRANSITION PROBLEM      ##############################
###################################################################################

function SolveTransition(V_R::Array, wgrid_R::Array, hgrid_R::Array, ygrid_R::Array,
                         wgrid::Array, hgrid::Array, agrid::Array, bgrid::Array,
                         ymedian::Float64,
                         wgridpoints::Int64, hgridpoints::Int64, agridpoints::Int64,
                         bgridpoints::Int64, zgridpoints::Int64, interpolateV::Function,
                         BellOpt_TRANS::Function, u::Function, R::Float64, δ::Float64,
                         λ::Float64, TW::Int64)

  @printf "5. Solving the problem for the last period of work\n"
  tic()
  wp = zeros(wgridpoints, hgridpoints, agridpoints, bgridpoints, zgridpoints, TW)
  V = zeros(wgridpoints, hgridpoints, agridpoints, bgridpoints, zgridpoints, TW)

  # INTERPOLATION
  valueRETIRE = interpolateV(V_R, wgrid_R, hgrid_R, ygrid_R, 1)

  # MAXIMIZATION
  wmin = wgrid_R[1, 1]

  for a = 1:agridpoints
    for b = 1:bgridpoints
      for z = 1:zgridpoints
        at = agrid[a]
        bt = bgrid[a]
        zt = 0

        yt = exp(at + bt*TW + zt)

        if yt < 0.3*ymedian
          yfixed = 0.9*yt
        elseif yt <= 2.0*ymedian
          yfixed = 0.27*ymedian + 0.32(yt - 0.3*ymedian)
        elseif yt <= 4.1*ymedian
          yfixed = 0.81*ymedian + 0.15*(yt - 2*ymedian)
        elseif yt > 4.1*ymedian
          yfixed = 1.1*yt
        end
        yfixed = 0.715*yfixed

        for w = 1: wgridpoints
          for h = 1:hgridpoints
            wt = wgrid[w, TW]
            ht = hgrid[h, TW]

            if wt + yfixed < wmin
              @printf "\tWarning: Cash in hand too low at w=%d, a=%d, b=%d, z=%d\n" w a b z
            end

            (wpopt, Vopt) = BellOpt_TRANS(wt, ht, yt, yfixed, wmin,
                                          valueRETIRE, u, R, δ, λ, TW)

            V[w, h, a, b, z, TW] = Vopt
            wp[w, h, a, b, z, TW] = wpopt
          end
        end
      end
    end
  end
  @printf "\tTransition period problem solved in %.1f seconds.\n" toq()

  return V, wp
end
