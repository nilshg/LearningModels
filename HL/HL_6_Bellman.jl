################################################################################
############################## BELLMAN RECURSION ###############################
################################################################################

function SolveWorkingLife(V::Array, wp::Array, wgrid::Array, hgrid::Array,
                          agrid::Array, bgrid::Array, zgrid::Array,
                          wgridpoints::Int64, hgridpoints::Int64,
                          agridpoints::Int64, bgridpoints::Int64,
                          zgridpoints::Int64, stdy::Array, BellOpt::Function,
                          interpolateV::Function, u::Function, R::Float64,
                          δ::Float64, λ::Float64, TW::Int64)

  @printf "6. Recursively solve for optimal decision rules\n"
  @printf "\ta = [%.2f, %.2f], b = [%.2f, %.2f], z = [%.2f, %.2f]\n" agrid[1] agrid[end] bgrid[1] bgrid[end] zgrid[1] zgrid[end]
  @printf "\tSolving the problem on %d gridpoints\n" wgridpoints*hgridpoints*agridpoints*bgridpoints*zgridpoints

  for t = (TW-1):-1:1
    @printf "     Period: %d out of %d\n" t TW
    @printf "\tw = [%.2f, %.2f], h = [%.2f, %.2f]\n" wgrid[1, t] wgrid[end, t] hgrid[1, t] hgrid[end, t]


    # INTERPOLATION
    V_interpol = interpolateV(V, wgrid, hgrid, agrid, bgrid, zgrid, t+1)

    # MAXIMIZATION
    wmin = wgrid[1, t+1]
    tic()
    for a = 1:agridpoints
      for b = 1:bgridpoints
        for z = 1:zgridpoints
          for h = 1:hgridpoints
            for w = 1:wgridpoints
              ht = hgrid[h, t]
              wt = wgrid[w, t]
              at = agrid[a]
              bt = bgrid[b]
              zt = zgrid[z]

              yt = exp(at + bt*t + zt)

              Yln = LogNormal(at + b*(t+1) + zt, stdy[t])

              if wt + yt - wmin/R < 0.01
                @printf "Error: Cash on hand is too low, at w=%d, a=%d, b=%d, z=%d" w a b z
              end

              (wpopt, Vopt) = BellOpt(wt, ht, yt, at, bt, zt, wmin,
                                      V_interpol, Yln, u, R, δ, λ, t)

              wp[w, h, a, b, z, t] = wpopt
              V[w, h, a, b, z, t] = Vopt
            end
          end
        end
      end
    end
    @printf "\tMaximization took %d seconds\n" toq()
  end

  return V, wp
end
