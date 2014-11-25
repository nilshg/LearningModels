################################################################################
###############################    SIMULATION    ###############################
################################################################################

function Simulate(wp::Array, wgrid::Array, hgrid::Array, agrid::Array, bgrid::Array,
                  zgrid::Array, Yit::Array, S_f_i::Array, R::Float64, interpolateV::Function,
                  λ::Float64, agents::Int64, bs::Int64, TW::Int64 )

  @printf "7. Simulate Consumption and Wealth Distribution\n"

  w_0 = 0
  h_0 = hgrid[round(hgridpoints/2), 1] + w_0/10 + rand(agents*bs)

  c_t = Array(Float64, (agents*bs, T))
  w_t = Array(Float64, (agents*bs, T))
  wp_t = Array(Float64, (agents*bs, T))
  h_t = Array(Float64, (agents*bs, T))

  w_t[:, 1] = w_0
  h_t[:, 2] = h_0

  for t = 2:T
    # INTERPOLATION
    wp_int = interpolateV(wp, wgrid, hgrid, agrid, bgrid, zgrid, t)

    # Bond Choice
    @printf "\tChoosing optimal consumption and savings for working life period %d\n" t
    for i = 1:agents*bs
      wt = w_t[i, t]
      ht = h_t[i, t]
      yt = Yit[i, t]
      (at, bt, zt) = S_f_i[:, i, t]

      xt = R*wt + yt

      wp_t[i, t] = wp_int[xt, ht, at, bt, zt]
      c_t[i, t] = xt - wp_t[i, t]

      if t < T
        w_t[i, t+1] = wp_t[i, t]
        h_t[i, t+1] = (1-λ)*ht + λ*c_t[i, t]
      end
    end
  end

  return c_t, h_t, w_t, wp_t
end
