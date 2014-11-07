################################################################################
################################### SIMULATION #################################
################################################################################

@printf "7. Simulate Consumption and Wealth Distribution\n"

w_0 = 0

c_t = zeros(agents*bs, T)
w_t = zeros(agents*bs, T)
wp_t = zeros(agents*bs, T)

w_t[:, 1] = w_0

for t = 1:T
  @printf "     Period: %d out of %d\n" t T

  # INTERPOLATION
  tic()
  wp_int = interpolateV(wp, wgrid, agrid, bgrid, zgrid, t)
  @printf "\tInterpolation took %.2f seconds\n" toq()

  # Bond Choice
  @printf "\tChoosing optimal consumption and savings\n"
  tic()
  for i = 1:agents*bs
    wt = w_t[i, t]
    yt = Yit[i, t]
    (at, bt, zt) = S_f_i[:, i, t]

    xt = R*wt + yt

    wp_t[i, t] = wp_int[xt, at, bt, zt]
    c_t[i, t] = xt - wpopt
    if t < T
      w_t[i, t+1] = wpopt
    end
  end
  @printf "\tConsumption choice took %.1f seconds" toq()
end
