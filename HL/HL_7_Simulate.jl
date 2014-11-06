############################
######## SIMULATION ########
############################

@printf "7. Simulate Consumption and Wealth Distribution"

w_0 = 0
h_0 = hgrid[round(hgridpoints/2), 1] + w_0/10 + rand(agents*bs)

c_t = zeros(agents*bs, T)
w_t = zeros(agents*bs, T)
wp_t = zeros(agents*bs, T)
h_t = zeros(agents*bs, T)

wp_coeffs = zeros(length(nreg), T);

w_t[:, 1] = w_0
h_t[:, 1] = h_0

for t = 35:40
  @printf "     Period: %d out of %d\n" t T

  wp_int = interpolatewp(wp, wgrid, hgrid, agrid, bgrid, zgrid)

  # Bond Choice

  for i = 1:agents*bs
    wt = w_t[i, t]
    ht = h_t[i, t]
    yt = Yit[i, t]
    (at, bt, zt) = S_f_i[:, i, t]



end
