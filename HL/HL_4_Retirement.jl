###################################################################################
#########################    RETIREMENT PROBLEM      ##############################
###################################################################################
@printf "4. Solving the retirement problem\n"

V_R = zeros(wgridpoints, hgridpoints, ygridpoints_R, TR)
wp_R = zeros(wgridpoints, hgridpoints, ygridpoints_R, TR)

# Value of last period of retirement
for w = 1:wgridpoints
  for h = 1:hgridpoints
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

  V_R_interpol = interpolateV_R(V_R, wgrid_R, hgrid_R, ygrid_R, t+1)

  for w = 1:wgridpoints
    for h = 1:hgridpoints
      for y = 1:ygridpoints_R

        wmin = wgrid_R[1, t+1]

        wt = wgrid_R[w, t]
        ht = hgrid_R[h, t]
        yt = ygrid_R[y]

        (wpopt, Vopt) = BellOpt_R(wt, ht, yt, wmin, V_R_interpol, u, R, dbeta, lambda, t)

        V_R[w, h, y, t-1] = Vopt
        wp_R[w, h, y, t-1] = wpopt
      end
    end
  end
end
