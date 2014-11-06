###################################################################################
#########################    RETIREMENT PROBLEM      ##############################
###################################################################################
@printf "4. Solving the retirement problem\n"
tic()

V_R = zeros(wgridpoints_R, ygridpoints_R, TR)
wp_R = zeros(wgridpoints_R, ygridpoints_R, TR)

# Value of last period of retirement
for w = 1:wgridpoints_R
  for y = 1:ygridpoints_R
    wt = wgrid_R[w, TR]
    yt = ygrid_R[y]
    xt = wt + yt
    V_R[w, y, TR] = u(xt)
  end
end

for t = (TR-1):-1:1
  # INTERPOLATION
  V_R_interpol = interpolateV(V_R, wgrid_R, ygrid_R, t+1)

  # MAXIMIZATION
  wmin = wgrid_R[1, t+1]
  for w = 1:wgridpoints_R
    for y = 1:ygridpoints_R
      wt = wgrid_R[w, t]
      yt = ygrid_R[y]

      (wpopt, Vopt) = BellOpt_R(wt, yt, wmin, V_R_interpol, u, R, dbeta, t)

      V_R[w, y, t] = Vopt
      wp_R[w, y, t] = wpopt
    end
  end
end
@printf "\tSolving retirement problem took %.1f seconds.\n" toq()
