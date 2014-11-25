using DataFrames, Dierckx, Distributions, Grid,  Optim, PyPlot, QuantEcon

# 1. Draw Income Distribution
(Yit, alpha, beta, ymedian) = IncomeDistribution(agents, bs, μₐ, μᵦ, var_a, var_b,
                                                   var_eps, ρ, var_eta, BR, TW)

# 2. Construct individual specific belief histories
(S_f_i, P_f, stdy) = Learning(agents, bs, TW, alpha, beta, Yit, ρ, var_eta, var_eps, init_z)

# 3. Construct Grids
(wgrid, hgrid, agrid, bgrid, zgrid, wgrid_R, hgrid_R, ygrid_R) =
    Grids(S_f_i, stdy, agents, bs, power, wgridpoints, hgridpoints, agridpoints,
          bgridpoints, zgridpoints, wgridpoints_R, hgridpoints_R, ygridpoints_R,
          wmaxR, R, TW, TR)

# 4. Solve Retirement Problem
(V_R, wp_R) = SolveRetirement(wgrid_R, hgrid_R, ygrid_R, wgridpoints_R, hgridpoints_R,
                              ygridpoints_R, u, R, δ, λ, TR, interpolateV, BellOpt_R)

# 5. Solve Transition Problem
(V, wp) = SolveTransition(V_R, wgrid_R, hgrid_R, ygrid_R, wgrid, hgrid, agrid, bgrid,
                          ymedian, wgridpoints, hgridpoints, agridpoints,
                          bgridpoints, zgridpoints, interpolateV, BellOpt_TRANS, u, R, δ, λ, TW)

# 6. Solve Working Life Problem
(V, wp) = SolveWorkingLife(V, wp, wgrid, hgrid, agrid, bgrid, zgrid,
                          wgridpoints, hgridpoints, agridpoints, bgridpoints,
                          zgridpoints, stdy, BellOpt, interpolateV, u, R, δ, λ, TW)

# 7. Simulate Economy
(c_t, h_t, w_t, wp_t) = Simulate(wp, wgrid, hgrid, agrid, bgrid, zgrid,
                                 Yit, S_f_i, R, interpolateV, λ, agents, bs, TW)
