################################################################################
############################# GRID CONStRUCTION ################################
################################################################################

function grids(s_f_i::Array{Float64, 3}, stdy::Array, wpoints::Int64,
  hpoints::Int64, apoints::Int64, bpoints::Int64, zpoints::Int64,
  wpoints_R::Int64, hpoints_R::Int64, ypoints_R, wmaxR::Float64, power::Float64,
  r::Float64, tR::Int64, guv_dist::Bool)

  println("Construct Grids")
  tW = size(s_f_i,3)

  ybelief = Array(Float64, (size(s_f_i,2), tW))
  ymaxbelief = Array(Float64, (tW, 1)); yminbelief = similar(ymaxbelief)
  for t = 1:tW
    for i = 1:size(s_f_i,2)
      ybelief[i, t] = exp(g_t[t] + [1.0, t, 1]'*s_f_i[:, i, t])[1]
    end
    ymaxbelief[t] = maximum(ybelief[:, t])*exp(3*stdy[t])
    yminbelief[t] = minimum(ybelief[:, t])/exp(3*stdy[t])
  end

  # Wealth and Habit Grids
  # Maximum wealth is given by three times the highest earnings
  # Minimum wealth is given by some ad hoc constraint

  xmin = Array(Float64, tW);
  xmax = 5*ymaxbelief
  xmin[tW] = -2*yminbelief[tW]
  for t = (tW-1):-1:1
    xmin[t] = xmin[t+1]/r - yminbelief[t]
  end

  xgrid = Array(Float64, (xpoints, tW)); xgridexp = similar(xgrid)
  hgrid = Array(Float64, (hpoints, tW)); hgridexp = similar(hgrid)
  for t = tW:-1:1
    xdistexp = (xmax[t] - xmin[t])^(1/power)
    hdistexp = (0.5*xmax[t] - 0.1)^(1/power)
    xinc = xdistexp/(xpoints-1)
    hinc = hdistexp/(hpoints-1)
    for i = 1:xpoints
      xgridexp[i, t] = (i-1)*xinc
    end
    for i = 1:hpoints
      hgridexp[i, t] = (i-1)*hinc
    end
    xgrid[:, t] = xgridexp[:, t].^power + xmin[t]
    hgrid[:, t] = hgridexp[:, t].^power + 0.1
  end

  # Belief grids
  agrid = collect(linspace(minimum(s_f_i[1, :, 2:tW]), maximum(s_f_i[1, :, :]), apoints))
  bgrid = collect(linspace(minimum(s_f_i[2, :, :]), maximum(s_f_i[2, :, :]), bpoints))
  zgrid = collect(linspace(minimum(s_f_i[3, :, :]), maximum(s_f_i[3, :, :]), zpoints))

  # Retirement grids
  wminR = 1.05*xgrid[1,end]
  wmaxR = 0.75*xgrid[end,end]
  wminR_phase = collect(linspace(wminR, 0.0, tR))
  wgrid_R = Array(Float64, (wpoints_R,tR))
  wgridexp_R = Array(Float64, wpoints_R)
  hgridexp_R = Array(Float64, hpoints_R)
  hgrid_R = Array(Float64, (hpoints_R, tR))

  for t = 1:tR
    wdistexp = (wmaxR - wminR_phase[t])^(1/power)
    hdistexp = (0.5*wmaxR - 0.1)^(1/power)
    winc = wdistexp/(wpoints_R-1)
    hinc = hdistexp/(hpoints_R-1)
    for i = 1:wpoints_R
      wgridexp_R[i] = (i-1)*winc
    end
    for i = 1:hpoints_R
      hgridexp_R[i] = (i-1)*hinc
    end
    wgrid_R[:, t] = wgridexp_R.^power + wminR_phase[t]
    hgrid_R[:, t] = hgridexp_R.^power + 0.1
  end

  ygrid_R = collect(linspace(0.8*minimum(pension), maximum(pension), ypoints_R))

  println("\tC-i-h grid: $(round([xgrid[1,1] xgrid[end,1]],2)) in period 1, "*
    "$(round([xgrid[1,end] xgrid[end,end]],2)) in period 40")
  println("\tBelief grids: α $(round([agrid[1] agrid[end]],2)), β "*
    "$(round([bgrid[1] bgrid[end]],2)), z $(round([zgrid[1] zgrid[end]],2))")
  println("\tRetirement grids: w_R $(round([wgrid_R[1,1] wgrid_R[end,1]],2)), "*
    "y_R $(round([ygrid_R[1,1] ygrid_R[end,1]],2))")

  return xgrid, hgrid, agrid, bgrid, zgrid, wgrid_R, hgrid_R, ygrid_R
end
