################################################################################
############################# GRID CONStRUCTION ################################
################################################################################

function grids{T<:Int}(s_f_i::Array{Float64,3}, stdy::Array{Float64,1},
  xpoints::T, apoints::T, bpoints::T, zpoints::T, wpoints_R::T, ypoints_R::T,
  power::Float64, r::Float64, tR::T, g_t::Array{Float64,1},
  pension::Array{Float64,1}, const_bel::Bool)

  println("Construct Grids")
  tW = size(s_f_i,3)

  ybelief = Array{Float64}(size(s_f_i,2), tW)
  ymaxbelief = Array{Float64}(tW, 1); yminbelief = similar(ymaxbelief)
  for t = 1:tW
    for i = 1:size(s_f_i,2)
      ybelief[i, t] = exp(g_t[t] + [1.0, t, 1]'*s_f_i[:, i, t])[1]
    end
    ymaxbelief[t] = maximum(ybelief[:, t])*exp(3*stdy[t])
    yminbelief[t] = minimum(ybelief[:, t])/exp(3*stdy[t])
  end

  # WEALTH GRID #
  # Maximum wealth is given by three times the highest earnings
  # Minimum wealth is given by some ad hoc constraint

  xmin = Array{Float64}(tW);
  xmax = 5*ymaxbelief
  xmin[tW] = -2*yminbelief[tW]
  for t = (tW-1):-1:1
    xmin[t] = xmin[t+1]/r - yminbelief[t]
  end

  xgrid = Array{Float64}(xpoints, tW); xgridexp = similar(xgrid)
  for t = tW:-1:1
    xdistexp = (xmax[t] - xmin[t])^(1/power)
    xinc = xdistexp/(xpoints-1)
    for i = 1:xpoints
      xgridexp[i, t] = (i-1)*xinc
    end
    xgrid[:, t] = xgridexp[:, t].^power + xmin[t]
  end

  # BELIEF GRIDS #
  (std(s_f_i[1,:,:]) > 0.01) || (apoints = 3)
  (std(s_f_i[2,:,:]) > 0.01) || (bpoints = 3)

  if const_bel
    agrid = collect(linspace(minimum(s_f_i[1, :, 2:tW]), maximum(s_f_i[1, :, :]), apoints))
    bgrid = collect(linspace(minimum(s_f_i[2, :, :]), maximum(s_f_i[2, :, :]), bpoints))
    zgrid = collect(linspace(minimum(s_f_i[3, :, :]), maximum(s_f_i[3, :, :]), zpoints))
  else
    agrid = Array{Float64}(apoints,tW)
    bgrid = Array{Float64}(bpoints,tW)
    zgrid = Array{Float64}(zpoints,tW)
    for t = 1:tW
      agrid[:, t] = collect(linspace(minimum(s_f_i[1, :, t]), maximum(s_f_i[1, :, t]), apoints))
      bgrid[:, t] = collect(linspace(minimum(s_f_i[2, :, t]), maximum(s_f_i[2, :, t]), bpoints))
      zgrid[:, t] = collect(linspace(minimum(s_f_i[3, :, t]), maximum(s_f_i[3, :, t]), zpoints))
    end
  end

  if std(s_f_i[1,:,:]) < 0.01
    agrid = collect(linspace(minimum(s_f_i[1, :, 2:tW])-0.001, maximum(s_f_i[1, :, :])+0.001, apoints))
    bgrid = collect(linspace(minimum(s_f_i[2, :, :])-0.001, maximum(s_f_i[2, :, :])+0.001, bpoints))
  end

  # RETIREMENT GRIDS #
  wminR = 1.05*xgrid[1,end]
  wmaxR = 0.75*xgrid[end,end]
  wgrid_R = Array{Float64}(wpoints_R)
  wgridexp_R = Array{Float64}(wpoints_R)
  wdistexp = (wmaxR - wminR)^(1/power)
  winc = wdistexp/(wpoints_R-1)
  for i = 1:wpoints_R
    wgridexp_R[i] = (i-1)*winc
  end
  wgrid_R = wgridexp_R.^power + wminR

  ygrid_R = collect(linspace(0.8*minimum(pension), maximum(pension), ypoints_R))

  println("\tC-i-h grid: $(round.([xgrid[1,1] xgrid[end,1]],2)) in period 1, "*
    "$(round.([xgrid[1,end] xgrid[end,end]],2)) in period 40")
  println("\tBelief grids: α $(round.([agrid[1] agrid[end]],2)), β "*
    "$(round.([bgrid[1] bgrid[end]],2)), z $(round.([zgrid[1] zgrid[end]],2))")
  println("\tRetirement grids: w_R $(round.([wgrid_R[1,1] wgrid_R[end,1]],2)), "*
    "y_R $(round.([ygrid_R[1,1] ygrid_R[end,1]],2))")

  return xgrid, agrid, bgrid, zgrid, wgrid_R, ygrid_R
end
