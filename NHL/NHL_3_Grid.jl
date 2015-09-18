################################################################################
############################# GRID CONStRUCTION ################################
################################################################################

function grids{T<:Int}(s_f_i::Array{Float64,3}, stdy::Array, xpoints::T,
  apoints::T, bpoints::T, zpoints::T, wpoints_R::T, ypoints_R::T,
  wmaxR::Float64, power::Float64, r::Float64, tR::T, const_bel::Bool)

  println("Construct Grids")
  tW = size(s_f_i,3)

  ybelief = Array(Float64, (size(s_f_i,2), tW))
  ymaxbelief = Array(Float64, (tW, 1)); yminbelief = similar(ymaxbelief)
  for t = 1:tW
    for i = 1:size(s_f_i,2)
      ybelief[i, t] = exp([1.0, t, 1]'*s_f_i[:, i, t] + 3*stdy[t])[1]
    end
    ymaxbelief[t] = maximum(ybelief[:, t])*exp(3*stdy[t])
    yminbelief[t] = minimum(ybelief[:, t])/exp(3*stdy[t])
  end

  # WEALTH GRID #
  # Maximum wealth is given by three times the highest earnings
  # Minimum wealth is given by some ad hoc constraint

  xmin = Array(Float64, tW); xmax = similar(xmin)
  xmin[tW] = -0.7*yminbelief[tW]
  for t = (tW-1):-1:1
    xmin[t] = xmin[t+1]/r - 0.5*yminbelief[t]
  end
  xmax = ymaxbelief

  xgrid = Array(Float64, (xpoints, tW)); xgridexp = similar(xgrid)
  for t = tW:-1:1
    xdistexp = (xmax[t] - xmin[t])^(1/power)
    xinc = xdistexp/(xpoints-1)
    for i = 1:xpoints
      xgridexp[i, t] = (i-1)*xinc
    end
    xgrid[:, t] = xgridexp[:, t].^power + xmin[t]
  end

  # BELIEF GRIDS #
  if const_bel
    agrid = convert(Array{Float64,1}, linspace(minimum(s_f_i[1, :, 2:tW]),
                             maximum(s_f_i[1, :, :]), apoints))
    bgrid = convert(Array{Float64,1}, linspace(minimum(s_f_i[2, :, :]),
                             maximum(s_f_i[2, :, :]), bpoints))
    zgrid = convert(Array{Float64,1}, linspace(minimum(s_f_i[3, :, :]),
                             maximum(s_f_i[3, :, :]), zpoints))
  else
    agrid = Array(Float64,(apoints,tW))
    bgrid = Array(Float64,(bpoints,tW))
    zgrid = Array(Float64,(zpoints,tW))
    for t = 1:tW
      agrid[:, t] = convert(Array{Float64,1}, linspace(minimum(s_f_i[1, :, t]),
                             maximum(s_f_i[1, :, t]), apoints))
      bgrid[:, t] = convert(Array{Float64,1}, linspace(minimum(s_f_i[2, :, t]),
                             maximum(s_f_i[2, :, t]), bpoints))
      zgrid[:, t] = convert(Array{Float64,1}, linspace(minimum(s_f_i[3, :, t]),
                             maximum(s_f_i[3, :, t]), zpoints))
    end
  end

  # RETIREMENT GRIDS #
  wminR = Array(Float64, (1, tR))
  wminR[tR] = 0.1059/(0.7)   # Directly from Guvenen's code
  for t = tR:-1:2
    wminR[t-1] = wminR[t]/r + wminR[tR] - 0.02
  end
  wminR = -0.7*wminR
  wgrid_R = Array(Float64, (wpoints_R, tR))
  wgridexp_R = Array(Float64, (wpoints_R, tR))
  for t = 1:tR
    wdistexp = (wmaxR - wminR[t])^(1/power)
    winc = wdistexp/(wpoints_R-1)
    for i = 1:wpoints_R
      wgridexp_R[i, t] = (i-1)*winc
    end
    wgrid_R[:, t] = wgridexp_R[:, t].^power + wminR[t]
  end

  ygrid_R = convert(Array, linspace(4., 15., ypoints_R))

  println("\tC-i-h grid: $(round([xgrid[1,1] xgrid[end,1]],2)) in period 1, "*
    "$(round([xgrid[1,end] xgrid[end,end]],2)) in period 40")
  println("\tBelief grids: α $(round([agrid[1] agrid[end]],2)), β "*
    "$(round([bgrid[1] bgrid[end]],2)), z $(round([zgrid[1] zgrid[end]],2))")
  println("\tRetirement grids: w_R $(round([wgrid_R[1,1] wgrid_R[end,1]],2)), "*
    "y_R $(round([ygrid_R[1,1] ygrid_R[end,1]],2))")

  return xgrid, agrid, bgrid, zgrid, wgrid_R, ygrid_R
end

################################################################################

function grids{T<:Int}(xpoints::T, apoints::T, bpoints::T, zpoints::T,
  wpoints_R::T, ypoints_R::T, power::Float64, user::AbstractString)

  println("Construct Grids (using Guvenen's data)")
  path="C:/Users/"*user*"/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/"
  # WEALTH GRID #
  wgrid_org = readdlm(path*"wealth.dat")'
  xgrid = Array(Float64, (xpoints, 40)); xgridexp = similar(xgrid)
  for t = 40:-1:1
    xdistexp = (wgrid_org[end, t] - wgrid_org[1, t])^(1./power)
    xinc = xdistexp/(xpoints-1)
    for i = 1: xpoints
      xgridexp[i, t] = (i-1)*xinc
    end
    xgrid[:, t] = xgridexp[:, t].^power + wgrid_org[1, t]
  end

  # BELIEF GRIDS #
  agrid = convert(Array{Float64,1}, linspace(1.92,2.09, apoints))
  bgrid = convert(Array{Float64,1}, linspace(-0.05,0.08, bpoints))
  zgrid = convert(Array{Float64,1}, linspace(-0.86,1.04, zpoints))

  # RETIREMENT GRIDS #
  guvgrid_R_org = readdlm(path*"wealthR.dat")'
  guvgrid_R = reshape(guvgrid_R_org, 1, 7200)
  guvgrid_R = unique(reshape(guvgrid_R_org, 12, 600), 2)
  wgrid_R = Array(Float64, (wpoints_R, 30)); wgridexp = similar(wgrid_R)

  for t = 1:30
    wdistexp = sqrt(guvgrid_R[end, t] - guvgrid_R[1, t])
    winc = wdistexp/(wpoints_R-1)
    for i = 1:wpoints_R
      wgridexp[i,t] = (i-1)*winc
    end
    wgrid_R[:, t] = wgridexp[:, t].^2 + guvgrid_R[1,t]
  end

  ygrid_R = convert(Array{Float64,1}, linspace(4.,15., ypoints_R))

  println("\tC-i-h grid: $(round([xgrid[1,1] xgrid[end,1]],2)) in period 1, "*
    "$(round([xgrid[1,end] xgrid[end,end]],2)) in period 40")
  println("\tBelief grids: α $(round([agrid[1] agrid[end]],2)), β "*
    "$(round([bgrid[1] bgrid[end]],2)), z $(round([zgrid[1] zgrid[end]],2))")
  println("\tRetirement grids: w_R $(round([wgrid_R[1,1] wgrid_R[end,1]],2)), "*
    "y_R $(round([ygrid_R[1,1] ygrid_R[end,1]],2))")

  return xgrid, agrid, bgrid, zgrid, wgrid_R, ygrid_R
end
