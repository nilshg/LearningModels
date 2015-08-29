################################################################################
############################# GRID CONStRUCTION ################################
################################################################################

function grids(s_f_i::Array{Float64, 3}, stdy::Array, wpoints::Int64,
  hpoints::Int64, apoints::Int64, bpoints::Int64, zpoints::Int64,
  wpoints_R::Int64, hpoints_R::Int64, ypoints_R, wmaxR::Float64, power::Float64,
  r::Float64, tR::Int64, guv_dist::Bool)

  @printf "3. Construct Grids\n"
  tW = size(s_f_i,3)
  hmat = [ones(1,tW); linspace(1,tW,tW)'; ones(1,tW)]
  ybelief = Array(Float64, (size(s_f_i,2), tW))
  ymaxbelief = Array(Float64, (tW, 1))
  yminbelief = similar(ymaxbelief)
  for t = 1:tW
    for i = 1:size(s_f_i,2)
      ybelief[i, t] = exp(hmat[:, t]'*s_f_i[:, i, t] + 3*stdy[t])[1]
    end
    ymaxbelief[t] = maximum(ybelief[:, t])*exp(3*stdy[t])
    yminbelief[t] = minimum(ybelief[:, t])/exp(3*stdy[t])
  end

  # Wealth grid
  # Maximum wealth is given by three times the highest earnings
  # Minimum wealth is given by some ad hoc constraint

  @printf "\t3.1 Wealth Grid\n"
  if guv_dist
    wgrid_org = readdlm("C:\\Users\\tew207\\Dropbox\\QMUL\\PhD\\Code\\Julia\\Guvenen\\wealth.dat")'
    wgrid = Array(Float64, (wpoints, tW)); wgridexp = similar(wgrid)

    for t = tW:-1:1
      wdistexp = (wgrid_org[end, t])^(1/power)# - wgrid_org[1, t])^(1/power)
      winc = wdistexp/(wpoints-1)
      for i = 1: wpoints
        wgridexp[i, t] = (i-1)*winc
      end
      wgrid[:, t] = wgridexp[:, t].^power# + wgrid_org[1, t]
    end
  else
    wmin = Array(Float64, tW)
    wmax = similar(wmin)

    wmin[tW] = 0.#-0.7*yminbelief[tW]
    wmax[tW] = 2*ymaxbelief[tW]
    for t = (tW-1):-1:1
      wmin[t] = wmin[t+1]/r # - 0.7*yminbelief[t]
      wmax[t] = 2*ymaxbelief[t]
    end

    wgrid = Array(Float64, (wpoints, tW))
    wgridexp = similar(wgrid)

    for t = tW:-1:1
      wdistexp = (wmax[t] - wmin[t])^(1/power)
      winc = wdistexp/(wpoints-1)
      for i = 1: wpoints
          wgridexp[i, t] = (i-1)*winc
      end
      wgrid[:, t] = wgridexp[:, t].^power + wmin[t]
    end
  end

  @printf "\t3.2 Habit Grid\n"
  hgrid = Array(Float64, (hpoints, tW))
  for t = 1:tW
      hgrid[:, t] = linspace(0.1, 1.5*ymaxbelief[t], hpoints)
  end

  @printf "\t3.3 Grids for α, β, z\n"
  agrid = convert(Array{Float64,1},
    linspace(minimum(s_f_i[1, :, 2:tW]), maximum(s_f_i[1, :, :]), apoints))
  bgrid = convert(Array{Float64,1},
    linspace(minimum(s_f_i[2, :, :]), maximum(s_f_i[2, :, :]), bpoints))
  zgrid = convert(Array{Float64,1},
    linspace(minimum(s_f_i[3, :, :]), maximum(s_f_i[3, :, :]), zpoints))

  @printf "\t3.4 Retirement Grids\n"
  if guv_dist # Use Guvenen's retirement grid
    guvgrid_R_org =
      readdlm("C:\\Users\\tew207\\Dropbox\\QMUL\\PhD\\Code\\Julia\\Guvenen\\wealthR.dat")'
    guvgrid_R = reshape(guvgrid_R_org, 1, 7200)
    guvgrid_R = unique(reshape(guvgrid_R_org, 12, 600), 2)
    wgrid_R = Array(Float64, (wpoints_R, tR))
    hgrid_R = Array(Float64, (hpoints_R, tR))
    for t = 1:30
      wgrid_R[:, t] = linspace(guvgrid_R[1, t], guvgrid_R[end, t], wpoints_R)
      hgrid_R[:, t] = linspace(0.1, 0.5*wgrid_R[end, t], hpoints_R)
    end
  else
    wminR = Array(Float64, (1, tR))
    wminR[tR] = 0.1059/(0.7)   # Directly from Guvenen's code

    for t = tR:-1:2
      wminR[t-1] = wminR[t]/r + wminR[tR] - 0.02
    end

    wminR = -0.7*wminR
    wgrid_R = Array(Float64, (wpoints_R, tR))
    wgridexp_R = similar(wgrid_R)
    hgrid_R = Array(Float64, (hpoints_R, tR))

    for t = 1:tR
      wdistexp = (wmaxR - wminR[t])^(1/power)
      winc = wdistexp/(wpoints_R-1)
      for i = 1:wpoints_R
        wgridexp_R[i, t] = (i-1)*winc
      end
      wgrid_R[:, t] = wgridexp_R[:, t].^power + wminR[t]
      hgrid_R[:, t] = linspace(0.1, 0.5*wgrid_R[end, t], hpoints_R)
    end
  end

  yminR = max(0.2*yminbelief[tW], 0.2)
  ymaxR = min(0.2*ymaxbelief[tW], 1000)
  ygrid_R = convert(Array{Float64, 1}, linspace(yminR, ymaxR, ypoints_R))

  @printf "Asset grid: [%.2f %.2f] in period 1, [%.2f %.2f] in period 40\n" wgrid[1,1] wgrid[end,1] wgrid[1,end] wgrid[end,end]
  @printf "Belief grids: α [%.2f %.2f], β [%.2f %.2f], z [%.2f %.2f]\n" agrid[1] agrid[end] bgrid[1] bgrid[end] zgrid[1] zgrid[end]
  @printf "Retirement grids: w_R [%.2f %.2f], y_R [%.2f %.2f]\n" wgrid_R[1, 5] wgrid_R[end, 5] ygrid_R[1] ygrid_R[end]

  return wgrid, hgrid, agrid, bgrid, zgrid, wgrid_R, hgrid_R, ygrid_R
end
