################################################################################
################################  DIAGNOSTICS  #################################
################################################################################

# Plot value and policy functions
function plot_value_policy(tw::Int64, tr::Int64)
  plotv(v, wgrid, agrid, bgrid, zgrid, "b", 3, 7, 1, tw,
      "Value Function, period "*string(tw))
  plotv(wp, wgrid, agrid, bgrid, zgrid, "b", 3, 5, 4, tw,
      "Policy Function, period "*string(tw))
  plotv(v_R, wgrid_R, ygrid_R, tr,
      "Retirement value function, period "*string(tr))
  plotv(wp_R, wgrid_R, ygrid_R, tr,
      "Retirement policy function, period "*string(tr))
end

## Simulation Results ##

function constrained_negative(w::Array{Float64,2}, wgrid::Array{Float64,2},
                              wgrid_R::Array{Float64,2})

  constrained = zeros(Float64, size(w,2))
  wmin = [wgrid[1, :] wgrid_R[1, :]]

  for t = 1:size(wgrid,2)
    for i = 1:size(wp,1)
      if abs(w[i, t] - wmin[t]) < 1e-3
        constrained[t] += 1/1000
      end
    end
  end

  neg_cons = similar(constrained)
  for t = 1:70
    for i = 1:size(w,1)
      if c_t[i, t] < 0
        neg_cons[t] += 1
      end
    end
  end
  neg_cons = neg_cons/1000

  return constrained, neg_cons
end

# Variance of consumption and asset series
function crosssec_stats(c::Array{Float64,2}, w::Array{Float64,2},
                            y::Array{Float64,2}, pension::Array)

  med_c = Array(Float64, size(c,2))
  med_w = similar(med_c)
  mean_w = similar(med_c)
  var_c = similar(med_c)
  var_w = similar(var_c)
  var_y = similar(var_c)

  for t = 1:size(c,2)
    med_c[t] = median(c[:, t])
    med_w[t] = median(w[:, t])
    mean_w[t] = mean(w[:, t])
    var_c[t] = var(c[:, t])
    var_w[t] = var(w[:, t])

    if t <= size(y,2)
      var_y[t] = var(y[:, t])
    else
      var_y[t] = var(pension)
    end
  end
  return med_c, med_w, mean_w, var_c, var_w, var_y
end

# Aggregate income, consumption, assets
function aggregates(c::Array{Float64,2}, w::Array{Float64,2},
                    y::Array{Float64,2}, pension::Array)
  agg_c = Array(Float64, size(c,2))
  agg_w = similar(agg_c)
  agg_y = similar(agg_c)

  for t = 1:size(c,2)
    agg_c[t] = sum(c[:, t])
    agg_w[t] = sum(w[:, t])
    if t <= tW
      agg_y[t] = sum(y[:, t])
    else
      agg_y[t] = sum(pension)
    end
  end
  return agg_c, agg_w, agg_y
end

# Beliefs vs grids vs realizations
function plot_beliefs_realizations()
  ybelief = Array(Float64, (agents*bs, tW))
  ymaxbelief = Array(Float64, (tW, 1))
  yminbelief = Array(Float64, (tW, 1))
  ymingridbelief = Array(Float64, tW)
  ymaxgridbelief = Array(Float64, tW)
  yminactual = Array(Float64, tW)
  ymaxactual = Array(Float64, tW)
  for t = 1:tW
    ymingridbelief[t] = exp(agrid[1] + t*bgrid[1] + zgrid[1])
    ymaxgridbelief[t] = exp(agrid[end] + t*bgrid[end] + zgrid[end])
    yminactual[t] = minimum(yit[:, t])
    ymaxactual[t] = maximum(yit[:, t])
    for i = 1:agents*bs
      ybelief[i, t] = exp([1, t, 1]'*s_f_i[:, i, t])[1]
    end
    ymaxbelief[t] = maximum(ybelief[:, t])
    yminbelief[t] = minimum(ybelief[:, t])
  end
  fig, ax = PyPlot.subplots(1, 2)
  ax[1,1][:plot](ymingridbelief, label = "Lowest belief on grid")
  ax[1,1][:plot](yminbelief, label = "Lowest actual belief")
  ax[1,1][:plot](yminactual, label = "Lowest actual income")
  ax[2,1][:plot](ymaxgridbelief, label = "Highest belief on grid")
  ax[2,1][:plot](ymaxbelief, label = "Highest actual belief")
  ax[2,1][:plot](ymaxactual, label = "Highest actual income")
  ax[1,1][:legend](loc="best")
  ax[2,1][:legend](loc="best")
  plt.show()
end

function plot2Dconfunc(c_x::Array{Float64, 5}, t::Int64,
                       figtext::String)
  α_mid = convert(Int64, round(size(c_x,2)/2))
  β_mid = convert(Int64, round(size(c_x,3)/2))
  z_mid = convert(Int64, round(size(c_x,4)/2))
  α_hi = size(c_x,2)
  β_hi = size(c_x,3)
  z_hi = size(c_x,4)

  fig, ax = PyPlot.subplots(2, 2)
  ax[1,1][:plot](c_x[:, 1, β_mid, z_mid, t], label = L"\alpha=1")
  ax[1,1][:plot](c_x[:, α_mid, β_mid, z_mid, t], label = L"\alpha=2")
  ax[1,1][:plot](c_x[:, end, β_mid, z_mid, t], label = L"\alpha="*string(α_hi))
  ax[1,1][:set_title](L"\beta="*string(β_mid)*", z="*string(z_mid))
  ax[1,1][:legend](loc = "best")
  ax[2,1][:plot](c_x[:, α_mid, 1, 4, t], label = L"\beta=1")
  ax[2,1][:plot](c_x[:, α_mid, 4, 4, t], label = L"\beta=4")
  ax[2,1][:plot](c_x[:, α_mid, 6, 4, t], label = L"\beta=6")
  ax[2,1][:plot](c_x[:, α_mid, end, 4, t], label = L"\beta="*string(β_hi))
  ax[2,1][:set_title](L"\alpha="*string(α_mid)*", z="*string(z_mid))
  ax[2,1][:legend](loc = "best")
  ax[1,2][:plot](c_x[:, α_mid, β_mid, 1, t], label = L"z=1")
  ax[1,2][:plot](c_x[:, α_mid, β_mid, 3, t], label = L"z=2")
  ax[1,2][:plot](c_x[:, α_mid, β_mid, 5, t], label = L"z=3")
  ax[1,2][:plot](c_x[:, α_mid, β_mid, end, t], label = L"z="*string(z_hi))
  ax[1,2][:set_title](L"\alpha="*string(α_mid)*", "*L"\beta="*string(β_mid))
  ax[1,2][:legend](loc = "best")
  ax[2,2][:plot](c_x[:, 1, 1, 1, t], label = L"\beta=1")
  ax[2,2][:plot](c_x[:, 1, 4, 1, t], label = L"\beta=4")
  ax[2,2][:plot](c_x[:, 1, 6, 1, t], label = L"\beta=6")
  ax[2,2][:plot](c_x[:, 1, 8, 1, t], label = L"\beta="*string(β_hi))
  ax[2,2][:set_title](L"\alpha=1, z=1")
  ax[2,2][:legend](loc = "best")
  fig[:suptitle]("Consumption over assets, period "*string(t))
  fig[:text](0.01, 0.01, figtext)
  plt.show()
end
