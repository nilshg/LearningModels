################################################################################
########################### INCOME DISTRIBUTION ################################
################################################################################

using Distributions

#################################################################################

function incomeDistribution(ypath::String, abpath::String)

  @printf "1. Import Guvenen's income distribution\n"
  yit = readdlm(ypath)
  alfabeta = readdlm(abpath)
  α = alfabeta[:, 1]
  α = reshape(repmat(α, 1, 100)', 100000, 1)
  β = alfabeta[:, 2]
  β = reshape(repmat(β, 1, 100)', 100000, 1)

  # Median income in last period for calculation of retirement benefits
  ymedian = median(yit[:, end])
  @printf "\tMedian income in period 40 is %.2f\n" ymedian

  ybari = mean(yit, 2)[:]
  (γ_0, γ_1) = linreg(yit[:, 40], ybari)

  function get_pension(y::Float64, k_0::Float64, k_1::Float64,
                       avgy::Float64)
      ytilde = (k_0 + k_1*y)/avgy
      rratio = 0.0

      if ytilde < 0.3
          rratio = 0.9*ytilde
      elseif ytilde <= 2.0
          rratio = 0.27 + 0.32*(ytilde - 0.3)
      elseif ytilde <= 4.1
          rratio = 0.814 + 0.15*(ytilde - 2.0)
      else
          rratio = 1.129
      end

      return rratio*avgy
  end

  pension = Array(Float64, size(yit,1))
  for i = 1:size(yit, 1)
    pension[i] = get_pension(yit[i, 40], γ_0, γ_1, mean(yit))
  end

  return yit, ymedian, pension, α, β
end

#################################################################################

function incomeDistribution(agents::Int64, bs::Int64, μₐ::Float64, μᵦ::Float64,
               var_α::Float64, var_β::Float64, cov_αβ::Float64, var_ɛ::Float64,
               var_η::Float64, ρ::Float64, y_adj::Float64, tW::Int64)

  min_β = max(-0.05, μᵦ-2.5*sqrt(var_β))

  @printf "1. Draw an income distribution\n"
  # Draw some alphas and betas
  ab = MvNormal([μₐ; μᵦ], [var_α cov_αβ; cov_αβ var_β])
  draw = rand(ab, bs)'
  for i = 1:bs
    while draw[i,2] < min_β
      draw[i,2] = min_β + abs(draw[i,2]-min_β)/50.0
    end
  end

  α = reshape(repmat(draw[:, 1],1,100)', agents*bs, 1)
  β = reshape(repmat(draw[:, 2],1,100)', agents*bs, 1)

  # Draw the income distribution:
  yit = zeros(bs*agents, tW)
  z = zeros(bs*agents, tW)
  z[:, 1] = sqrt(var_η/(1-ρ^2))*randn(agents*bs)

  for t = 1:tW, i = 1:bs*agents
    yit[i, t] = y_adj + exp(α[i] + β[i]*t + z[i, t] + sqrt(var_ɛ)*randn())
    if t < tW
      z[i, t+1] = ρ*z[i, t] + sqrt(var_η)*randn()
    end
  end

  # Median income in last period for calculation of retirement benefits
  ymedian = median(yit[:, end])
  @printf "\tMedian income in period 40 is %.2f\n" ymedian

  ybari = mean(yit, 2)[:]
  (γ_0, γ_1) = linreg(yit[:, 40], ybari)

  function get_pension(y::Float64, k_0::Float64, k_1::Float64,
                       avgy::Float64)
      ytilde = (k_0 + k_1*y)/avgy
      rratio = 0.0

      if ytilde < 0.3
          rratio = 0.9*ytilde
      elseif ytilde <= 2.0
          rratio = 0.27 + 0.32*(ytilde - 0.3)
      elseif ytilde <= 4.1
          rratio = 0.814 + 0.15*(ytilde - 2.0)
      else
          rratio = 1.129
      end

      return rratio*avgy
  end

  pension = Array(Float64, size(yit,1))
  for i = 1:size(yit, 1)
    pension[i] = get_pension(yit[i, 40], γ_0, γ_1, mean(yit))
  end

  return yit, ymedian, pension, α, β
end

################################################################################

function incomeDistribution(α::Float64, β::Float64, var_η_RIP::Float64,
                            var_ɛ_RIP::Float64, zgridpoints_RIP::Int64,
                            epsgridpoints::Int64, tW::Int64)

  @printf "1. Drawing an RIP-consistent income distribution\n"
  @printf "\tσ²(η) = %.3f,  σ²(ɛ) =  %.3f\n" var_η_RIP var_ɛ_RIP

  zdisc = [-(1/2)*(zgridpoints_RIP-1)*sqrt(var_η_RIP)+(i-1)*sqrt(var_η_RIP)
               for i = 1:zgridpoints_RIP]

  epsdisc = [-(1/2)*(epsgridpoints-1)*sqrt(var_ɛ_RIP)+(i-1)*sqrt(var_ɛ_RIP)
               for i = 1:epsgridpoints]

  yit = Array(Float64, (zgridpoints_RIP, epsgridpoints, tW))

  for t = 1:tW, z = 1:zgridpoints_RIP, ɛ = 1:epsgridpoints
    if t < 2
      yit[:, ɛ, t] = exp(α + β*t + 0.8*zdisc[z] + 0.6*epsdisc[ɛ])
    elseif t < 3
      yit[:, ɛ, t] = exp(α + β*t + zdisc[z] + epsdisc[ɛ])
    elseif t < 5
      yit[:, ɛ, t] = exp(α + β*t + zdisc[z] + epsdisc[ɛ])
    else
      yit[:, ɛ, t] = exp(α + β*t + zdisc[z] + epsdisc[ɛ])
    end
  end

  return yit
end
