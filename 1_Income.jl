################################################################################
########################### INCOME DISTRIBUTION ################################
################################################################################
#
# Contains:
#
# (yit, alpha, beta, ymedian) = incomeDistribution(ypath, abpath):
#   - Get yit from Guvenen's data
#
# (yit, alpha, beta, ymedian) = incomeDistribution(agents, μₐ, μᵦ, var_a, var_b,
#                                                  var_eps, ρ, var_eta, br, tW)
#   - Draw a new income distribution yit with break parameters br
#
#################################################################################

using DataFrames

#################################################################################

function incomeDistribution(ypath::String, abpath::String)

  @printf "1. Draw an Income Distribution\n"
  @printf "\tWe are working with Guvenen's data\n"
  yit = readdlm(ypath)
  alfabeta = readdlm(abpath)
  α = alfabeta[:, 1]
  α = reshape(repmat(α, 1, agents)', 100000, 1)
  β = alfabeta[:, 2]
  β = reshape(repmat(β, 1, agents)', 100000, 1)

  # Calculate median income in last period for calculation of retirement benefits
  ymedian = median(yit[:, end])
  @printf "\tMedian income in period 40 is %.2f\n" ymedian

  pension = Array(Float64, agents*bs)
  for i = 1:100000 # Directly copied out of Guvenen's code
    ytemp = yit[i, end]
    if (ytemp<0.3*ymedian)
      yfixed=0.9*ytemp
    elseif (ytemp<=2.0*ymedian)
      yfixed=0.27*ymedian+0.32*(ytemp-0.3*ymedian)
    elseif (ytemp<4.1*ymedian)
      yfixed = 0.81*ymedian+0.15*(ytemp-2.0*ymedian)
    else
      yfixed = 1.1*ymedian
    end
    pension[i] = 0.715*yfixed
  end

  return yit, α, β, ymedian, pension
end

#################################################################################

function incomeDistribution(agents::Int64, bs::Int64, μₐ::Float64, μᵦ::Float64,
                            var_a::Float64, var_b::Float64, var_ɛ::Float64,
                            ρ::Float64, var_η::Float64, br::Int64, tW::Int64)

  # Draw some alphas and betas
  α = Array(Float64, bs)
  β_1 = Array(Float64, bs)
  for i = 1:bs
      α[i] = μₐ + sqrt(var_a)*randn()
      β_1[i] = μᵦ + sqrt(var_b)*randn()
  end

  # Sort
  sort!(β_1);

  β_2 = Array(Float64, length(β_1))
  for i in [1:length(β_2)]
      if i < length(β_2)/2
          β_2[i] = β_1[i] + sqrt(var_b)*randn()
      else
          β_2[i] = β_1[i] + sqrt(var_b)*randn()+0.015
      end
  end

  # Create one beta matrix that holds each agents beta for each year
  β = Array(Float64, agents*bs, tW)

  for t = 1:tW
      for i = 1:bs
          if t < br
              β[(i-1)*agents+1:agents*i,t] = β_1[i]
          else
              β[(i-1)*agents+1:agents*i,t] = β_2[i]
          end
      end
  end

  α = reshape(repmat(α,1,100)', agents*bs, 1)

  # Draw the income distribution:
  yit = zeros(bs*agents, tW)
  z = zeros(bs*agents, tW)

  for t = 1:tW
      for i = 1:bs*agents
          yit[i, t] = exp(α[i] + β[i, t]*t + z[t] + sqrt(var_ɛ)*randn())
          if t < tW
              z[i, t+1] = ρ*z[i, t] + sqrt(var_η)*randn()
          end
      end
  end
  @printf "\tβ=[%.2f, %.2f], β_2=[%.2f, %.2f]\n" minimum(β_1) maximum(β_1) minimum(β_2) maximum(β_2)

  # Calculate median income in last period for calculation of retirement benefits
  ymedian = median(yit[:, end])
  @printf "\tMedian income in period 40 is %.2f\n" ymedian

  pension = Array(Float64, agents*bs)
  for i = 1:agents*bs # Directly copied out of Guvenen's code
    ytemp = yit[i, end]
    if (ytemp<0.3*ymedian)
      yfixed=0.9*ytemp
    elseif (ytemp<=2.0*ymedian)
      yfixed=0.27*ymedian+0.32*(ytemp-0.3*ymedian)
    elseif (ytemp<4.1*ymedian)
      yfixed = 0.81*ymedian+0.15*(ytemp-2.0*ymedian)
    else
      yfixed = 1.1*ymedian
    end
    pension[i] = 0.715*yfixed
  end

  return yit, α, β, ymedian, pension
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
  for t = 1:tW
    for z = 1:zgridpoints_RIP
      for ɛ = 1:epsgridpoints
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
    end
  end

  return yit
end
