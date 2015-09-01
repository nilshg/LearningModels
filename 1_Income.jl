################################################################################
########################### INCOME DISTRIBUTION ################################
################################################################################

using Distributions

################################################################################

function incomeDistribution(user::AbstractString)

  @printf "1. Import Guvenen's income distribution\n"
  path="C:/Users/"*user*"/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/"
  yit = readdlm(path*"LaborReal.dat")
  println("\tMedian income in period 40 is $(median(yit[:, end]))")

  # Retirement benefits
  ybar_i = mean(yit, 2)[:]
  (γ_0, γ_1) = linreg(yit[:, 40], ybar_i)

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

  pension = Array(Float64, 100000)
  ybar = mean(yit)
  for i = 1:100000
    pension[i] = get_pension(yit[i, 40], γ_0, γ_1, ybar)
  end

  return yit, pension
end

################################################################################

function incomeDistribution{T<:AbstractFloat}(agents::Int64, bs::Int64, μₐ::T,
  μᵦ::T,var_α::T,var_β::T,cov_αβ::T,var_ɛ::T,var_η::T,ρ::T,y_adj::T,tW::Int64)

  min_β = max(-0.05, μᵦ-2.5*sqrt(var_β))

  @printf "1. Draw an income distribution\n"
  # Draw some alphas and betas
  ab = MvNormal([μₐ; μᵦ], [var_α cov_αβ; cov_αβ var_β])
  draw1 = rand(ab, bs)'
  draw2 = rand(ab, bs)'
  for i = 1:bs
    draw1[i,2] < min_β ? draw1[i,2] = min_β + abs(draw1[i,2]-min_β)/50.0 : 0
    draw2[i,2] < min_β ? draw2[i,2] = min_β + abs(draw2[i,2]-min_β)/50.0 : 0
  end
  α = (draw1[:,1] + draw2[:,1])/2.
  α = reshape(repmat(α,1,100)', agents*bs)
  β_u = reshape(repmat(draw1[:, 2],1,100)', agents*bs)
  β_k = reshape(repmat(draw2[:, 2],1,100)', agents*bs)
  β = (1-fpu)*β_k + fpu*β_u

  # Draw the income distribution:
  yit = zeros(bs*agents, tW); z = similar(yit)

  z[:, 1] = sqrt(var_η/(1-ρ^2))*randn(agents*bs)
  for t = 1:tW, i = 1:bs*agents
    yit[i, t] = y_adj + exp(α[i] + β[i]*t + z[i, t] + sqrt(var_ɛ)*randn())
    t < tW ? z[i, t+1] = ρ*z[i, t] + sqrt(var_η)*randn() : 0
  end

  # Median income in last period for calculation of retirement benefits
  ymedian =
  println("\tMedian income in period 40 is $(median(yit[:, end]))")

  ybar_i = mean(yit, 2)[:]
  (γ_0, γ_1) = linreg(yit[:, 40], ybar_i)

  function get_pension(y::T, k_0::T, k_1::T, avgy::T)
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
  ybar = mean(yit)
  for i = 1:size(yit, 1)
    pension[i] = get_pension(yit[i, 40], γ_0, γ_1, ybar)
  end

  return yit, pension, α, β, β_k
end

################################################################################

function incomeDistribution{T<:AbstractFloat}(α::T, β::T, var_η_RIP::T,
  var_ɛ_RIP::T, zpoints_RIP::Int64, epspoints::Int64, tW::Int64)

  @printf "1. Drawing an RIP-consistent income distribution\n"
  println("\tσ²(η) = $var_η_RIP,  σ²(ɛ) =  $var_ɛ_RIP")

  zdisc = [-(1/2)*(zpoints_RIP-1)*sqrt(var_η_RIP)+(i-1)*sqrt(var_η_RIP)
               for i = 1:zpoints_RIP]

  epsdisc = [-(1/2)*(epspoints-1)*sqrt(var_ɛ_RIP)+(i-1)*sqrt(var_ɛ_RIP)
               for i = 1:epspoints]

  yit = Array(Float64, (zpoints_RIP, epspoints, tW))

  for t = 1:tW, z = 1:zpoints_RIP, ɛ = 1:epspoints
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
