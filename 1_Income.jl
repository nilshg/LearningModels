################################################################################
########################### INCOME DISTRIBUTION ################################
################################################################################

using Distributions

################################################################################

function incomeDistribution(user::AbstractString)

  @printf "Import Guvenen's income distribution\n"
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

function incomeDistribution{T<:Int}(agents::T, bs::T, tW::Int64; profile="none")

  # Life cycle profiles
  if profile == "psid"
    # PSID log real labour income, 1968-1996 (median income at 40: $29703)
    g_t = [8.36 + 0.07*t - 0.15*t^2/100 for t = 1:tW]
    μₐ = 1.17; μᵦ = 0.009
    var_α = 0.03; var_β = 0.00031; corr_αβ = -0.3
    cov_αβ = corr_αβ*sqrt(var_β*var_α)
    var_η = 0.013; var_ɛ = 0.03
    ρ = 0.853
  elseif profile == "psid_68_86"
    # PSID log real labour income, 1968-1996 (median income at 40: $29703)
    g_t = [8.36 + 0.07*t - 0.15*t^2/100 for t = 1:tW]
    μₐ = 1.17; μᵦ = 0.009
    var_α = 0.11; var_β = 0.00001; corr_αβ = -0.42
    cov_αβ = corr_αβ*sqrt(var_β*var_α)
    var_η = 0.013; var_ɛ = 0.043
    ρ = 0.885
  elseif profile == "psid_87_13"
    # PSID log real labour income, 1968-1996 (median income at 40: $29703)
    g_t = [8.36 + 0.07*t - 0.15*t^2/100 for t = 1:tW]
    μₐ = 1.17; μᵦ = 0.009
    var_α = 0.097; var_β = 0.00025; corr_αβ = -0.31
    cov_αβ = corr_αβ*sqrt(var_β*var_α)
    var_η = 0.032; var_ɛ = 0.085
    ρ = 0.854
  elseif profile == "bhps_grosslab"
    # BHPS log real gross labour income, 1992-2008 (median £25k mean £28.7k)
    g_t = [9.65 + 0.024*t + 0.019*t^2/100 - 0.014*t^3/1000 for t = 1:tW]
    μₐ = .0; μᵦ = 0.009
    var_α = 0.036; var_β = 0.00032; corr_αβ = -0.51
    cov_αβ = corr_αβ*sqrt(var_β*var_α)
    var_η = 0.106; var_ɛ = 0.08
    ρ = 0.719
  elseif profile == "bhps_netlab"
    # BHPS log real net labour income, 1992-2008 (median £19k, mean £21.6k)
    g_t = [9.41 + 0.019*t + 0.029*t^2/100 - 0.014*t^3/1000 for t = 1:tW]
    μₐ = .0;  μᵦ = 0.009
    var_α = 0.032; var_β = 0.00019; corr_αβ = -0.59
    cov_αβ = corr_αβ*sqrt(var_β*var_α)
    var_η = 0.073; var_ɛ = 0.07
    ρ = 0.808
  elseif profile == "bhps_netinc"
    # BHPS log real net houseold income, 1992-2008 (median £22.6k, mean £26k)
    g_t = [9.49 + 0.026*t - 0.008*t^2/100 - 0.007*t^3/1000 for t = 1:tW]
    μₐ = .0; μᵦ = 0.009
    var_α = 0.052; var_β = 0.00011; corr_αβ = -0.42
    cov_αβ = corr_αβ*sqrt(var_β*var_α);
    var_η = 0.027; var_ɛ = 0.056
    ρ = 0.857
  elseif profile == "bhps_netincdef"
    # BHP log real net household income, deflated & equivalized, 1992-2008
    # (median £16.9k, mean £19.4k)
    g_t = [9.64 - 0.002*t + 0.027*t^2/100 - 0.004*t^3/1000 for t = 1:tW]
    μₐ = .0; μᵦ = 0.009
    var_α = 0.1; var_β = 0.0; corr_αβ = -1.
    cov_αβ = corr_αβ*sqrt(var_β*var_α)
    var_η = 0.039; var_ɛ = 0.042
    ρ = 0.812
  elseif profile == "baseline"
    # PSID Profile
    g_t = [9.64 - 0.002*t + 0.027*t^2/100 - 0.004*t^3/1000 for t = 1:tW]
    μₐ = .0; μᵦ = 0.009
    var_α = 0.03; var_β = 0.00038; corr_αβ = -0.2
    cov_αβ = corr_αβ*sqrt(var_β*var_α);
    var_η = 0.03; var_ɛ = 0.05;
    ρ = 0.82
  elseif profile == "low_beta"
    g_t = [9.64 - 0.002*t + 0.027*t^2/100 - 0.004*t^3/1000 for t = 1:tW]
    μₐ = .0; μᵦ = 0.009;
    var_α = 0.03; var_β = 0.00015; corr_αβ = -0.2
    cov_αβ = corr_αβ*sqrt(var_β*var_α);
    var_η = 0.03; var_ɛ = 0.05;
    ρ = 0.82
  elseif profile == "high_beta"
    g_t = [9.64 - 0.002*t + 0.027*t^2/100 - 0.004*t^3/1000 for t = 1:tW]
    μₐ = .0; μᵦ = 0.009;
    var_α = 0.03; var_β = 0.0007; corr_αβ = -0.2
    cov_αβ = corr_αβ*sqrt(var_β*var_α);
    var_η = 0.03; var_ɛ = 0.05;
    ρ = 0.82
  elseif profile == "low_alpha"
    g_t = [9.64 - 0.002*t + 0.027*t^2/100 - 0.004*t^3/1000 for t = 1:tW]
    μₐ = .0; μᵦ = 0.009;
    var_α = 0.01; var_β = 0.00038; corr_αβ = -0.2
    cov_αβ = corr_αβ*sqrt(var_β*var_α);
    var_η = 0.03; var_ɛ = 0.05;
    ρ = 0.82
  elseif profile == "high_alpha"
    g_t = [9.64 - 0.002*t + 0.027*t^2/100 - 0.004*t^3/1000 for t = 1:tW]
    μₐ = .0; μᵦ = 0.009;
    var_α = 0.1; var_β = 0.00038; corr_αβ = -0.2
    cov_αβ = corr_αβ*sqrt(var_β*var_α);
    var_η = 0.03; var_ɛ = 0.05;
    ρ = 0.82
  elseif profile == "low_cov"
    # PSID Profile
    g_t = [9.64 - 0.002*t + 0.027*t^2/100 - 0.004*t^3/1000 for t = 1:tW]
    μₐ = .0; μᵦ = 0.009
    var_α = 0.03; var_β = 0.00038; corr_αβ = -0.8
    cov_αβ = corr_αβ*sqrt(var_β*var_α);
    var_η = 0.03; var_ɛ = 0.05;
    ρ = 0.82
  elseif profile == "high_cov"
    # PSID Profile
    g_t = [9.64 - 0.002*t + 0.027*t^2/100 - 0.004*t^3/1000 for t = 1:tW]
    μₐ = .0; μᵦ = 0.009
    var_α = 0.03; var_β = 0.00038; corr_αβ = 0.
    cov_αβ = corr_αβ*sqrt(var_β*var_α);
    var_η = 0.03; var_ɛ = 0.05;
    ρ = 0.82
  elseif profile == "low_eta"
    g_t = [9.64 - 0.002*t + 0.027*t^2/100 - 0.004*t^3/1000 for t = 1:tW]
    μₐ = .0; μᵦ = 0.009;
    var_α = 0.03; var_β = 0.00038; corr_αβ = -0.2
    cov_αβ = corr_αβ*sqrt(var_β*var_α);
    var_η = 0.01; var_ɛ = 0.05;
    ρ = 0.82
  elseif profile == "high_eta"
    g_t = [9.64 - 0.002*t + 0.027*t^2/100 - 0.004*t^3/1000 for t = 1:tW]
    μₐ = .0; μᵦ = 0.009;
    var_α = 0.03; var_β = 0.00038; corr_αβ = -0.2
    cov_αβ = corr_αβ*sqrt(var_β*var_α);
    var_η = 0.1; var_ɛ = 0.05;
    ρ = 0.82
  elseif profile == "low_eps"
    g_t = [9.64 - 0.002*t + 0.027*t^2/100 - 0.004*t^3/1000 for t = 1:tW]
    μₐ = .0; μᵦ = 0.009;
    var_α = 0.03; var_β = 0.00038; corr_αβ = -0.2
    cov_αβ = corr_αβ*sqrt(var_β*var_α);
    var_η = 0.03; var_ɛ = 0.02;
    ρ = 0.82
  elseif profile == "high_eps"
    g_t = [9.64 - 0.002*t + 0.027*t^2/100 - 0.004*t^3/1000 for t = 1:tW]
    μₐ = .0; μᵦ = 0.009;
    var_α = 0.03; var_β = 0.00038; corr_αβ = -0.2
    cov_αβ = corr_αβ*sqrt(var_β*var_α);
    var_η = 0.03; var_ɛ = 0.13;
    ρ = 0.82
  elseif profile == "low_rho"
    g_t = [9.64 - 0.002*t + 0.027*t^2/100 - 0.004*t^3/1000 for t = 1:tW]
    μₐ = .0; μᵦ = 0.009;
    var_α = 0.03; var_β = 0.00038; corr_αβ = -0.2
    cov_αβ = corr_αβ*sqrt(var_β*var_α);
    var_η = 0.03; var_ɛ = 0.05;
    ρ = 0.71
  elseif profile == "high_rho"
    g_t = [9.64 - 0.002*t + 0.027*t^2/100 - 0.004*t^3/1000 for t = 1:tW]
    μₐ = .0; μᵦ = 0.009;
    var_α = 0.03; var_β = 0.00038; corr_αβ = -0.2
    cov_αβ = corr_αβ*sqrt(var_β*var_α);
    var_η = 0.03; var_ɛ = 0.05;
    ρ = 0.91
  elseif profile == "RIP"
    g_t = [9.64 - 0.002*t + 0.027*t^2/100 - 0.004*t^3/1000 for t = 1:tW]
    μₐ = .0; μᵦ = 0.009;
    var_α = 0.0; var_β = 0.0; corr_αβ = 0.
    cov_αβ = corr_αβ*sqrt(var_β*var_α);
    var_η = 0.05; var_ɛ = 0.03;
    ρ = 0.95
  elseif profile == "guvenen"
    g_t = [0.0 for t = 1:tW]
    μₐ = 2.0; μᵦ = 0.009;
    var_α = 0.005; var_β = 0.00038; corr_αβ = -0.23
    cov_αβ = corr_αβ*sqrt(var_β*var_α);
    var_η = 0.029; var_ɛ = 0.047;
    ρ = 0.82
  else
    g_t = [0. for t = 1:tW]
    μₐ = 1.17; μᵦ = 0.009
    var_α = 0.03; var_β = 0.00031; corr_αβ = -0.3
    cov_αβ = corr_αβ*sqrt(var_β*var_α)
    var_η = 0.013; var_ɛ = 0.03
    ρ = 0.853
    profile = "no"
  end

  println("Draw an income distribution using the $profile profile")
  println("Parameters of the income process are:")
  println("ρ=$ρ, var_α=$var_α, var_β=$var_β, var_η=$var_η, var_ɛ=$var_ɛ")
  # Draw some alphas and betas
  if (var_α>0.) & (var_β > 0.)
    min_β = max(-0.05, μᵦ-2.5*sqrt(var_β))
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
  else
    α = [μₐ for i = 1:agents*bs]
    β = [μᵦ for i = 1:agents*bs]
    β_k = [μᵦ for i = 1:agents*bs]
  end

  # Draw the income distribution:
  yit = zeros(bs*agents, tW); z = similar(yit)

  z[:, 1] = sqrt(var_η/(1-ρ^2))*randn(agents*bs)
  for t = 1:tW, i = 1:bs*agents
    yit[i, t] = exp(g_t[t] + α[i] + β[i]*t + z[i, t] + sqrt(var_ɛ)*randn()) +0.4
    t < tW ? z[i, t+1] = ρ*z[i, t] + sqrt(var_η)*randn() : 0
  end

  # Median income in last period for calculation of retirement benefitss
  println("\tMedian and mean income in period 40 are "*
                    "$(round(median(yit[:, end]),2)) and "*
                    "$(round(mean(yit[:, end]),2))")

  ybar_i = mean(yit, 2)[:]
  (γ_0, γ_1) = linreg(yit[:, 40], ybar_i)

  function get_pension{T<:AbstractFloat}(y::T, k_0::T, k_1::T, avgy::T)
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

  return yit, pension, α, β, β_k, g_t, ρ, var_α, var_β, cov_αβ, var_η, var_ɛ
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
