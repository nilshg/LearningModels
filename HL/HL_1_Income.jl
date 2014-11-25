################################################################################
########################### INCOME DISTRIBUTION ################################
################################################################################

function IncomeDistribution{T<:Int64}(agents::T, bs::T)

  @printf "1. Draw an Income Distribution\n"
  @printf "\tWe are working with Guvenen's data\n"
  Yit = readcsv("C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/LaborReal.csv")
  alfabeta = readcsv("C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/alfabeta.csv")
  alpha = alfabeta[:, 1]
  alpha = reshape(repmat(alpha, 1, agents)', agents*bs, 1)
  beta = alfabeta[:, 2]
  beta = reshape(repmat(beta, 1, agents)', agents*bs, 1)

  # Calculate median income in last period for calculation of retirement benefits
  ymedian = median(Yit[:, 40])
  @printf "\tMedian income in period 40 is %.2f\n" ymedian

  return Yit, alpha, beta, ymedian
end

function IncomeDistribution{T<:Float64}(agents::Int64, bs::Int64, μₐ::T, μᵦ::T,
                            var_a::T, var_b::T, var_eps::T,
                            ρ::T, var_eta::T, BR::Int64, TW::Int64)

  # Draw some alphas and betas
  alpha = Array(Float64, bs)
  beta1 = Array(Float64, bs)
  for i = 1:bs
      alpha[i] = μₐ + sqrt(var_a)*randn()
      beta1[i] = μᵦ + sqrt(var_b)*randn()
  end

  # Sort
  sort!(beta1);

  beta2 = Array(Float64, (length(beta1), 1))
  for i in [1:length(beta2)]
      if i < length(beta2)/2
          beta2[i] = beta1[i] + sqrt(var_b)*randn()
      else
          beta2[i] = beta1[i] + sqrt(var_b)*randn()+0.03
      end
  end

  # Create one beta matrix that holds each agents beta for each year
  beta = zeros(agents*bs, TW)

  for t = 1:TW
      for i = 1:bs
          if t < BR
              beta[(i-1)*agents+1:agents*i,t] = beta1[i]
          else
              beta[(i-1)*agents+1:agents*i,t] = beta2[i]
          end
      end
  end

  # Reshape alpha into 100,000x1 vector with each of the 1,000 unique values repeated 100 times
  alpha = reshape(repmat(alpha,1,100)', agents*bs, 1)

  # Draw the income distribution:
  Yit = zeros(bs*agents, TW)
  z = zeros(bs*agents, TW)

  for t = 1:TW
      for i = 1:bs*agents
          Yit[i, t] = exp(alpha[i] + beta[i, t]*t + z[t] + sqrt(var_eps)*randn())
          if t < TW
              z[i, t+1] = ρ*z[i, t] + sqrt(var_eta)*randn()
          end
      end
  end
  @printf "\tβ is between %.2f and %.2f, β_2 is between %.2f and %.2f\n" minimum(beta1)  maximum(beta1) minimum(beta2) maximum(beta2)

  # Calculate median income in last period for calculation of retirement benefits
  ymedian = median(Yit[:, 40])
  @printf "\tMedian income in period 40 is %.2f\n" ymedian

  return Yit, alpha, beta, ymedian
end
