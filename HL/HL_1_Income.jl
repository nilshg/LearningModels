################################################################################
############################ INCOME DISTRIBUTION ###############################
################################################################################
@printf "2. Draw an Income Distribution\n"

guvenen_distribution = true

if guvenen_distribution # Use distribution from Guvenen's paper
    @printf "\tWe are working with Guvenen's data\n"
    Yit = readcsv("C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/LaborReal.csv")
    alfabeta = readcsv("C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/alfabeta.csv")
    alpha = alfabeta[:, 1]
    alpha = reshape(repmat(alpha, 1, agents)', agents*bs, 1)
    beta = alfabeta[:, 2]
    beta = reshape(repmat(beta, 1, agents)', agents*bs, 1)
else # Draw a new distribution
    alpha = zeros(bs)
    beta1 = zeros(bs)
    for i = 1:bs
        alpha[i] = mu_a + sqrt(var_a)*randn()
        beta1[i] = mu_b + sqrt(var_b)*randn()
    end

    # Sort
    sort!(beta1);

    beta2 = zeros(length(beta1))
    for i in [1:length(beta2)]
        if i < length(beta2)/2
            beta2[i] = beta1[i] + sqrt(var_b)*randn()
        else
            beta2[i] = beta1[i] + sqrt(var_b)*randn()+0.03
        end
    end
    # Set break point at which growth rates switches from beta to beta2
    BR = 30;

    # Create one beta matrix that holds each agents beta for each year
    beta = zeros(agents*bs, T)

    for t = 1:T
        for i = 1:bs
            if t < BR
                beta[(i-1)*agents+1:agents*i,t] = beta1[i]
            else
                beta[(i-1)*agents+1:agents*i,t] = beta2[i]
            end
        end
    end

    alpha = reshape(repmat(alpha,1,100)', agents*bs, 1)

    # Draw the income distribution:
    Yit = zeros(bs*agents, T)
    z = zeros(bs*agents, T)

    for t = 1:T
        for i = 1:bs*agents
            Yit[i, t] = exp(alpha[i] + beta[i, t]*t + z[t] + sqrt(var_eps)*randn())
            if t < T
                z[i, t+1] = ar1rho*z[i, t] + sqrt(var_eta)*randn()
            end
        end
    end
    @printf "\tβ is between %.2f and %.2f, β_2 is between %.2f and %.2f\n" minimum(beta1)  maximum(beta1) minimum(beta2) maximum(beta2)

end;

# Calculate median income in last period for calculation of retirement benefits
ymedian = median(Yit[:, 40])
@printf "\tMedian income in period 40 is %.2f\n" ymedian
