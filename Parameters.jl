####################################################################################
################################# PARAMETERS #######################################
####################################################################################

tol = 1e-9             # Tolerance Level for Iterations
TW = 40                # Number of working life periods
TR = 30                # Number of retirement periods
r = 0.015              # Real interest rate
R = 1 + r              # Handy for notation
agents = 100           # Number of agents sharing an alpha/beta combination
bs = 1000              # Different alpha/beta combinations in the population;

# Preference Parameters
δ = 0.966              # Discount factor (Note: 1/beta = 1.035)
γ = 0.8                # Strenght of habits (0 = no habits, 1 = only relative consumption matters)
λ = 0.8                # Persistence of habits (lambda=0 --> no persistence)
σ = 2.0                # Coefficient of relative risk aversion

# Parameters of income distribution
μₐ = 2.0               # Mean of the income intercept
μᵦ = 0.005             # Mean of income growth rates
var_a = 0.005          # Cross-sectional variance of alpha
var_b = 0.00037        # Cross-sectional variance of beta
var_eta = 0.029        # Variance of persistent shock
var_eps = 0.047        # Variance of transitory shock
ρ = 0.82               # Persistence of AR(1) component of income
BR = 30                # Period of structural break (has to be <T)

# Paremeters for belief calculation
FPU = 0.65             # Proportion of variance of beta that is known
init_beta = 0.04       # Initial belief about mean of beta
init_z = 0.0           # Initial belief about mean of z
init_varbeta = 0.01    # Initial belief about variance of beta
init_varz = 0.01       # Initial belief about ariance of z

# Parameters for grid construction
Sgridpoints = 15       # Partitions for belief cube
wgridpoints = 10       # Wealth grid points (working life)
wgridpoints_R = 60     # Wealth grid points (retirement)v
hgridpoints = 7        # Habit grid points (working life)
hgridpoints_R = 40     # Habit grid points (retirement)
ygridpoints = 6        # Income gridpoints (working life)
ygridpoints_R = 50     # Pension gridpoints
wmaxR = 1000.0         # Maximum retirement wealth
agridpoints = 4        # Grid points for beliefs about α
bgridpoints = 8        # Grid points for beliefs about β
zgridpoints= 7         # Grid points for beliefs about z
power = 1.0            # Grid Curvature
retirement_grid_irregular = false

# Utility function
function u{T<:Float64}(c::T, h::T, γ=γ, σ=σ)
    ((c/(h^γ))^(1-σ))/(1-σ)
end

function u{T<:Float64}(c::T, σ=σ)
     (c^(1-σ))/(1-σ)
end
