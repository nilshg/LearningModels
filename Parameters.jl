####################################################################################
################################# PARAMETERS #######################################
####################################################################################

tW = 40                # Number of working life periods
tR = 30                # Number of retirement periods
r = 1/0.96             # Real Return
agents = 100           # Number of agents sharing an alpha/beta combination
bs = 1000              # Different alpha/beta combinations in the population;

# Preference Parameters
δ = 0.966              # Discount factor (Note: 1/beta = 1.035)
γ = 0.8                # Strenght of habits (0 = no habits, 1 = only relative consumption matters)
λ = 0.8                # Persistence of habits (lambda=0 --> no persistence)
σ = 2.0                # Coefficient of relative risk aversion

# Parameters of income process (HIP)
μₐ = 2.0               # Mean of the income intercept
μᵦ = 0.009              # Mean of income growth rates
var_α = 0.005          # Cross-sectional variance of alpha
var_β = 0.00037        # Cross-sectional variance of beta
var_η = 0.029          # Variance of persistent shock
var_ɛ = 0.047          # Variance of transitory shock
ρ = 0.82               # Persistence of AR(1) component of income
br = 30                # Period of structural break (has to be <T)

# Parameters of income process (RIP)
var_η_RIP = 0.015      # σ²(η)
var_ɛ_RIP = 0.061      # σ²(ɛ)
ρ_RIP = 0.988          # AR(1) persistence

# Paremeters for belief calculation
fpu = 0.65             # Proportion of variance of beta that is known
init_β = 0.029         # Initial belief about mean of beta
init_z = 0.0           # Initial belief about mean of z
init_var_β = 0.01      # Initial belief about variance of beta
init_var_z = 0.01      # Initial belief about ariance of z

# Parameters for grid construction
wpoints = 40           # Wealth grid points (working life)
wpoints_R = 130        # Wealth grid points (retirement)
hpoints = 7            # Habit grid points (working life)
hpoints_R = 29         # Habit grid points (retirement)
ypoints_R = 120        # Pension points
wmaxR = 1000.0         # Maximum retirement wealth
apoints = 4            # Grid points for beliefs about α
bpoints = 8           # Grid points for beliefs about β
zpoints = 7            # Grid points for beliefs about z
zpoints_RIP = 32       # Grid points for RIP persistent shock
epspoints = 2          # Grid points for RIP transitory shock
power = 1.0            # Grid Curvature

# Utility function
function u(c::Float64, h::Float64, γ=γ, σ=σ)
    ((c/(h^γ))^(1-σ))/(1-σ)
end

function u(c::Float64, σ=σ)
     (c^(1-σ))/(1-σ)
end
