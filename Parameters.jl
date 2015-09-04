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
γ = 0.8                # Strenght of habits (0 = no habits)
λ = 0.8                # Persistence of habits (lambda=0 --> no persistence)
σ = 2.0                # Coefficient of relative risk aversion

# Paremeters for belief calculation
fpu = 0.65             # Proportion of variance of beta that is known
init_α = 2.0           # Initial belief about mean of beta
init_β = 0.009         # Initial belief about mean of beta
init_z = 0.0           # Initial belief about mean of z
init_var_β = 0.01      # Initial belief about variance of beta
init_var_z = 0.01      # Initial belief about ariance of z

# Parameters of income process (HIP)
μₐ = 2.0               # Mean of the income intercept
μᵦ = 0.009              # Mean of income growth rates
var_α = 0.005          # Cross-sectional variance of α (std = 0.07)
var_β = 0.00037        # Cross-sectional variance of β (std = 0.019)
corr_αβ = -0.25        # Cross-sectional correlation of α, β (cov = -0.0003)
cov_αβ = corr_αβ*sqrt(var_β*var_α)
var_η = 0.029          # Variance of persistent shock (std = 0.17)
var_ɛ = 0.047          # Variance of transitory shock (std = 0.22)
ρ = 0.82               # Persistence of AR(1) component of income
br = 40                # Period of structural break (has to be <T)
y_adj = 0.4            # For comparability with Guvenen's code

# Parameters of income process (RIP)
var_η_RIP = 0.015      # σ²(η) (std = 0.122)
var_ɛ_RIP = 0.061      # σ²(ɛ) (std = 0.247)
ρ_RIP = 0.988          # AR(1) persistence

# Parameters for grid construction
wpoints = 25           # Wealth grid points (working life)
wpoints_R = 160        # Wealth grid points (retirement)
hpoints = 6            # Habit grid points (working life)
hpoints_R = 35         # Habit grid points (retirement)
ypoints_R = 150        # Pension points
wmaxR = 1000.0         # Maximum retirement wealth
apoints = 3            # Grid points for beliefs about α
bpoints = 11           # Grid points for beliefs about β
zpoints = 8            # Grid points for beliefs about z
zpoints_RIP = 32       # Grid points for RIP persistent shock
epspoints = 2          # Grid points for RIP transitory shock
power = 2.0            # Wealth Grid Curvature

# Utility function
function u_h(c::Float64, h::Float64, γ=γ, σ=σ)
    ((c/(h^γ))^(1-σ))/(1-σ)
end

function u(c::Float64, σ=σ)
     c^(1-σ)/(1-σ)
end
