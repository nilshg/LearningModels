####################################################################################
################################# PARAMETERS #######################################
####################################################################################

tW = 40                # Number of working life periods
tR = 30                # Number of retirement periods
r = 1/0.96             # Real Return
agents = 100           # Number of agents sharing an alpha/beta combination
bs = 1000              # Different alpha/beta combinations in the population;

# Preference Parameters
δ = 0.955              # Discount factor (Note: 1/beta = 1.035)
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
br = 40                # Period of structural break (has to be <T)

# Parameters of income process (RIP)
var_η_RIP = 0.015      # σ²(η) (std = 0.122)
var_ɛ_RIP = 0.061      # σ²(ɛ) (std = 0.247)
ρ_RIP = 0.988          # AR(1) persistence
ξ = 0.00               # Probability of zero income shock

# Parameters for grid construction
xpoints = 40           # Cash-in-hand grid points (working life)
wpoints_R = 90         # Wealth grid points (retirement)
hpoints = 6            # Habit grid points (working life)
hpoints_R = 15         # Habit grid points (retirement)
ypoints_R = 80         # Pension points
wmaxR = 1000.0         # Maximum retirement wealth
apoints = 3            # Grid points for beliefs about α
bpoints = 11           # Grid points for beliefs about β
zpoints = 7            # Grid points for beliefs about z
zpoints_RIP = 32       # Grid points for RIP persistent shock
epspoints = 2          # Grid points for RIP transitory shock
power = 3.5            # Wealth Grid Curvature

# Utility functions
function u(c::Float64, h::Float64, γ::Float64, σ::Float64)
    ((c/(h^γ))^(1-σ))/(1-σ)
end

function u(c::Float64, σ::Float64)
  if c > 0.00001
     ut = c^(1-σ)/(1-σ)
   else
     ut = -10000. - 100*c^2.
   end
end

# Bequest motive
υ = (-1/0.544)
function bq(w::Float64, υ::Float64)
  υ==0.0 ? 0.0 :57726.^υ*w
end

# Probability of death
 ψ = hcat([.01368 .01493 .01628 .01767 .01911 .02059 .02216 .02389 .02585],
          [.02806 .03052 .03315 .03593 .03882 .04184 .04507 .04867  .05274],
          [.05742 .06277 .06882 .07552 .08278 .09041 .09842 .10725 .11712],
          [.12717 .13708 .147281])
