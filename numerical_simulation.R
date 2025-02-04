# ====================================================
# Numerical Simulations of Predator-Prey Dynamics with NTIs
# ====================================================

# Load necessary libraries
library(deSolve)         # Solving differential equations
library(tidyverse)       # Data manipulation and visualization
library(plyr)            # List-based operations
library(stringr)         # String manipulation (if needed)
# ====================================================
# 1. DEFINE PREDATOR-PREY MODEL WITH NON-TROPHIC INTERACTIONS (NTIs)
# ====================================================

# Function for calculating the rate of change of prey (dV) and predator (dP)
PVmod_nti <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    # Prey population growth
    dV <- j*Lv*(A - Wv*V) - q*V - ((lambda*V)/(1+lambda*h*V)) * P  # Prey growth and predation impact
    
    # Predator population growth
    dP <- ((g + (g_prime - g) * (V^a / (V^a + Vor^a))) * Lp * (A - Wp * P)) - 
      ((mmax - (mmax - mi) * (V^d / (V^d + Vos^d))) * P)  # Recruitment facilitation & refuge provision
    
    # Return rate of change as a list
    list(c(dV, dP))
  })
}

# Vectorized version for large-scale simulations
PVmod_nti_vectorized <- function(Time, State, Pars) { 
  with(as.list(c(Pars)), {
    
    # Separate prey and predator populations
    V <- State[seq(1, length(State)/2)]
    P <- State[seq(1+length(State)/2, length(State))]
    
    # Prey growth and predator effect
    dV <- j*Lv*(A - Wv*V) - q*V - ((lambda*V)/(1+lambda*h*V)) * P
    
    # Predator population change with NTI effects
    dP <- ((g + (g_prime - g) * (V^a / (V^a + Vor^a))) * Lp * (A - Wp * P)) - 
      ((mmax - (mmax - mi) * (V^d / (V^d + Vos^d))) * P)
    
    # Return updated state
    list(c(dV, dP))
  })
}

# ====================================================
# 2. INITIAL CONDITIONS & TIMEFRAME
# ====================================================

# Initial prey (V) and predator (P) abundances
yiniPV <- c(V = 3000, P = 8.217147) 

# Define simulation timeframe
times <- seq(0, 1000, by = 0.1)

# ====================================================
# 3. DEFINE PARAMETER SPACE FOR SIMULATIONS
# ====================================================

# Generate all possible parameter combinations
simulation_plan <- expand.grid(
  Lv      = seq(0, 5000, by=50),      # Prey recruitment rate
  Lp      = c(seq(0, 40, by=2), 200, 2000),  # Predator larval availability
  mi      = c(0.005, 0.155, 0.255),   # Minimum predator mortality
  mmax    = c(0.005, 0.155, 0.255),   # Maximum predator mortality
  g       = c(0.1),                   # Baseline predator recruitment rate
  g_prime = c(0.1, 0.25, 0.5)          # Enhanced predator recruitment with facilitation
)

# Remove biologically unrealistic parameter combinations
simulation_plan <- simulation_plan %>% 
  mutate(delta_g = g_prime - g, delta_m_positive = mmax - mi) %>%
  filter(delta_m_positive >= 0 & delta_g >= 0)

# ====================================================
# 4. ADD FIXED MODEL PARAMETERS TO THE SIMULATION PLAN
# ====================================================

# Default model parameters
parameters <- list(
  a = 2,         # Recruitment facilitation exponent
  d = 10,        # Refuge provision exponent
  j = 0.5,       # Prey intrinsic growth rate
  A = 1,         # Available space
  Wv = 0.0001,   # Space occupied by prey
  Wp = 0.022,    # Space occupied by predators
  q = 0.005,     # Prey mortality rate
  lambda = 0.923725, # Predation coefficient
  h = 0.030306,  # Handling time
  Vor = 3000,    # Threshold for recruitment facilitation
  Vos = 7000     # Threshold for refuge provision
)

# Append fixed parameters to all rows in the simulation plan
for (n in names(parameters)) {
  if (!n %in% names(simulation_plan)) {
    simulation_plan[, n] <- parameters[[n]]
  }
}

# ====================================================
# 5. RUN SIMULATIONS USING ODE SOLVER
# ====================================================

# Initial state vector, repeated for each parameter combination
yiniPV_vectorized <- c(rep(yiniPV, each = nrow(simulation_plan)))

# Choose between full or chunked simulation execution
RUN_BY_CHUNK <- TRUE
CHUNK_SIZE <- 256  # Chunk size should be even

if (!RUN_BY_CHUNK) { 
  # Run full simulation in a single step
  yiniPV_vectorized <- c(rep(yiniPV, each = nrow(simulation_plan)))
  timings <- system.time({
    out <- ode(y = yiniPV_vectorized, times = times, func = PVmod_nti_vectorized, parms = simulation_plan)
  })
  var_si <- c(out[, -1])
} else { 
  # Split simulations into chunks for efficiency
  splits <- seq.int(nrow(simulation_plan)) %/% CHUNK_SIZE 
  splits <- llply(unique(splits), function(spl) which(splits == spl))
  
  # Run simulations in chunks (parallelization supported)
  results <- llply(splits, function(sub_sp) { 
    yiniPV_vectorized <- c(rep(yiniPV, each = length(sub_sp)))
    out <- ode(y = yiniPV_vectorized, times = times, func = PVmod_nti_vectorized, parms = simulation_plan[sub_sp, ])
    
    # Store results in two-column format (prey, predator)
    out <- c(out[, -1])
    cbind(head(out, length(out)/2), tail(out, length(out)/2))
  }, .progress = "time")  # Progress bar enabled
  
  var_si <- c(do.call(rbind, results))
}

# ====================================================
# 6. ORGANIZE OUTPUT DATA
# ====================================================

# Create dataframe with simulation results
temp_traj <- data.frame(
  simulation_plan[rep(seq_len(nrow(simulation_plan)), each = length(times)), ],
  time = rep(times, times = length(simulation_plan[, 1])),
  Prey = head(var_si, length(var_si) / 2),
  Predator = tail(var_si, length(var_si) / 2)
)

# Assign interaction types based on parameter values
temp_traj <- temp_traj %>%
  mutate(interaction_type = case_when(
    delta_m_positive == 0 & delta_g == 0 ~ 'Only trophic',
    delta_m_positive != 0 & delta_g == 0 ~ 'Refuge provision',
    delta_m_positive != 0 & delta_g != 0 ~ 'Both NTI',
    TRUE ~ 'Recruitment facilitation'
  ))

# Extract steady-state values (final time point)
last_time = filter(temp_traj, time == max(unique(temp_traj$time)))

# ====================================================
# 7. SAVE OUTPUT FILES
# ====================================================

# Save full time-series and final equilibrium values
# write.csv(temp_traj, file="temporal_trajectories.csv", row.names=FALSE)
# write.csv(last_time, "last_time.csv")
