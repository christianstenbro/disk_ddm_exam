# ===== SET UP ======

# installing packages
install.packages("pacman")
pacman::p_load(tidyverse,
               brms,
               RWiener)

# load data <- change to ucloud path!
model_data_downsampled <- read_csv("/work/DATA/model_data_downsampled.csv")

# loading previous m4_fit
m4_fit <- readRDS("/work/fit_wiener_model_4.rds")

# get seed/iteration from command line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  seed_entropy <- 0
} else {
  seed_entropy <- as.numeric(args[1])
}

# set the same seed across bash iterations for generating payoff matrix
set.seed(1983)

# setting working directory - make ucloud friendly
setwd('/work')

# ==== DEFINING MODEL =====

# formula
m4_formula <- bf(RT | dec(response) ~ abs_H_difference * H_sum + (1|ID) + (1 | stimulus), 
                 bs ~ 1 + (1|ID), 
                 ndt ~ 1 + (1|ID), 
                 bias ~ 1 + (1|ID))
# priors
m4_prior <- c(
  prior(normal(0, 1), class = "b", coef = "abs_H_difference"),
  prior(normal(0, 1), class = "b", coef = "H_sum"),
  prior(normal(1.5, 0.5), class = "Intercept", dpar = "bs"),
  prior(normal(-3, 0.5), class = "Intercept", dpar = "ndt")
)

# init function
m4_initfun <- function() {
  list(
    b = as.array(c(0, 0, 0)), 
    Intercept = 0,
    Intercept_bs = 1.5,
    Intercept_ndt = -3.5, # this is in log!
    Intercept_bias = 0.5,
    # random effect SDs (for ID and stimulus)
    sd_1 = as.array(0.1), 
    sd_2 = as.array(0.1),
    # SDs for the dpars since they have (1|ID)
    sd_3 = as.array(0.1), # for bs
    sd_4 = as.array(0.1), # for ndt
    sd_5 = as.array(0.1)  # for bias
  )
}

# ===== SETTING UP LOOP FOR PARAMETER RECOVERY =====

# defining number of fits for parameter recovery
n_fit <- 1 # change for testing

# setting up a grid of entropy difference scores
min_H_diff <- 0
max_H_diff <- 9
H_diff_grid <- seq(min_H_diff, max_H_diff, 1)

# check that we have enough values in the grid
length(H_diff_grid) == n_fit

# defining storage structure
true_drift_intercept <- array(NA, c(n_fit))
true_drift_abs_H_diff <- array(NA, c(n_fit))
true_drift_H_sum <- array(NA, c(n_fit))
true_drift_interaction <- array(NA, c(n_fit))

infer_drift_intercept <- array(NA, c(n_fit))
infer_drift_abs_H_diff <- array(NA, c(n_fit))
infer_drift_H_sum <- array(NA, c(n_fit))
infer_drift_interaction <- array(NA, c(n_fit))

true_bs_intercept <- array(NA, c(n_fit))
true_ndt_intercept <- array(NA, c(n_fit))
true_bias_intercept <- array(NA, c(n_fit))

infer_bs_intercept <- array(NA, c(n_fit))
infer_ndt_intercept <- array(NA, c(n_fit))
infer_bias_intercept <- array(NA, c(n_fit))

# ====== STARTING THE LOOP ; not actually looping anymore ; this is handled in terminal =======

# starting loop
cat("Starting parameter recovery batch with seed offset:", seed_entropy, "\n")
start_time <- Sys.time()

# defining true parameters from fitted model (m4)
# varying abs_H_diff systematically using the grid

# setting fixed effect parameters
fixed_drift_intercept <- 1.77
fixed_drift_abs_H <- seed_entropy
fixed_drift_H_sum <- 0.04
fixed_drift_interaction <- -3.69

fixed_bs <- 3.58
fixed_ndt <- -1.03
fixed_bias <-  0.45

# Extract random effects
id_ranef <- ranef(m4_fit)$ID
stimulus_ranef <- ranef(m4_fit)$stimulus

# Create simulation data
sim_data <- model_data_downsampled %>%
  select(ID, stimulus, abs_H_difference, H_sum) %>%
  rowwise() %>%
  mutate(
    # Drift rate with random intercepts
    drift_intercept_total = fixed_drift_intercept + 
      id_ranef[ID, "Estimate", "Intercept"] + 
      stimulus_ranef[stimulus, "Estimate", "Intercept"],
    
    delta = drift_intercept_total + 
      fixed_drift_abs_H * abs_H_difference +
      fixed_drift_H_sum * H_sum +
      fixed_drift_interaction * (abs_H_difference * H_sum),
    
    # Other parameters with ID random effects
    alpha = fixed_bs + id_ranef[ID, "Estimate", "bs_Intercept"],
    tau = exp(fixed_ndt + id_ranef[ID, "Estimate", "ndt_Intercept"]),
    beta = fixed_bias + id_ranef[ID, "Estimate", "bias_Intercept"]
  ) %>%
  ungroup()

# Simulate responses trial-by-trial
set.seed(123)

# Initialize empty vectors
RT_vec <- numeric(nrow(sim_data))
response_vec <- numeric(nrow(sim_data))

# Loop through each trial
for (i in 1:nrow(sim_data)) {
  sim_result <- rwiener(n = 1,
                        alpha = sim_data$alpha[i],
                        tau = sim_data$tau[i],
                        beta = sim_data$beta[i],
                        delta = sim_data$delta[i])
  
  RT_vec[i] <- abs(sim_result$q)
  response_vec[i] <- as.numeric(sim_result$resp == "upper")
}

# Add to data
sim_data$RT <- RT_vec
sim_data$response <- response_vec

### ==== COMMENTED OUT FOR TESTING =====

# ===== fitting m4 on simulated data =====
fit_sim <- brm(m4_formula,
               data = sim_data,
               family = wiener(link_bs = "identity",
                               link_ndt = "log",
                               link_bias = "identity"),
               prior = m4_prior,
               init = m4_initfun,
               iter = 1000, # 1000
               warmup = 500, # 500
               chains = 4,
               cores = 4,
               refresh = 1,
               # file = "fit_wiener_model_4_sim_3" # we are not saving everything here
               control = list(max_treedepth = 15, adapt_delta = 0.95))

# extracting fixed effect coefficients
sim_coefficients <- fixef(fit_sim)
sim_coefficients_df <- data.frame(t(sim_coefficients))

# adding inferred parameter values to storage structure
infer_drift_intercept <- sim_coefficients_df$Intercept
infer_drift_abs_H_diff <- sim_coefficients_df$abs_H_difference
infer_drift_H_sum <- sim_coefficients_df$H_sum
infer_drift_interaction <- sim_coefficients_df$abs_H_difference.H_sum

infer_bs_intercept <- sim_coefficients_df$bs_Intercept
infer_ndt_intercept <- sim_coefficients_df$ndt_Intercept
infer_bias_intercept <- sim_coefficients_df$bias_Intercept

cat(sprintf("Batch %d, Iteration %d complete\n", seed_entropy, i))


end_time <- Sys.time()
cat("Time elapsed:", format(end_time - start_time), "\n")

# ===== SAVING RESULTS ======
output_file <- sprintf('results/param_recovery_batch_%03d.RData', seed_entropy)
dir.create('results', showWarnings = FALSE)

true_drift_intercept <- fixed_drift_intercept
true_drift_abs_H_diff <- fixed_drift_abs_H
true_drift_H_sum <- fixed_drift_H_sum
true_drift_interaction <- fixed_drift_interaction

true_bs_intercept <- fixed_bs
true_ndt_intercept <- fixed_ndt
true_bias_intercept <- fixed_bias

save(true_drift_intercept, 
     true_drift_abs_H_diff, 
     true_drift_H_sum, 
     true_drift_interaction,
     infer_drift_intercept, 
     infer_drift_abs_H_diff, 
     infer_drift_H_sum, 
     infer_drift_interaction,
     true_bs_intercept, 
     true_ndt_intercept, 
     true_bias_intercept, 
     infer_bs_intercept,
     infer_ndt_intercept, 
     infer_bias_intercept,
     file = output_file)

cat("Results saved to:", output_file, "\n")