#### FIT STAN MODEL - for Snake River, Wild DPS ####
# CHAIN 4: Seed 104

library(cmdstanr)
library(posterior)
library(tidyverse)
library(lubridate)

# Step 0: load the data
load(file = "SRW_T_model_data.rda")



# Fit stan model using cmdstan
# Step 1: load the model
mod <- cmdstan_model("SRW_T_model.stan", compile = FALSE)


# Step 2: Compile the model, set up to run in parallel
print("Compiling model")
mod$compile(cpp_options = list(stan_threads = TRUE))

print("model compiled at")
print(Sys.time())





# Step 3: Run MCMC (HMC) - One chain, using all of the cores on a node
# each R script runs one chain, using different seed/init values; these are joined later
fit <- mod$sample(
  data = data,
  seed = 104,
  # seed = 456,
  # chains = 3,
  chains = 1,
  parallel_chains = 1,
  # parallel_chains = 3,
  refresh = 100, # print update every 10 iter
  iter_warmup = 2000,
  iter_sampling = 2000,
  # thin = 10,
  max_treedepth = 12,
  # adapt_delta = 0.95,
  init = 1, # I believe that this needs to not be zero to ensure that each chain starts in a different place
  # from mc-stan.org: A real number x>0. This initializes all parameters randomly between [-x,x] on the unconstrained parameter space.
  threads_per_chain = 28
)

fit$save_object(file = "SRW_T_chain4_fit.rds")
