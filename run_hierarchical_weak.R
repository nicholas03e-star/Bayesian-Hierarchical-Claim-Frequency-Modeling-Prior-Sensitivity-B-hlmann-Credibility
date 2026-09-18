library(tidyverse)
library(cmdstanr)

freMTPL2freq <- read_csv("C:/Users/Nicholas/OneDrive/Documents/freMTPL2freq.csv")

df <- freMTPL2freq %>%
  filter(Exposure > 0, Exposure <= 1) %>%
  mutate(
    log_density = log(Density),
    VehGas_bin  = ifelse(VehGas == "Diesel", 1, 0),
    Region_idx  = as.integer(factor(Region))
  ) %>%
  mutate(across(
    c(DrivAge, BonusMalus, VehPower, VehAge, log_density),
    scale
  ))

N <- nrow(df)
R <- n_distinct(df$Region_idx)

stan_data_weak <- list(
  N          = N,
  R          = R,
  y          = as.integer(df$ClaimNb),
  exposure   = df$Exposure,
  X          = as.matrix(df %>% select(DrivAge, BonusMalus, VehPower, VehAge, log_density, VehGas_bin)),
  K          = 6,
  region     = as.integer(df$Region_idx),
  sigma_scale = 0.5
)

hier_model <- cmdstan_model(
  "C:/Users/Nicholas/OneDrive/Documents/hierarchical_poisson.stan",
  cpp_options = list(stan_threads = TRUE)
)

fit_hier_weak <- hier_model$sample(
  data              = stan_data_weak,
  chains            = 4,
  parallel_chains   = 4,
  threads_per_chain = 2,
  iter_warmup       = 1000,
  iter_sampling     = 1000,
  seed              = 3,
  refresh           = 200,
  output_dir = "D:/stan_output"
)
fit_hier_weak$save_object("C:/Users/Nicholas/OneDrive/Documents/fit_hier_weak.rds")
cat("Done.\n")