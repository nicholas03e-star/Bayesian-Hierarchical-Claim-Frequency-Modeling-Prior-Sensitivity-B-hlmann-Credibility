library(tidyverse)
library(cmdstanr)

# data
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

stan_data <- list(
  N        = N,
  R        = R,
  y        = as.integer(df$ClaimNb),
  exposure = df$Exposure,
  X        = as.matrix(df %>% select(DrivAge, BonusMalus, VehPower, VehAge, log_density, VehGas_bin)),
  K        = 6,
  region   = as.integer(df$Region_idx)
)

stan_data_inform <- c(stan_data, list(sigma_scale = 0.1))

# compile
flat_model <- cmdstan_model(
  "C:/Users/Nicholas/OneDrive/Documents/flat_poisson.stan",
  cpp_options = list(stan_threads = TRUE)
)

hier_model <- cmdstan_model(
  "C:/Users/Nicholas/OneDrive/Documents/hierarchical_poisson.stan",
  cpp_options = list(stan_threads = TRUE)
)

#  flat
cat("Starting flat model...\n")
fit_flat <- flat_model$sample(
  data              = stan_data,
  chains            = 4,
  parallel_chains   = 4,
  threads_per_chain = 2,
  iter_warmup       = 1000,
  iter_sampling     = 1000,
  seed              = 3,
  refresh           = 200,
  output_dir        = "D:/stan_output_flat"
)

flat_files <- list.files("D:/stan_output_flat",
                         pattern = "\\.csv$",
                         full.names = TRUE)
flat_files <- flat_files[!grepl("config|metric", flat_files)]

fit_flat <- read_cmdstan_csv(
  files               = flat_files,
  variables           = c("alpha", "beta"),
  sampler_diagnostics = c("divergent__", "treedepth__")
)
saveRDS(fit_flat, "C:/Users/Nicholas/OneDrive/Documents/fit_flat.rds")
cat("Flat done.\n")

# informative
cat("Starting informative hierarchical...\n")
fit_hier_inform <- hier_model$sample(
  data              = stan_data_inform,
  chains            = 4,
  parallel_chains   = 4,
  threads_per_chain = 2,
  iter_warmup       = 1000,
  iter_sampling     = 1000,
  seed              = 3,
  refresh           = 200,
  output_dir        = "D:/stan_output_inform"
)

inform_files <- list.files("D:/stan_output_inform",
                           pattern = "\\.csv$",
                           full.names = TRUE)
inform_files <- inform_files[!grepl("config|metric", inform_files)]

fit_hier_inform <- read_cmdstan_csv(
  files               = inform_files,
  variables           = c("alpha", "beta", "sigma_u", "u"),
  sampler_diagnostics = c("divergent__", "treedepth__")
)
saveRDS(fit_hier_inform, "C:/Users/Nicholas/OneDrive/Documents/fit_hier_inform.rds")
cat("All done.\n")