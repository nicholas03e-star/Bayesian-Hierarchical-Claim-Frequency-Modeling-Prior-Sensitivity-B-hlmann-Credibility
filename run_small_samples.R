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

R <- n_distinct(df$Region_idx)

# proportional n by region
prop_targets <- function(df, n_total) {
  df %>%
    group_by(Region_idx) %>%
    summarise(n_region = n(), .groups = "drop") %>%
    mutate(n_target = pmax(1, round(n_region / sum(n_region) * n_total)))
}

# heir vs flat
hier_model <- cmdstan_model(
  "C:/Users/Nicholas/OneDrive/Documents/hierarchical_poisson.stan",
  cpp_options = list(stan_threads = TRUE)
)

flat_model <- cmdstan_model(
  "C:/Users/Nicholas/OneDrive/Documents/flat_poisson.stan",
  cpp_options = list(stan_threads = TRUE)
)

# sample and fit
run_sample <- function(n_total, label, sigma_scale = NULL) {
  
  is_flat   <- is.null(sigma_scale)
  prior_tag <- if (is_flat) "flat" else if (sigma_scale == 0.5) "weak" else "inform"
  
  cat("Starting:", label, "| model:", prior_tag, "\n")
  
  set.seed(3)
  targets <- prop_targets(df, n_total)
  
  df_s <- df %>%
    left_join(targets, by = "Region_idx") %>%
    group_by(Region_idx) %>%
    group_modify(~ slice_sample(.x, n = .x$n_target[1])) %>%
    ungroup()
  
  saveRDS(df_s, paste0("C:/Users/Nicholas/OneDrive/Documents/df_prop_", label, ".rds"))
  
  stan_data_s <- list(
    N        = nrow(df_s),
    R        = R,
    y        = as.integer(df_s$ClaimNb),
    exposure = df_s$Exposure,
    X        = as.matrix(dplyr::select(df_s, DrivAge, BonusMalus, VehPower,
                                       VehAge, log_density, VehGas_bin)),
    K        = 6,
    region   = as.integer(df_s$Region_idx)
  )

  if (!is_flat) stan_data_s <- c(stan_data_s, list(sigma_scale = sigma_scale))
  
  model      <- if (is_flat) flat_model else hier_model
  output_dir <- paste0("D:/stan_output_prop_", label, "_", prior_tag)
  dir.create(output_dir, showWarnings = FALSE)
  
  fit <- model$sample(
    data              = stan_data_s,
    chains            = 4,
    parallel_chains   = 4,
    threads_per_chain = 2,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    seed              = 3,
    refresh           = 200,
    output_dir        = output_dir
  )
  
  variables <- if (is_flat) c("alpha", "beta") else c("alpha", "beta", "sigma_u", "u")
  
  files <- list.files(output_dir, pattern = "\\.csv$", full.names = TRUE)
  files <- files[!grepl("config|metric", files)]
  
  fit_out <- read_cmdstan_csv(
    files               = files,
    variables           = variables,
    sampler_diagnostics = c("divergent__", "treedepth__")
  )
  
  saveRDS(fit_out, paste0("C:/Users/Nicholas/OneDrive/Documents/fit_prop_",
                          label, "_", prior_tag, ".rds"))
  
  cat("Done:", label, "| model:", prior_tag, "\n")
}

# combos
sizes <- list(c(500, "500"), c(2000, "2k"), c(10000, "10k"))

for (s in sizes) {
  run_sample(as.integer(s[1]), s[2], sigma_scale = NULL)  # flat
  run_sample(as.integer(s[1]), s[2], sigma_scale = 0.5)   # weak
  run_sample(as.integer(s[1]), s[2], sigma_scale = 0.1)   # informative
}

cat("All done.\n")