###############################################################################
# ADEMP: Bayesian DYNAMITE Method (single-scenario core, with checkpoints)
###############################################################################

setwd("/data/cephfs-1/home/users/tato10_c/work/causal_inference_panel_data")

library(dynamite)
library(tidyverse)
library(data.table)
library(purrr)
library(stringr)

## ---------------------------------------------------------------------------
## SCENARIO must be defined by the calling script (1, 2, or 3)
## ---------------------------------------------------------------------------

if (!exists("SCENARIO")) {
  stop("SCENARIO not set. Define SCENARIO in the wrapper script before sourcing dynamite_common.R.")
}

## ---------------------------------------------------------------------------
## Directories
## ---------------------------------------------------------------------------

data_dir    <- file.path(getwd(), "data")
method_dir  <- file.path(getwd(), "methods", "dynamite")
fits_dir    <- file.path(method_dir, "fits")
post_dir    <- file.path(method_dir, "posteriors")
plots_dir   <- file.path(method_dir, "plots")      # kept, but unused now
perf_dir    <- file.path(method_dir, "performance")

dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fits_dir,   showWarnings = FALSE)
dir.create(post_dir,   showWarnings = FALSE)
dir.create(plots_dir,  showWarnings = FALSE)
dir.create(perf_dir,   showWarnings = FALSE)

## ---------------------------------------------------------------------------
## Config (lighter + fewer reps)
## ---------------------------------------------------------------------------

nsim             <- 2000L          
scenarios        <- SCENARIO      # single scenario per run
checkpoint_every <- 100L

# MCMC settings: lighter but still reasonable for CI
mcmc_chains   <- 2L
mcmc_iter     <- 2000L
mcmc_warmup   <- 1000L
mcmc_adapt    <- 0.95
mcmc_treedepth <- 12L

###############################################################################
# 1. True causal estimands
###############################################################################

true_estimands <- function(scenario) {
  tibble(
    ACE_XY = ifelse(scenario >= 2, 0.8, 0),
    ACE_YX = ifelse(scenario == 3, 0.5, 0)
  )
}

###############################################################################
# 2. Dataset loading
###############################################################################

get_dataset <- function(scenario, rep_id) {
  f <- file.path(data_dir, sprintf("dgm%i_rep%i.rds", scenario, rep_id))
  if (!file.exists(f)) stop("Dataset missing: ", f)
  readRDS(f)
}

###############################################################################
# 3. Standardization
###############################################################################

standardize_panel <- function(dat) {
  dat %>%
    mutate(
      X_sd   = sd(X, na.rm = TRUE),
      Y_sd   = sd(Y, na.rm = TRUE),
      c_sd   = sd(c, na.rm = TRUE),
      X_mean = mean(X),
      Y_mean = mean(Y),
      c_mean = mean(c)
    ) %>%
    mutate(
      X_std = (X - X_mean) / X_sd,
      Y_std = (Y - Y_mean) / Y_sd,
      c_std = (c - c_mean) / c_sd
    )
}

###############################################################################
# 4. Fit one replicate (no plotting)
###############################################################################

fit_dynamite_scenario <- function(scenario, rep_id, seed = NULL) {
  
  dat_raw <- get_dataset(scenario, rep_id)
  dat_std <- standardize_panel(dat_raw)
  setDT(dat_std)
  
  if (is.null(seed))
    seed <- 20260000L + scenario * 100000L + rep_id
  
  set.seed(seed)
  
  model_formula <-
    obs(Y_std ~ lag(Y_std) + lag(X_std) + lag(c_std), family = "gaussian") +
    obs(X_std ~ lag(X_std) + lag(Y_std) + lag(c_std), family = "gaussian")
  
  dyn_fit <- dynamite(
    dformula = model_formula,
    data     = dat_std,
    time     = "time",
    group    = "id",
    chains   = mcmc_chains,
    iter     = mcmc_iter,
    warmup   = mcmc_warmup,
    refresh  = 0,
    control  = list(adapt_delta = mcmc_adapt,
                    max_treedepth = mcmc_treedepth),
    verbose  = FALSE
    # you can add: cores = mcmc_chains, threads_per_chain = 1L, ...
  )
  
  ## Save fit
  fit_path <- file.path(fits_dir, sprintf("dyn_fit_s%i_rep%i.rds", scenario, rep_id))
  saveRDS(dyn_fit, fit_path)
  
  ## Posterior summaries for standardized betas
  post_sum <- as.data.frame(dyn_fit, types = "beta", summary = TRUE)
  post_path <- file.path(post_dir, sprintf("dyn_post_std_s%i_rep%i.csv", scenario, rep_id))
  readr::write_csv(post_sum, post_path)
  
  ## Plotting removed to speed things up
  
  ## Extract standardized key params
  beta_xy <- post_sum[post_sum$parameter == "beta_Y_std_X_std_lag1", ]
  beta_yx <- post_sum[post_sum$parameter == "beta_X_std_Y_std_lag1", ]
  
  sd_X <- unique(dat_std$X_sd)
  sd_Y <- unique(dat_std$Y_sd)
  
  k_xy <- sd_Y / sd_X
  k_yx <- sd_X / sd_Y
  
  beta_xy_hat <- beta_xy$mean * k_xy
  beta_xy_sd  <- beta_xy$sd   * k_xy
  beta_xy_q05 <- beta_xy$q5   * k_xy
  beta_xy_q95 <- beta_xy$q95  * k_xy
  
  beta_yx_hat <- beta_yx$mean * k_yx
  beta_yx_sd  <- beta_yx$sd   * k_yx
  beta_yx_q05 <- beta_yx$q5   * k_yx
  beta_yx_q95 <- beta_yx$q95  * k_yx
  
  truth <- true_estimands(scenario)
  
  tibble(
    scenario     = scenario,
    rep_id       = rep_id,
    mcmc_seed    = seed,
    beta_xy_hat  = beta_xy_hat,
    beta_xy_sd   = beta_xy_sd,
    beta_xy_q05  = beta_xy_q05,
    beta_xy_q95  = beta_xy_q95,
    beta_yx_hat  = beta_yx_hat,
    beta_yx_sd   = beta_yx_sd,
    beta_yx_q05  = beta_yx_q05,
    beta_yx_q95  = beta_yx_q95,
    beta_xy_true = truth$ACE_XY,
    beta_yx_true = truth$ACE_YX
  )
}

###############################################################################
# 5. Mapping global index → (scenario, rep)
###############################################################################

build_global_map <- function(scenarios, nsim) {
  expand.grid(
    scenario = scenarios,
    rep_id   = seq_len(nsim)
  ) %>%
    arrange(scenario, rep_id) %>%
    mutate(global_id = row_number())
}

global_map <- build_global_map(scenarios, nsim)

###############################################################################
# 6. Auto-Resume Logic (scenario-specific checkpoints)
###############################################################################

find_resume_point <- function(scenario_id) {
  pat <- sprintf("^dynamite_results_s%i_checkpoint_\\d+\\.rds$", scenario_id)
  ckpts <- list.files(method_dir, pattern = pat)
  
  if (length(ckpts) == 0L) return(0L)
  
  ids <- sub(sprintf("dynamite_results_s%i_checkpoint_(\\d+)\\.rds", scenario_id),
             "\\1", ckpts)
  max(as.integer(ids), na.rm = TRUE)
}

resume_global_id <- find_resume_point(SCENARIO)

###############################################################################
# 7. Load previous partial results (if any)
###############################################################################

dyn_results_list <- list()

if (resume_global_id > 0) {
  
  ckpt_file <- file.path(
    method_dir,
    sprintf("dynamite_results_s%i_checkpoint_%05d.rds", SCENARIO, resume_global_id)
  )
  
  message("RESUMING scenario ", SCENARIO, " from checkpoint: ", ckpt_file)
  
  prev <- readRDS(ckpt_file)
  dyn_results_list <- as.list(split(prev, seq_len(nrow(prev)))) 
}

###############################################################################
# 8. Continue computation
###############################################################################

message("Starting/continuing Monte Carlo for scenario ", SCENARIO, "…")

start_at <- resume_global_id + 1L
end_at   <- nrow(global_map)

for (gid in start_at:end_at) {
  
  row <- global_map[gid, ]
  s   <- row$scenario
  r   <- row$rep_id
  
  message("Global ", gid, "/", end_at, " — Scenario ", s, ", Rep ", r)
  
  ## **Skip if fit exists**  
  fit_out <- file.path(fits_dir, sprintf("dyn_fit_s%i_rep%i.rds", s, r))
  if (file.exists(fit_out)) {
    message("  Found existing fit — SKIPPING")
    next
  }
  
  res <- fit_dynamite_scenario(s, r)
  dyn_results_list[[length(dyn_results_list)+1]] <- res
  
  ## Save checkpoint (scenario-specific)
  if (gid %% checkpoint_every == 0L) {
    current <- bind_rows(dyn_results_list)
    ckpt_rds <- file.path(
      method_dir,
      sprintf("dynamite_results_s%i_checkpoint_%05d.rds", SCENARIO, gid)
    )
    ckpt_csv <- file.path(
      method_dir,
      sprintf("dynamite_results_s%i_checkpoint_%05d.csv", SCENARIO, gid)
    )
    saveRDS(current, ckpt_rds)
    readr::write_csv(current, ckpt_csv)
    message("  -> CHECKPOINT SAVED @ global ", gid,
            " (scenario ", SCENARIO, ")")
  }
}

###############################################################################
# 9. Final results assembly
###############################################################################

dyn_results <- bind_rows(dyn_results_list)

saveRDS(
  dyn_results,
  file.path(method_dir, sprintf("dynamite_results_s%i_all.rds", SCENARIO))
)
readr::write_csv(
  dyn_results,
  file.path(method_dir, sprintf("dynamite_results_s%i_all.csv", SCENARIO))
)

###############################################################################
# 10. Monte Carlo performance metrics
###############################################################################

mc_se <- function(estimates, theta, ci_low, ci_high) {
  nsim <- length(estimates)
  bias <- mean(estimates) - theta
  bias_mcse <- sqrt(var(estimates) / nsim)
  empSE <- sd(estimates)
  empSE_mcse <- empSE * sqrt(1 / (2 * (nsim - 1)))
  MSE <- mean((estimates - theta)^2)
  MSE_mcse <- sqrt(sum(((estimates - theta)^2 - MSE)^2) /
                     (nsim * (nsim - 1)))
  cover <- mean(ci_low <= theta & ci_high >= theta)
  cover_mcse <- sqrt(cover * (1 - cover) / nsim)
  data.frame(
    bias, bias_mcse,
    empSE, empSE_mcse,
    MSE, MSE_mcse,
    cover, cover_mcse
  )
}

setDT(dyn_results)

dyn_perf <- dyn_results[, {
  theta_xy <- unique(beta_xy_true)
  theta_yx <- unique(beta_yx_true)
  xy <- mc_se(beta_xy_hat, theta_xy, beta_xy_q05, beta_xy_q95)
  yx <- mc_se(beta_yx_hat, theta_yx, beta_yx_q05, beta_yx_q95)
  cbind(effect = c("XY", "YX"), rbind(xy, yx))
}, by = scenario]

dyn_perf <- dyn_perf %>%
  mutate(
    scenario_label = recode(
      as.character(scenario),
      "1" = "S1: No causal effect",
      "2" = "S2: X → Y",
      "3" = "S3: X ↔ Y"
    )
  ) %>%
  select(scenario, scenario_label, effect,
         bias, bias_mcse, empSE, empSE_mcse,
         MSE, MSE_mcse, cover, cover_mcse)

saveRDS(
  dyn_perf,
  file.path(perf_dir, sprintf("dynamite_performance_s%i.rds", SCENARIO))
)
readr::write_csv(
  dyn_perf,
  file.path(perf_dir, sprintf("dynamite_performance_s%i.csv", SCENARIO))
)

message("FINISHED scenario ", SCENARIO, ". Outputs in: ", method_dir)
