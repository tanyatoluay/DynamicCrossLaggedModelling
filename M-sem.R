###############################################################################
# SEM MONTE CARLO: CROSS-LAGGED PANEL WITH TIME-VARYING CONFOUNDER
###############################################################################

# 0. Setup --------------------------------------------------------------------

setwd("/data/cephfs-1/home/users/tato10_c/work/causal_inference_panel_data")

library(lavaan)
library(tidyverse)
library(data.table)

# Directories -----------------------------------------------------------------
data_dir   <- file.path(getwd(), "data")           
method_dir <- file.path(getwd(), "methods", "sem")
dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)

# Monte Carlo setup -----------------------------------------------------------
nsim <- 2000

###############################################################################
# 1. Helper: true causal estimands per scenario
###############################################################################
true_estimands <- function(scenario) {
  if (scenario == 1) {
    tibble(ACE_XY = 0,   ACE_YX = 0)
  } else if (scenario == 2) {
    tibble(ACE_XY = 0.8, ACE_YX = 0)
  } else if (scenario == 3) {
    tibble(ACE_XY = 0.8, ACE_YX = 0.5)
  } else {
    stop("Scenario must be 1, 2, or 3")
  }
}

###############################################################################
# 2. Helper: load dataset
###############################################################################
get_dataset <- function(scenario, rep_id) {
  f <- file.path(data_dir, sprintf("dgm%i_rep%i.rds", scenario, rep_id))
  if (!file.exists(f)) {
    stop(glue::glue("File not found: {f}"))
  }
  readRDS(f)
}

###############################################################################
# 3. SEM model syntax
###############################################################################
sem_model <- "
  X_1 ~ gamma_cx*c_1
  Y_1 ~ gamma_cy*c_1

  c_2 ~ 1*c_1
  c_3 ~ 1*c_2
  c_4 ~ 1*c_3

  X_2 ~ phi_x*X_1 + beta_yx*Y_1 + gamma_cx*c_1
  Y_2 ~ phi_y*Y_1 + beta_xy*X_1 + gamma_cy*c_1

  X_3 ~ phi_x*X_2 + beta_yx*Y_2 + gamma_cx*c_2
  Y_3 ~ phi_y*Y_2 + beta_xy*X_2 + gamma_cy*c_2

  X_4 ~ phi_x*X_3 + beta_yx*Y_3 + gamma_cx*c_3
  Y_4 ~ phi_y*Y_3 + beta_xy*X_3 + gamma_cy*c_3

  X_1 ~~ Y_1
  X_2 ~~ Y_2
  X_3 ~~ Y_3
  X_4 ~~ Y_4
"

###############################################################################
# 4. Fit SEM for one (scenario, rep)
###############################################################################
fit_sem_scenario <- function(scenario, rep_id) {
  
  dat <- get_dataset(scenario, rep_id)
  setDT(dat)
  
  df_wide <- dat %>%
    select(id, time, X, Y, c) %>%
    pivot_wider(
      names_from  = time,
      values_from = c(X, Y, c),
      names_sep   = "_"
    )
  
  fit <- sem(
    model   = sem_model,
    data    = df_wide,
    fixed.x = FALSE,
    missing = "ML"
  )
  
  # Save full model object
  saveRDS(
    fit,
    file.path(method_dir, sprintf("sem_fit_s%i_rep%i.rds", scenario, rep_id))
  )
  
  # Extract raw (unstandardized) estimates
  est_table <- parameterEstimates(fit, standardized = FALSE, ci = TRUE)
  fit_meas  <- fitMeasures(fit)
  
  # beta_xy / beta_yx appear 3× each → take the first
  est_xy_yx <- est_table %>%
    filter(label %in% c("beta_xy", "beta_yx")) %>%
    group_by(label) %>%
    slice(1) %>%
    ungroup() %>%
    select(label, est, se, ci.lower, ci.upper)
  
  truth <- true_estimands(scenario)
  
  tibble(
    scenario = scenario,
    rep_id   = rep_id,
    
    beta_xy_hat  = est_xy_yx$est[est_xy_yx$label == "beta_xy"],
    beta_xy_se   = est_xy_yx$se[est_xy_yx$label == "beta_xy"],
    beta_xy_low  = est_xy_yx$ci.lower[est_xy_yx$label == "beta_xy"],
    beta_xy_high = est_xy_yx$ci.upper[est_xy_yx$label == "beta_xy"],
    beta_xy_true = truth$ACE_XY,
    
    beta_yx_hat  = est_xy_yx$est[est_xy_yx$label == "beta_yx"],
    beta_yx_se   = est_xy_yx$se[est_xy_yx$label == "beta_yx"],
    beta_yx_low  = est_xy_yx$ci.lower[est_xy_yx$label == "beta_yx"],
    beta_yx_high = est_xy_yx$ci.upper[est_xy_yx$label == "beta_yx"],
    beta_yx_true = truth$ACE_YX,
    
    par_table    = list(est_table),
    fit_measures = list(as.list(fit_meas))
  )
}

###############################################################################
# 5. Run SEM across scenarios + replicates
###############################################################################
set.seed(20251026)

sem_results <- map_dfr(1:3, function(s) {
  message("Running SEM — Scenario ", s)
  map_dfr(1:nsim, function(r) {
    message("  Replicate ", r, "/", nsim)
    fit_sem_scenario(s, r)
  })
})

###############################################################################
# 6. MC performance
###############################################################################
sem_performance <- sem_results %>%
  summarise(
    bias_XY = mean(beta_xy_hat - beta_xy_true),
    bias_YX = mean(beta_yx_hat - beta_yx_true),
    var_XY  = var(beta_xy_hat),
    var_YX  = var(beta_yx_hat),
    mse_XY  = mean((beta_xy_hat - beta_xy_true)^2),
    mse_YX  = mean((beta_yx_hat - beta_yx_true)^2)
  )

###############################################################################
# 7. Human-readable summary table
###############################################################################
sem_output <- sem_results %>%
  mutate(
    scenario_label = recode(
      scenario,
      `1` = "S1: No causal effect",
      `2` = "S2: X → Y",
      `3` = "S3: X ↔ Y"
    )
  ) %>%
  select(
    scenario_label, rep_id,
    beta_xy_hat, beta_xy_se, beta_xy_low, beta_xy_high, beta_xy_true,
    beta_yx_hat, beta_yx_se, beta_yx_low, beta_yx_high, beta_yx_true
  )

###############################################################################
# 8. Persist outputs
###############################################################################
saveRDS(sem_results,     file.path(method_dir, "sem_results_all.rds"))
saveRDS(sem_performance, file.path(method_dir, "sem_performance.rds"))
saveRDS(sem_output,      file.path(method_dir, "sem_output_table.rds"))

readr::write_csv(sem_output,
                 file.path(method_dir, "sem_output_table.csv"))
readr::write_csv(sem_performance,
                 file.path(method_dir, "sem_performance.csv"))
