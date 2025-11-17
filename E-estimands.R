###############################################################################
# ADEMP — E: ESTIMANDS
# Quantify the true causal effects embedded in the DGM and validate them
# analytically and via Monte Carlo.
#
# This script:
#   - Sets working directory 
#   - Sources the data-generating mechanism (simulate_panel)
#   - Computes true ACEs (X→Y and Y→X) per scenario analytically
#   - Validates them via Monte Carlo under the specified DGM
#   - Saves all outputs to /estimands (no printing)
###############################################################################

# =============================================================================
# 0. PROJECT SETUP — WORKING DIRECTORY, LIBRARIES, FOLDERS
# =============================================================================

# Set working directory 
setwd("/data/cephfs-1/home/users/tato10_c/work/causal_inference_panel_data")

# Load libraries
library(tidyverse)  # dplyr, tibble, tidyr, ggplot2, etc.
library(data.table)
library(purrr)
library(glue)
library(readr)

# Top-level folders (flat structure)
save_dir   <- "data"        # simulated datasets (RDS)
meta_dir   <- "meta"        # general metadata (if needed elsewhere)
estim_dir  <- "estimands"   # estimands CSV outputs

dir.create(save_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(meta_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(estim_dir, showWarnings = FALSE, recursive = TRUE)

# Source the data-generating mechanism (must define simulate_panel)
# This file should live in the project root as:
#   D-data_generating_mechanism.R
source("D-data_generating_mechanism.R")


# =============================================================================
# 1. TRUE STRUCTURAL COEFFICIENTS (GROUND TRUTH)
# =============================================================================
# From the DGM:
#   X_{t+1} = φ_X * X_t + β_YX * Y_t + γ_cX * c_t + ε_X
#   Y_{t+1} = φ_Y * Y_t + β_XY * X_t + γ_cY * c_t + ε_Y
#
# We define:
#   ACE_{X→Y}(t) = β_XY
#   ACE_{Y→X}(t) = β_YX
#
# Scenario logic encoded as:
#   Scenario 1: No causal effect
#       β_XY = 0, β_YX = 0
#   Scenario 2: X → Y only
#       β_XY = 0.8, β_YX = 0
#   Scenario 3: X ↔ Y bidirectional
#       β_XY = 0.8, β_YX = 0.5
###############################################################################

true_estimands <- function(scenario, beta_xy = NULL, beta_yx = NULL) {
  # Use defaults from DGM if user does not override
  bxy <- if (is.null(beta_xy)) ifelse(scenario >= 2, 0.8, 0) else beta_xy
  byx <- if (is.null(beta_yx)) ifelse(scenario == 3, 0.5, 0) else beta_yx
  
  list(
    ACE_XY = bxy,   # Average Causal Effect of X on Y
    ACE_YX = byx    # Average Causal Effect of Y on X
  )
}


# =============================================================================
# 2. ANALYTIC ESTIMANDS (CLOSED-FORM "TRUTH")
# =============================================================================
# These are the exact ACE values implied by the DGM, without simulation.
#
# For a given scenario:
#   ACE_XY_true(t) = β_XY       (time-invariant under your DGM)
#   ACE_YX_true(t) = β_YX
#
# We also allow for an intervention size Δ (delta) if you want to think in
# terms of moving X_t or Y_t by Δ units.
###############################################################################

validate_estimands_analytic <- function(scenario, T_obs = 4, delta = 1) {
  par <- true_estimands(scenario)
  
  tibble(
    scenario         = scenario,
    time_intervened  = 0:(T_obs - 1),   # t where X_t or Y_t is perturbed
    delta            = delta,
    ACE_XY_true      = par$ACE_XY,      # true structural parameter
    ACE_YX_true      = par$ACE_YX,
    ACE_XY_estimated = par$ACE_XY * delta,  # trivial identity here
    ACE_YX_estimated = par$ACE_YX * delta,
    method           = "analytic"
  )
}


# =============================================================================
# 3. MONTE CARLO ESTIMANDS (NUMERICAL VALIDATION)
# =============================================================================
# Purpose:
#   Numerically check that the ACEs embedded in the DGM behave as expected.
#
# Strategy:
#   1) Simulate a large panel dataset under a given scenario.
#   2) Take values at time t and t+1.
#   3) Compute potential outcomes under a Δ intervention on X_t (for ACE_XY)
#      and on Y_t (for ACE_YX).
#   4) Average differences to get empirical ACE estimates.
###############################################################################

validate_estimands_mc <- function(
    scenario = 3,
    N       = 100000,   # large N for Monte Carlo precision
    T_obs   = 4,
    t       = 2,        # intervene at time t
    delta   = 1,
    phi_x   = 0.5,
    phi_y   = 0.5,
    gamma_cx = 0.3,
    gamma_cy = 0.4,
    sigma_x  = 1,
    sigma_y  = 1,
    sigma_c  = 0.5
) {
  # Retrieve true coefficients
  par <- true_estimands(scenario)
  bxy <- par$ACE_XY
  byx <- par$ACE_YX
  
  # Step 1: simulate large panel
  base <- simulate_panel(
    scenario = scenario,
    N        = N,
    T_obs    = T_obs,
    phi_x    = phi_x,
    phi_y    = phi_y,
    gamma_cx = gamma_cx,
    gamma_cy = gamma_cy,
    sigma_x  = sigma_x,
    sigma_y  = sigma_y
  )
  
  # Step 2: keep time t and t+1 and go wide
  d <- base %>%
    filter(time %in% c(t, t + 1)) %>%
    select(id, time, X, Y, c) %>%
    pivot_wider(
      names_from  = time,
      values_from = c(X, Y, c),
      names_sep   = "_"
    )
  
  Xt <- glue("X_{t}")
  Yt <- glue("Y_{t}")
  Ct <- glue("c_{t}")   # confounder at time t
  
  # Step 3: Potential outcomes for Y_{t+1} under X_t and X_t + Δ
  Y_next0 <- phi_y * d[[Yt]] + bxy * d[[Xt]]           + gamma_cy * d[[Ct]]
  Y_next1 <- phi_y * d[[Yt]] + bxy * (d[[Xt]] + delta) + gamma_cy * d[[Ct]]
  
  # Potential outcomes for X_{t+1} under Y_t and Y_t + Δ
  X_next0 <- phi_x * d[[Xt]] + byx * d[[Yt]]           + gamma_cx * d[[Ct]]
  X_next1 <- phi_x * d[[Xt]] + byx * (d[[Yt]] + delta) + gamma_cx * d[[Ct]]
  
  tibble(
    scenario         = scenario,
    time_intervened  = t,
    delta            = delta,
    ACE_XY_true      = bxy,
    ACE_YX_true      = byx,
    ACE_XY_estimated = mean(Y_next1 - Y_next0),
    ACE_YX_estimated = mean(X_next1 - X_next0),
    method           = "montecarlo"
  )
}


# =============================================================================
# 4. RUN ANALYTIC ESTIMANDS FOR ALL SCENARIOS
# =============================================================================

results_estimands_full <- map_dfr(
  1:3,
  ~validate_estimands_analytic(.x, T_obs = 4, delta = 1)
)

# Persist analytic table to disk
write_csv(
  results_estimands_full,
  file.path(estim_dir, "analytic_estimands.csv")
)


# =============================================================================
# 5. RUN MONTE CARLO VALIDATION (HEAVIER, OPTIONAL)
# =============================================================================
# For each scenario:
#   - Evaluate ACEs at time t = 1 and t = 3
#   - N = 20,000 per check (tune N up/down for precision vs. compute cost)
###############################################################################

mc_checks <- map_dfr(
  1:3,
  ~bind_rows(
    validate_estimands_mc(scenario = .x, t = 1, N = 20000),
    validate_estimands_mc(scenario = .x, t = 3, N = 20000)
  )
)

write_csv(
  mc_checks,
  file.path(estim_dir, "montecarlo_estimands.csv")
)


# =============================================================================
# 6. FINAL COMBINED ESTIMANDS TABLE
# =============================================================================
# This is the human-readable “truth table” we will compare estimators against
# later in the ADEMP pipeline.
###############################################################################

final_estimands <- results_estimands_full %>%
  mutate(
    scenario_label = recode(
      as.character(scenario),
      "1" = "S1: No causal effect",
      "2" = "S2: X → Y",
      "3" = "S3: X ↔ Y"
    )
  ) %>%
  select(
    scenario_label,
    time_intervened,
    delta,
    ACE_XY_true,
    ACE_XY_estimated,
    ACE_YX_true,
    ACE_YX_estimated,
    method
  )

# Save combined summary
write_csv(
  final_estimands,
  file.path(estim_dir, "final_estimands_summary.csv")
)

###############################################################################
# END OF ADEMP — E: ESTIMANDS
###############################################################################
