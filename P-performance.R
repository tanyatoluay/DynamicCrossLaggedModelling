###############################################################################
# Master performance + plotting script for SEM vs DYNAMITE (1000 runs)
# Project root: /data/cephfs-1/home/users/tato10_c/work/causal_inference_panel_data
#
# Functionality:
#   - Load SEM 1000-run results (single file)
#   - Load and merge DYNAMITE 1000-run results (multiple checkpoint files)
#   - Build combined long-format dataset
#   - Compute full set of performance metrics (Morris-style)
#   - Generate tables (SEM / DYNAMITE / combined)
#   - Generate plots:
#       * Density of estimates
#       * Boxplots of estimates
#       * Jittered estimates by replicate
#       * Bias / RMSE / Coverage barplots
#       * ZIP plots
#       * Lollipop performance summaries
#       * SE(model) density
#       * SE(model) vs estimate
#       * SEM vs DYNAMITE limits-of-agreement
###############################################################################

## 0. Setup -------------------------------------------------------------------

setwd("/data/cephfs-1/home/users/tato10_c/work/causal_inference_panel_data")

library(tidyverse)
library(data.table)
library(lubridate)
library(ggplot2)

root_dir      <- getwd()
sem_dir       <- file.path(root_dir, "methods", "sem")
dyn_dir       <- file.path(root_dir, "methods", "dynamite")

perf_root     <- file.path(root_dir, "performance")
tables_dir    <- file.path(perf_root, "tables")
plots_dir     <- file.path(perf_root, "plots")
combined_dir  <- file.path(perf_root, "combined")

dir.create(perf_root,   showWarnings = FALSE)
dir.create(tables_dir,  showWarnings = FALSE)
dir.create(plots_dir,   showWarnings = FALSE)
dir.create(combined_dir,showWarnings = FALSE)

###############################################################################
# 1. Utility: labels, Morris colour scheme, MC performance function
###############################################################################

# Scenario labels for plotting
scenario_labels <- c(
  "1" = "S1: No causal effect",
  "2" = "S2: X → Y",
  "3" = "S3: X ↔ Y"
)

# Morris-style colour palette
morris_cols <- c(
  "SEM"      = "#1B63A0",   # deep blue
  "DYNAMITE" = "#7F7F7F"    # neutral grey
)

morris_fill <- c(
  "SEM"      = "#1B63A0",
  "DYNAMITE" = "#B1B1B1"
)

# Monte Carlo performance summary (extended with ModSE etc.)
mc_se <- function(estimates, theta, ci_low, ci_high, se_mod) {
  nsim <- length(estimates)
  
  # Basic quantities
  bias <- mean(estimates) - theta
  bias_mcse <- sqrt(var(estimates) / nsim)
  
  empSE <- sd(estimates)
  empSE_mcse <- empSE * sqrt(1 / (2 * (nsim - 1)))
  
  MSE <- mean((estimates - theta)^2)
  MSE_mcse <- sqrt(sum(((estimates - theta)^2 - MSE)^2) /
                     (nsim * (nsim - 1)))
  
  cover <- mean(ci_low <= theta & ci_high >= theta)
  cover_mcse <- sqrt(cover * (1 - cover) / nsim)
  
  # Model-based SE from CI width
  s2 <- se_mod^2
  m_s2 <- mean(s2, na.rm = TRUE)
  ModSE <- sqrt(m_s2)
  
  # Delta-method MC SE for ModSE
  var_m <- var(s2, na.rm = TRUE) / nsim
  ModSE_mcse <- if (m_s2 > 0) sqrt(var_m / (4 * m_s2)) else NA_real_
  
  # Bias-eliminated coverage (coverage of CI for mean(est))
  theta_hat <- mean(estimates)
  cover_be  <- mean(ci_low <= theta_hat & ci_high >= theta_hat)
  cover_be_mcse <- sqrt(cover_be * (1 - cover_be) / nsim)
  
  # Relative error in ModSE vs empirical SE
  rel_ModSE <- ModSE / empSE - 1
  # Approx MC SE using independence assumption between ModSE and empSE
  rel_ModSE_mcse <- sqrt(
    (ModSE_mcse / empSE)^2 +
      (ModSE * empSE_mcse / empSE^2)^2
  )
  
  data.frame(
    bias, bias_mcse,
    empSE, empSE_mcse,
    MSE, MSE_mcse,
    cover, cover_mcse,
    ModSE, ModSE_mcse,
    cover_be, cover_be_mcse,
    rel_ModSE, rel_ModSE_mcse
  )
}

###############################################################################
# 2. Load SEM results (per-replicate, 1000 runs)
###############################################################################

sem_results_file <- file.path(sem_dir, "sem_results_all_1000.rds")
if (!file.exists(sem_results_file)) {
  stop("SEM results file not found: ", sem_results_file)
}

sem_results <- readRDS(sem_results_file)

# Expected columns in sem_results:
# scenario, rep_id, beta_xy_hat, beta_xy_low, beta_xy_high, beta_xy_true,
# beta_yx_hat, beta_yx_low, beta_yx_high, beta_yx_true

# Long format: one row per method × scenario × rep × effect
sem_long <- bind_rows(
  sem_results %>%
    transmute(
      method = "SEM",
      scenario,
      rep_id,
      effect = "XY",
      est    = beta_xy_hat,
      ci_low = beta_xy_low,
      ci_high= beta_xy_high,
      theta  = beta_xy_true
    ),
  sem_results %>%
    transmute(
      method = "SEM",
      scenario,
      rep_id,
      effect = "YX",
      est    = beta_yx_hat,
      ci_low = beta_yx_low,
      ci_high= beta_yx_high,
      theta  = beta_yx_true
    )
)

###############################################################################
# 3. Load DYNAMITE results (1000 runs, merging checkpoint files)
###############################################################################

# Helper: merge all dynamite_results_sX_*.rds using file modification times
load_dyn_results_for_scenario <- function(s) {
  pattern <- sprintf("^dynamite_results_s%i_.*\\.rds$", s)
  files <- list.files(dyn_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0L) {
    stop("No DYNAMITE results files found for scenario ", s)
  }
  
  file_info <- tibble(
    path  = files,
    mtime = file.info(files)$mtime
  ) %>%
    arrange(mtime)
  
  dyn_list <- purrr::map2(
    file_info$path,
    seq_along(file_info$path),
    ~ {
      df <- readRDS(.x)
      df$file_order <- .y
      df
    }
  )
  
  dyn_all <- bind_rows(dyn_list)
  
  # Keep the latest entry per scenario × replicate (handles re-runs/checkpoints)
  dyn_all %>%
    arrange(scenario, rep_id, file_order) %>%
    group_by(scenario, rep_id) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    select(-file_order)
}

dyn_results <- bind_rows(
  load_dyn_results_for_scenario(1),
  load_dyn_results_for_scenario(2),
  load_dyn_results_for_scenario(3)
)

# Optional sanity check: report number of reps per scenario/method
message("DYNAMITE reps per scenario (after merge):")
print(
  dyn_results %>%
    count(scenario) %>%
    arrange(scenario)
)

# DYNAMITE to long format
dyn_long <- bind_rows(
  dyn_results %>%
    transmute(
      method = "DYNAMITE",
      scenario,
      rep_id,
      effect = "XY",
      est    = beta_xy_hat,
      ci_low = beta_xy_q05,
      ci_high= beta_xy_q95,
      theta  = beta_xy_true
    ),
  dyn_results %>%
    transmute(
      method = "DYNAMITE",
      scenario,
      rep_id,
      effect = "YX",
      est    = beta_yx_hat,
      ci_low = beta_yx_q05,
      ci_high= beta_yx_q95,
      theta  = beta_yx_true
    )
)

###############################################################################
# 4. Combined long dataset for plotting
###############################################################################

all_long <- bind_rows(sem_long, dyn_long) %>%
  mutate(
    scenario_label = recode(as.character(scenario), !!!scenario_labels),
    method = factor(method, levels = c("SEM", "DYNAMITE")),
    effect = factor(effect, levels = c("XY", "YX")),
    # Model-based SE from CI width (assuming 95% equal-tailed intervals)
    se_mod = (ci_high - ci_low) / (2 * qnorm(0.975)),
    covered = ci_low <= theta & ci_high >= theta
  )

saveRDS(all_long, file.path(combined_dir, "all_methods_long_results_1000.rds"))
readr::write_csv(all_long,
                 file.path(combined_dir, "all_methods_long_results_1000.csv"))

###############################################################################
# 5. Performance metrics per method × scenario × effect (Morris-style)
###############################################################################

# Use data.table for clean MC performance calculation
all_dt <- as.data.table(all_long)

perf_dt <- all_dt[, mc_se(est, unique(theta), ci_low, ci_high, se_mod),
                  by = .(method, scenario, effect)]

# Add labels + RMSE and its MC SE
perf_dt[, `:=`(
  scenario_label = scenario_labels[as.character(scenario)],
  RMSE = sqrt(MSE),
  RMSE_mcse = ifelse(MSE > 0, 0.5 * MSE_mcse / sqrt(MSE), NA_real_)
)]

setcolorder(
  perf_dt,
  c("method","scenario","scenario_label","effect",
    "bias","bias_mcse",
    "empSE","empSE_mcse",
    "MSE","MSE_mcse","RMSE","RMSE_mcse",
    "cover","cover_mcse",
    "ModSE","ModSE_mcse",
    "cover_be","cover_be_mcse",
    "rel_ModSE","rel_ModSE_mcse")
)

perf_all <- as_tibble(perf_dt)

# Split out per-method tables
sem_performance  <- perf_all %>% filter(method == "SEM")
dyn_performance  <- perf_all %>% filter(method == "DYNAMITE")

# Persist tables (CSV & RDS)
readr::write_csv(sem_performance,
                 file.path(tables_dir, "sem_performance_by_scenario_1000.csv"))
readr::write_csv(dyn_performance,
                 file.path(tables_dir, "dynamite_performance_by_scenario_1000.csv"))
readr::write_csv(perf_all,
                 file.path(combined_dir, "combined_performance_by_scenario_1000.csv"))

saveRDS(sem_performance,
        file.path(tables_dir, "sem_performance_by_scenario_1000.rds"))
saveRDS(dyn_performance,
        file.path(tables_dir, "dynamite_performance_by_scenario_1000.rds"))
saveRDS(perf_all,
        file.path(combined_dir, "combined_performance_by_scenario_1000.rds"))

###############################################################################
# 6. Plot factory: standard plots + ZIP + lollipop + SE diagnostics
###############################################################################

# Helper to save both PNG and PDF
save_plot_dual <- function(plot, filename_base, width = 8, height = 4) {
  png_file <- file.path(plots_dir, paste0(filename_base, ".png"))
  pdf_file <- file.path(plots_dir, paste0(filename_base, ".pdf"))
  ggsave(png_file, plot, width = width, height = height, dpi = 300)
  ggsave(pdf_file, plot, width = width, height = height)
}

# Iterate over scenarios
for (s in sort(unique(all_long$scenario))) {
  
  scen_data <- all_long %>% filter(scenario == s)
  scen_label <- unique(scen_data$scenario_label)
  
  # 6.1 Density of estimates by method, faceted by effect --------------------
  p_density <- ggplot(scen_data,
                      aes(x = est, colour = method, fill = method)) +
    geom_density(alpha = 0.25) +
    geom_vline(aes(xintercept = theta),
               linetype = "dashed", colour = "black") +
    facet_wrap(~effect, scales = "free", nrow = 1) +
    labs(
      title = paste0("Distribution of estimates: ", scen_label),
      x = "Estimate", y = "Density",
      colour = "Method", fill = "Method"
    ) +
    theme_bw() +
    scale_colour_manual(values = morris_cols) +
    scale_fill_manual(values = morris_fill)
  
  save_plot_dual(p_density,
                 sprintf("s%i_density_estimates_1000", s),
                 width = 9, height = 4)
  
  # 6.2 Boxplots of estimates per method × effect ----------------------------
  p_box <- ggplot(scen_data,
                  aes(x = method, y = est, fill = method)) +
    geom_boxplot(outlier.alpha = 0.4) +
    geom_hline(aes(yintercept = theta),
               linetype = "dashed", colour = "black") +
    facet_wrap(~effect, scales = "free_y", nrow = 1) +
    labs(
      title = paste0("Boxplots of estimates: ", scen_label),
      x = "Method", y = "Estimate"
    ) +
    theme_bw() +
    scale_fill_manual(values = morris_fill)
  
  save_plot_dual(p_box,
                 sprintf("s%i_boxplots_estimates_1000", s),
                 width = 7, height = 4)
  
  # 6.3 Estimates by replicate (jittered) ------------------------------------
  p_scatter <- ggplot(scen_data,
                      aes(x = rep_id, y = est, colour = method)) +
    geom_point(alpha = 0.6, size = 0.8,
               position = position_jitter(width = 0.15, height = 0)) +
    geom_hline(aes(yintercept = theta),
               linetype = "dashed", colour = "black") +
    facet_wrap(~effect, scales = "free_y", nrow = 2) +
    labs(
      title = paste0("Estimates across replicates: ", scen_label),
      x = "Replication", y = "Estimate",
      colour = "Method"
    ) +
    theme_bw() +
    scale_colour_manual(values = morris_cols)
  
  save_plot_dual(p_scatter,
                 sprintf("s%i_scatter_estimates_by_rep_1000", s),
                 width = 9, height = 5)
  
  # 6.4 Performance barplots: bias, RMSE, coverage ---------------------------
  scen_perf <- perf_all %>% filter(scenario == s)
  
  # Bias
  p_bias <- ggplot(scen_perf,
                   aes(x = method, y = bias, fill = method)) +
    geom_col(position = position_dodge(width = 0.7)) +
    geom_errorbar(aes(ymin = bias - bias_mcse,
                      ymax = bias + bias_mcse),
                  width = 0.2, position = position_dodge(width = 0.7)) +
    facet_wrap(~effect, nrow = 1) +
    labs(
      title = paste0("Bias by method: ", scen_label),
      x = "Method", y = "Bias"
    ) +
    theme_bw() +
    scale_fill_manual(values = morris_fill)
  
  save_plot_dual(p_bias,
                 sprintf("s%i_bias_by_method_1000", s),
                 width = 7, height = 4)
  
  # RMSE
  p_rmse <- ggplot(scen_perf,
                   aes(x = method, y = RMSE, fill = method)) +
    geom_col(position = position_dodge(width = 0.7)) +
    geom_errorbar(aes(ymin = RMSE - RMSE_mcse,
                      ymax = RMSE + RMSE_mcse),
                  width = 0.2, position = position_dodge(width = 0.7)) +
    facet_wrap(~effect, nrow = 1, scales = "free_y") +
    labs(
      title = paste0("RMSE by method: ", scen_label),
      x = "Method", y = "RMSE"
    ) +
    theme_bw() +
    scale_fill_manual(values = morris_fill)
  
  save_plot_dual(p_rmse,
                 sprintf("s%i_rmse_by_method_1000", s),
                 width = 7, height = 4)
  
  # Coverage
  p_cover <- ggplot(scen_perf,
                    aes(x = method, y = cover, fill = method)) +
    geom_col(position = position_dodge(width = 0.7)) +
    geom_errorbar(aes(ymin = cover - cover_mcse,
                      ymax = cover + cover_mcse),
                  width = 0.2, position = position_dodge(width = 0.7)) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    scale_y_continuous(limits = c(0, 1)) +
    facet_wrap(~effect, nrow = 1) +
    labs(
      title = paste0("Coverage by method: ", scen_label),
      x = "Method", y = "Coverage"
    ) +
    theme_bw() +
    scale_fill_manual(values = morris_fill)
  
  save_plot_dual(p_cover,
                 sprintf("s%i_coverage_by_method_1000", s),
                 width = 7, height = 4)
  
  # 6.5 ZIP plots (per effect) -----------------------------------------------
  for (eff in levels(scen_data$effect)) {
    
    scen_eff <- scen_data %>% filter(effect == eff)
    
    zip_data <- scen_eff %>%
      filter(!is.na(se_mod), se_mod > 0) %>%
      mutate(
        z = (est - theta) / se_mod,
        rank = rank(abs(z))
      )
    
    p_zip <- ggplot(zip_data,
                    aes(x = rank, y = z, colour = method)) +
      geom_point(size = 1.2, alpha = 0.7) +
      geom_hline(yintercept = c(-1.96, 1.96),
                 linetype = "dashed") +
      labs(
        title = paste0("ZIP plot: ", scen_label, " (Effect ", eff, ")"),
        x = "Rank by |z|",
        y = "z = (est - theta) / SE(model)"
      ) +
      theme_bw() +
      scale_colour_manual(values = morris_cols)
    
    save_plot_dual(
      p_zip,
      sprintf("s%i_zip_%s_1000", s, eff),
      width = 7, height = 4
    )
  }
  
  # 6.6 Lollipop performance summary (per scenario) --------------------------
  perf_scen <- perf_all %>% filter(scenario == s)
  
  perf_long <- perf_scen %>%
    select(method, effect,
           bias, bias_mcse,
           RMSE, RMSE_mcse,
           cover, cover_mcse,
           cover_be, cover_be_mcse,
           ModSE, ModSE_mcse,
           rel_ModSE, rel_ModSE_mcse) %>%
    pivot_longer(
      cols = c(bias, RMSE, cover, cover_be, ModSE, rel_ModSE),
      names_to = "metric",
      values_to = "estimate"
    ) %>%
    mutate(
      mcse = case_when(
        metric == "bias"      ~ bias_mcse,
        metric == "RMSE"      ~ RMSE_mcse,
        metric == "cover"     ~ cover_mcse,
        metric == "cover_be"  ~ cover_be_mcse,
        metric == "ModSE"     ~ ModSE_mcse,
        metric == "rel_ModSE" ~ rel_ModSE_mcse,
        TRUE ~ NA_real_
      ),
      metric = factor(
        metric,
        levels = c("bias","RMSE","cover","cover_be","ModSE","rel_ModSE"),
        labels = c("Bias","RMSE","Coverage","Bias-elim coverage",
                   "Model-based SE","Rel. error ModSE")
      )
    )
  
  p_lolli <- ggplot(perf_long,
                    aes(x = metric, y = estimate, colour = method)) +
    geom_point(position = position_dodge(width = 0.6), size = 2) +
    geom_errorbar(aes(ymin = estimate - mcse,
                      ymax = estimate + mcse),
                  position = position_dodge(width = 0.6),
                  width = 0) +
    coord_flip() +
    facet_wrap(~effect, scales = "free_x") +
    theme_bw() +
    labs(
      title = paste0("Performance summary (lollipop): ", scen_label),
      x = "Metric",
      y = "Value"
    ) +
    scale_colour_manual(values = morris_cols)
  
  save_plot_dual(
    p_lolli,
    sprintf("s%i_lollipop_summary_1000", s),
    width = 10, height = 6
  )
  
  # 6.7 SE(model) density ----------------------------------------------------
  p_se_density <- ggplot(scen_data,
                         aes(x = se_mod, colour = method, fill = method)) +
    geom_density(alpha = 0.25) +
    facet_wrap(~effect, nrow = 1, scales = "free") +
    theme_bw() +
    labs(title = paste0("Model-based SE distribution: ", scen_label),
         x = "SE(model)", y = "Density") +
    scale_colour_manual(values = morris_cols) +
    scale_fill_manual(values = morris_fill)
  
  save_plot_dual(p_se_density,
                 sprintf("s%i_se_density_1000", s),
                 width = 8, height = 4)
  
  # 6.8 SE(model) vs estimate -----------------------------------------------
  p_se_vs_est <- ggplot(scen_data,
                        aes(x = est, y = se_mod, colour = method)) +
    geom_point(alpha = 0.5) +
    facet_wrap(~effect, scales = "free") +
    theme_bw() +
    labs(title = paste0("SE(model) vs estimate: ", scen_label),
         x = "Estimate", y = "SE(model)") +
    scale_colour_manual(values = morris_cols)
  
  save_plot_dual(p_se_vs_est,
                 sprintf("s%i_se_vs_est_1000", s),
                 width = 8, height = 4)
  
  # 6.9 SEM vs DYNAMITE limits-of-agreement ----------------------------------
  wide <- scen_data %>%
    select(method, rep_id, effect, est) %>%
    pivot_wider(names_from = method, values_from = est) %>%
    filter(!is.na(SEM), !is.na(DYNAMITE))
  
  if (nrow(wide) > 0) {
    p_loa <- ggplot(wide,
                    aes(x = (SEM + DYNAMITE)/2,
                        y = SEM - DYNAMITE)) +
      geom_point(alpha = 0.5) +
      geom_hline(yintercept = mean(wide$SEM - wide$DYNAMITE),
                 linetype = "dashed") +
      facet_wrap(~effect, scales = "free") +
      labs(
        title = paste0("Limits of agreement: ", scen_label),
        x = "Mean(SEM, DYNAMITE)",
        y = "SEM − DYNAMITE"
      ) +
      theme_bw()
    
    save_plot_dual(p_loa,
                   sprintf("s%i_limits_of_agreement_1000", s),
                   width = 8, height = 4)
  }
}

###############################################################################
# 7. Done ---------------------------------------------------------------------

message("Performance tables and plots (1000 runs) written under: ", perf_root)
