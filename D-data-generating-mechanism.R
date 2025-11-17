# =====================================================================
# Causal Inference Panel Data — Simulation + Visualization Pipeline
# All outputs saved to disk (no printing)
# Working directory fixed to project folder on Charité HPC
# =====================================================================

# =====================================================================
# 0. WORKING DIRECTORY + LIBRARIES
# =====================================================================
setwd("/data/cephfs-1/home/users/tato10_c/work/causal_inference_panel_data")

library(data.table)
library(tidyverse)
library(glue)
library(patchwork)

# =====================================================================
# 1. DIRECTORY STRUCTURE (SAFE, IDEMPOTENT)
# =====================================================================
save_dir <- "data"
meta_dir <- "meta"
fig_dir  <- file.path(meta_dir, "figures")

dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(meta_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,  showWarnings = FALSE, recursive = TRUE)

# =====================================================================
# 2. GLOBAL SETTINGS
# =====================================================================
scenarios <- 1:3

# RNG log constructor
make_rng_log <- function(nsim, base_seed = 20251001L) {
  tibble(
    scenario = rep(scenarios, each = nsim),
    rep_id   = rep(seq_len(nsim), times = length(scenarios)),
    seed     = base_seed + 1e6L * rep(scenarios, each = nsim) + rep_id
  )
}

# Reproducibility metadata
write_session_info <- function(meta_dir, nsim, base_seed) {
  writeLines(c(
    paste("nsim:", nsim),
    paste("base_seed:", base_seed),
    paste("R.version:", R.version.string),
    paste("data.table:", as.character(packageVersion("data.table"))),
    paste("tidyverse:", as.character(packageVersion("tidyverse")))
  ), file.path(meta_dir, "session_info.txt"))
}

# =====================================================================
# 3. DATA-GENERATING MECHANISM
# =====================================================================
simulate_panel <- function(
    scenario = 1,
    N = 1000,
    T_obs = 4,
    phi_x = 0.5,
    phi_y = 0.5,
    drift_c = 1.0,
    beta_xy = NULL,
    beta_yx = NULL,
    gamma_cx = 0.3,
    gamma_cy = 0.4,
    sigma_x = 1,
    sigma_y = 1,
    mu_c0 = 30,
    sd_c0 = 5
) {
  
  # Scenario-specific causal effects
  if (is.null(beta_xy)) beta_xy <- ifelse(scenario >= 2, 0.8, 0)
  if (is.null(beta_yx)) beta_yx <- ifelse(scenario == 3, 0.5, 0)
  
  times <- 0:T_obs
  n_time <- length(times)
  
  C <- X <- Y <- matrix(NA_real_, nrow = N, ncol = n_time)
  
  # Baseline
  C[, 1] <- rnorm(N, mu_c0, sd_c0)
  X[, 1] <- gamma_cx * C[, 1] + rnorm(N, 0, sigma_x)
  Y[, 1] <- gamma_cy * C[, 1] + rnorm(N, 0, sigma_y)
  
  # Positive increments (aging mechanism)
  mean_between_mo <- 2
  mean_within_mo  <- 1
  k_b <- 2; theta_b <- (mean_between_mo/12)/k_b
  k_w <- 2; theta_w <- (mean_within_mo/12)/k_w
  
  b_i <- rgamma(N, shape = k_b, scale = theta_b)
  
  # Follow-ups
  for (t in 2:n_time) {
    e_it <- rgamma(N, shape = k_w, scale = theta_w)
    inc  <- drift_c + b_i + e_it
    
    C[, t] <- C[, t - 1] + inc
    
    X[, t] <- phi_x * X[, t - 1] +
      beta_yx * Y[, t - 1] +
      gamma_cx * C[, t - 1] +
      rnorm(N, 0, sigma_x)
    
    Y[, t] <- phi_y * Y[, t - 1] +
      beta_xy * X[, t - 1] +
      gamma_cy * C[, t - 1] +
      rnorm(N, 0, sigma_y)
  }
  
  data.frame(
    id   = rep(seq_len(N), each = n_time),
    time = rep(times, times = N),
    c    = as.vector(t(C)),
    X    = as.vector(t(X)),
    Y    = as.vector(t(Y)),
    scenario = scenario
  )
}

# =====================================================================
# 4. THEMES + COLORS
# =====================================================================
theme_pub <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold"),
    axis.title       = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    legend.position  = "top"
  )

col_c <- "#E69F00"
col_X <- "#0072B2"
col_Y <- "#009E73"

scen_labs <- c(
  "Scenario 1: No causal effect",
  "Scenario 2: X → Y causal effect",
  "Scenario 3: Bidirectional feedback (X ↔ Y)"
)

# =====================================================================
# 5. PLOTTING FUNCTIONS 
# =====================================================================

# ------------------------------- Trajectories
plot_trajectories <- function(data, scenario_id, n_lines = 30,
                              alpha_lines = 0.4, outfile = NULL) {
  
  sampled_ids <- sample(unique(data$id), min(n_lines, length(unique(data$id))))
  
  plot_data <- data %>%
    filter(id %in% sampled_ids) %>%
    pivot_longer(cols = c(c, X, Y), names_to = "variable", values_to = "value") %>%
    mutate(variable = factor(variable,
                             levels = c("c","X","Y"),
                             labels = c("Confounder","Exposure","Outcome")))
  
  means <- data %>%
    pivot_longer(cols = c(c,X,Y), names_to = "variable", values_to = "value") %>%
    mutate(variable = factor(variable,
                             levels=c("c","X","Y"),
                             labels=c("Confounder","Exposure","Outcome"))) %>%
    group_by(time, variable) %>%
    summarise(mean = mean(value), .groups="drop")
  
  p <- ggplot() +
    geom_line(data = plot_data,
              aes(time, value, group=id, color=variable),
              alpha = alpha_lines, linewidth=0.4) +
    geom_line(data = means,
              aes(time, mean, color=variable),
              linewidth = 1.4) +
    facet_wrap(~variable, scales="free_y") +
    scale_color_manual(values = c("Confounder"=col_c,
                                  "Exposure"=col_X,
                                  "Outcome"=col_Y)) +
    labs(title = scen_labs[scenario_id]) +
    theme_pub +
    theme(legend.position="none")
  
  if (!is.null(outfile)) {
    ggsave(outfile, p, width = 10, height = 4, dpi = 150)
  }
  
  invisible(p)
}

# ------------------------------- Mean trajectories
plot_mean_trajectories <- function(data, scenario_id, outfile=NULL) {
  
  summary_data <- data %>%
    pivot_longer(cols = c(c,X,Y), names_to="variable", values_to="value") %>%
    mutate(variable=factor(variable,
                           levels=c("c","X","Y"),
                           labels=c("Confounder","Exposure","Outcome"))) %>%
    group_by(time, variable) %>%
    summarise(mean=mean(value), sd=sd(value), .groups="drop")
  
  p <- ggplot(summary_data,
              aes(time, mean, color=variable, fill=variable)) +
    geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd), alpha=0.15, color=NA) +
    geom_line(linewidth=1.2) +
    labs(title=scen_labs[scenario_id], subtitle="Mean ± SD") +
    scale_color_manual(values=c("Confounder"=col_c,
                                "Exposure"=col_X,
                                "Outcome"=col_Y)) +
    scale_fill_manual(values=c("Confounder"=col_c,
                               "Exposure"=col_X,
                               "Outcome"=col_Y)) +
    theme_pub
  
  if (!is.null(outfile)) {
    ggsave(outfile, p, width = 9, height = 4.5, dpi = 150)
  }
  
  invisible(p)
}

# ------------------------------- Scatterplot matrix
plot_scatter_matrix <- function(data, scenario_id, outfile=NULL) {
  
  max_time <- max(data$time)
  
  scatter <- data %>%
    filter(time %in% c(0, max_time)) %>%
    mutate(timepoint = ifelse(time==0,"Baseline",
                              paste0("Final (t=",max_time,")"))) %>%
    select(id, timepoint, c, X, Y) %>%
    pivot_longer(cols=c(c,X,Y), names_to="variable", values_to="value")
  
  wide <- scatter %>%
    pivot_wider(names_from = variable, values_from=value)
  
  var_names <- c("c","X","Y")
  pairs <- list()
  
  for (vx in var_names) {
    for (vy in var_names) {
      if (vx != vy) {
        pairs[[length(pairs)+1]] <- wide %>%
          select(id, timepoint, all_of(c(vx,vy))) %>%
          rename(val_x = !!vx, val_y=!!vy) %>%
          mutate(var_x=vx, var_y=vy)
      }
    }
  }
  
  pairs_data <- bind_rows(pairs)
  
  p <- ggplot(pairs_data, aes(val_x, val_y)) +
    geom_point(alpha=0.25, size=1, color="grey40") +
    geom_smooth(method="lm", se=TRUE, color=col_c) +
    facet_grid(var_y ~ var_x + timepoint, scales="free") +
    labs(title=scen_labs[scenario_id]) +
    theme_pub +
    theme(strip.text=element_text(size=9))
  
  if (!is.null(outfile)) {
    ggsave(outfile, p, width=12, height=9, dpi=150)
  }
  
  invisible(p)
}

# =====================================================================
# 6. GENERATE DATA + SAVE ALL FIGURES
# =====================================================================

# Example: N=500 for visualization
data_s1 <- simulate_panel(1, N=500)
data_s2 <- simulate_panel(2, N=500)
data_s3 <- simulate_panel(3, N=500)

# Save figures
plot_trajectories(data_s1,1, outfile=file.path(fig_dir,"S1_trajectories.png"))
plot_mean_trajectories(data_s1,1, outfile=file.path(fig_dir,"S1_means.png"))
plot_scatter_matrix(data_s1,1, outfile=file.path(fig_dir,"S1_scatter.png"))

plot_trajectories(data_s2,2, outfile=file.path(fig_dir,"S2_trajectories.png"))
plot_mean_trajectories(data_s2,2, outfile=file.path(fig_dir,"S2_means.png"))
plot_scatter_matrix(data_s2,2, outfile=file.path(fig_dir,"S2_scatter.png"))

plot_trajectories(data_s3,3, outfile=file.path(fig_dir,"S3_trajectories.png"))
plot_mean_trajectories(data_s3,3, outfile=file.path(fig_dir,"S3_means.png"))
plot_scatter_matrix(data_s3,3, outfile=file.path(fig_dir,"S3_scatter.png"))

# =====================================================================
# 7. FULL DATASET GENERATION PIPELINE (SAVES TO /data)
# =====================================================================

nsim <- 2000
rng_log <- make_rng_log(nsim)

write_csv(rng_log, file.path(meta_dir, "rng_log.csv"))
write_session_info(meta_dir, nsim, 20251001L)

message("Starting dataset generation…")

for (s in 1:3) {
  for (r in 1:nsim) {
    
    outfile <- file.path(save_dir, sprintf("dgm%i_rep%i.rds", s, r))
    
    if (!file.exists(outfile)) {
      seed_r <- rng_log$seed[rng_log$scenario==s & rng_log$rep_id==r]
      set.seed(seed_r)
      
      dat <- simulate_panel(s, N=1000, T_obs=4)
      saveRDS(dat, outfile)
    }
    
    if (r %% 100 == 0) {
      message("Scenario ", s, " — ", r, " / ", nsim)
    }
  }
}

message("Dataset generation complete.")

###############################################################################
# END OF ADEMP — D: DATA
###############################################################################

