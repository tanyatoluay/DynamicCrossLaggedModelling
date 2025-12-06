###############################################################################
# Parametric g-Formula Method (scenario × direction × treatment safe version)
###############################################################################

setwd("/data/cephfs-1/home/users/tato10_c/work/causal_inference_panel_data")

library(gfoRmula)
library(data.table)
library(tidyverse)

##----------------------------------------------------------------------------
## Required globals (must be set before sourcing)
##----------------------------------------------------------------------------

if (!exists("SCENARIO"))        stop("SCENARIO not set.")
if (!exists("TREATMENT_LEVEL")) stop("TREATMENT_LEVEL not set.")
if (!exists("DIRECTION"))       stop("DIRECTION not set. Use 'forward' or 'reverse'.")

scenario         <- SCENARIO
treat_val        <- TREATMENT_LEVEL
direction        <- DIRECTION
nsim             <- if (!exists("NSIM")) 2000L else NSIM
checkpoint_every <- 100L

##----------------------------------------------------------------------------
## Directory structure (new clean architecture)
##----------------------------------------------------------------------------

root_dir  <- file.path(getwd(), "methods", "g-formula")

iter_dir  <- file.path(root_dir, "iterations",
                       paste0("s", scenario), direction,
                       paste0("treat", treat_val))

ckpt_dir  <- file.path(root_dir, "checkpoints_clean",
                       paste0("s", scenario), direction,
                       paste0("treat", treat_val))

final_dir <- file.path(root_dir, "final",
                       paste0("s", scenario), direction)

data_dir  <- file.path(getwd(), "data")

dir.create(iter_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(ckpt_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

##----------------------------------------------------------------------------
## High-level run info
##----------------------------------------------------------------------------

cat("\n----------------------------------------------\n")
cat("Running Parametric g-Formula\n")
cat("Scenario:       ", scenario, "\n")
cat("Direction:      ", direction, "\n")
cat("Treatment:      ", treat_val, "\n")
cat("NSIM:           ", nsim, "\n")
cat("Iteration dir:  ", iter_dir, "\n")
cat("Checkpoint dir: ", ckpt_dir, "\n")
cat("Final dir:      ", final_dir, "\n")
cat("----------------------------------------------\n\n")

##----------------------------------------------------------------------------
## Helpers
##----------------------------------------------------------------------------

iter_file <- function(r) {
  file.path(iter_dir, sprintf("rep%05d.rds", r))
}

ckpt_file <- function(r) {
  file.path(ckpt_dir, sprintf("checkpoint_%05d.rds", r))
}

final_file <- file.path(final_dir,
                        sprintf("treat%s.rds", treat_val))

##----------------------------------------------------------------------------
## Dataset loader
##----------------------------------------------------------------------------

get_dataset <- function(scenario, rep_id) {
  f <- file.path(data_dir, sprintf("dgm%i_rep%i.rds", scenario, rep_id))
  if (!file.exists(f)) stop("Missing dataset: ", f)
  readRDS(f)
}

##----------------------------------------------------------------------------
## Fit one iteration
##----------------------------------------------------------------------------

fit_gformula_once <- function(scenario, rep_id, treat_val, direction) {
  
  dat <- get_dataset(scenario, rep_id)
  setDT(dat)
  K <- max(dat$time)
  
  if (direction == "forward") {
    dat[time < K, Y := NA_real_]
    exposure     <- "X"
    outcome_name <- "Y"
  } else {
    dat[time < K, X := NA_real_]
    exposure     <- "Y"
    outcome_name <- "X"
  }
  
  outcome_type <- "continuous_eof"
  covnames     <- c("c", exposure)
  covtypes     <- c("normal", "normal")
  
  id          <- "id"
  time_name   <- "time"
  time_points <- length(unique(dat[[time_name]]))
  
  histvars  <- list(c(exposure), c("c"))
  histories <- c(lagged, lagged)
  
  covparams <- list(
    covmodels = c(
      c ~ lag1_c + time,
      as.formula(paste0(exposure, " ~ lag1_", exposure, " + lag1_c + c + time"))
    )
  )
  
  intvars       <- list(exposure)
  interventions <- list(list(c(static, rep(treat_val, time_points))))
  int_desc      <- paste0("treat_", treat_val)
  
  ymodel <- as.formula(
    paste0(outcome_name, " ~ ", exposure, " + lag1_", exposure, " + c + lag1_c")
  )
  
  seed_r <- 20260000L + scenario * 100000L + rep_id * 10L +
    as.integer(treat_val * 1000) + ifelse(direction == "forward", 1L, 2L)
  
  res <- gformula(
    obs_data      = dat,
    id            = id,
    time_name     = time_name,
    time_points   = time_points,
    covnames      = covnames,
    covtypes      = covtypes,
    outcome_name  = outcome_name,
    outcome_type  = outcome_type,
    covparams     = covparams,
    ymodel        = ymodel,
    intvars       = intvars,
    interventions = interventions,
    int_descript  = int_desc,
    histories     = histories,
    histvars      = histvars,
    nsimul        = 10000,
    nsamples      = 200,
    seed          = seed_r
  )
  
  out <- as.data.table(res$result)
  out[, scenario  := scenario]
  out[, rep_id    := rep_id]
  out[, treat     := treat_val]
  out[, direction := direction]
  
  out
}

##----------------------------------------------------------------------------
## Resume logic (per scenario × direction × treatment)
##----------------------------------------------------------------------------

existing_ckpts <- list.files(ckpt_dir, pattern = "^checkpoint_\\d+\\.rds$")

if (length(existing_ckpts) == 0L) {
  
  resume_at        <- 0L
  gf_results_list  <- list()
  cat("No checkpoint found. Starting at iteration 1.\n\n")
  
} else {
  
  ck_ids   <- as.integer(gsub("checkpoint_(\\d+)\\.rds", "\\1", existing_ckpts))
  resume_at <- max(ck_ids)
  
  latest_ckpt <- ckpt_file(resume_at)
  cat("Found checkpoint: ", latest_ckpt, "\n", sep = "")
  cat("Resuming at iteration: ", resume_at + 1L, " of ", nsim, "\n\n", sep = "")
  
  prev <- readRDS(latest_ckpt)
  gf_results_list <- as.list(split(prev, seq_len(nrow(prev))))
}

##----------------------------------------------------------------------------
## Main compute loop
##----------------------------------------------------------------------------

if (resume_at < nsim) {
  for (r in (resume_at + 1L):nsim) {
    
    outfile <- iter_file(r)
    if (file.exists(outfile)) next
    
    res <- fit_gformula_once(scenario, r, treat_val, direction)
    
    # Save iteration
    saveRDS(res, outfile)
    
    # Store in memory
    gf_results_list[[length(gf_results_list) + 1L]] <- res
    
    # Checkpoint
    if (r %% checkpoint_every == 0L) {
      ck <- data.table::rbindlist(gf_results_list)
      saveRDS(ck, ckpt_file(r))
      cat("Saved checkpoint: ", ckpt_file(r), "\n", sep = "")
    }
  }
} else {
  cat("All ", nsim, " iterations already completed for this config.\n\n", sep = "")
}

##----------------------------------------------------------------------------
## Final output
##----------------------------------------------------------------------------

final <- data.table::rbindlist(gf_results_list)
saveRDS(final, final_file)

cat("\nCompleted all ", nsim, " iterations.\n", sep = "")
cat("Final output saved to:\n", final_file, "\n\n")
