#!/usr/bin/env Rscript
# =====================================================================
# 03_solid_fit.R  -- a properly-sampled fit (not the capped feasibility one)
#
# The capped feasibility run (treedepth=10) gave bad diagnostics, but that was
# an artifact of the cap: the treedepth=13 Stage-2a fit converged cleanly
# (div=2, Rhat=1.008, bulk-ESS=764). So the model is sound; it just needs deep
# NUTS trees (the t0/k/nu ridge on resorption-truncated series). This script
# samples the scientifically central model -- shape (nu) pooled per species --
# at full treedepth, and compares it to the global-nu fit by LOO.
#
# Reuses outputs/fit_stage2_pooled.rds (clean global-nu Stage 2a).
# Runtime warning: rstan + treedepth=13 is ~1-1.5 h here. Off-sandbox use
# cmdstanr (much faster). Saves incrementally.
# =====================================================================

suppressMessages({ library(brms); library(dplyr); library(readr); library(ggplot2) })
options(mc.cores = 1, brms.backend = "rstan")
set.seed(1)
root <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
clean <- file.path(root, "data_clean"); out <- file.path(root, "outputs")

con <- file(file.path(out, "solid_fit_report.txt"), open = "wt")
sink(con, split = TRUE); on.exit({ sink(); close(con) }, add = TRUE)
rule <- function(s) cat("\n", strrep("=", 70), "\n", s, "\n", strrep("=", 70), "\n", sep = "")

dat <- read_csv(file.path(clean, "lag_femur_long.csv"), show_col_types = FALSE)
nlexpr <- circ_mm ~ logA - exp(-lognu) * log1p(exp(lognu) * exp(-exp(logk) * (age - t0)))

# known-good config (the one that gave clean Stage 2a)
CTRL <- list(adapt_delta = 0.97, max_treedepth = 13)
ITER <- 1500; WARM <- 750; CH <- 4

diag_block <- function(fit, label) {
  rule(paste("DIAGNOSTICS:", label)); print(summary(fit))
  np <- brms::nuts_params(fit)
  cat(sprintf(">> divergences: %d | max treedepth hits: %d | max Rhat: %.3f | min bulk-ESS: %.0f\n",
    sum(subset(np, Parameter == "divergent__")$Value),
    sum(subset(np, Parameter == "treedepth__")$Value >= 13),
    max(brms::rhat(fit), na.rm = TRUE),
    min(posterior::summarise_draws(as_draws(fit), "ess_bulk")$ess_bulk, na.rm = TRUE)))
}

fit2a <- readRDS(file.path(out, "fit_stage2_pooled.rds"))   # clean global-nu

# ---- solid Stage 2b: nu ~ 1 + (1|species) -------------------------
rule("SOLID FIT: shape pooled per species  (treedepth=13)")
f2b <- bf(nlexpr,
          logA  ~ 1 + (1 | species) + (1 | specimen),
          logk  ~ 1 + (1 | species),
          lognu ~ 1 + (1 | species),
          t0    ~ 1 + (1 | specimen),
          nl = TRUE)
p2b <- c(prior(normal(5.5, 1),  nlpar = "logA"),
         prior(normal(-1, 1),   nlpar = "logk"),
         prior(normal(0, 0.5),  nlpar = "lognu"),
         prior(normal(0, 2),    nlpar = "t0"),
         prior(student_t(3, 0, 0.7), class = "sd", nlpar = "logA"),
         prior(student_t(3, 0, 0.5), class = "sd", nlpar = "logk"),
         prior(student_t(3, 0, 0.3), class = "sd", nlpar = "lognu"),
         prior(student_t(3, 0, 1),   class = "sd", nlpar = "t0"))
fit2b <- brm(f2b, data = dat, family = lognormal(), prior = p2b,
             chains = CH, iter = ITER, warmup = WARM, init = 0.2, seed = 1,
             control = CTRL, refresh = 250)
saveRDS(fit2b, file.path(out, "fit_solid_2b.rds"))
diag_block(fit2b, "solid 2b (nu by species)")

# ---- LOO comparison: global-nu vs nu-by-species -------------------
tryCatch({
  rule("MODEL COMPARISON (LOO): global nu vs nu-by-species")
  l2a <- loo(fit2a); l2b <- loo(fit2b)
  print(loo_compare(l2a, l2b))
  cat("\n(positive elpd_diff favors the first row; SE in 2nd column)\n")
}, error = function(e) cat("[loo skipped] ", conditionMessage(e), "\n"))

# ---- posterior predictive check -----------------------------------
tryCatch({
  ggsave(file.path(out, "fig_solid_ppcheck.png"),
         pp_check(fit2b, ndraws = 80) + ggtitle("PP check: solid nu-by-species"),
         width = 6, height = 4, dpi = 120)
}, error = function(e) cat("[ppcheck skipped] ", conditionMessage(e), "\n"))

# ---- nu shrinkage figure (solid) ----------------------------------
tryCatch({
  sp <- sort(unique(dat$species))
  d1 <- as_draws_df(readRDS(file.path(out, "fit_stage1_trex.rds")))
  d2a <- as_draws_df(fit2a); d2b <- as_draws_df(fit2b)
  df <- rbind(data.frame(group = "T. rex only (1 series)", nu = exp(d1$b_lognu_Intercept)),
              data.frame(group = "GLOBAL pooled (2a)",      nu = exp(d2a$b_lognu_Intercept)))
  for (s in sp) { c0 <- paste0("r_species__lognu[", s, ",Intercept]")
    if (c0 %in% names(d2b)) df <- rbind(df, data.frame(group = paste0(s, " (2b solid)"),
                                                       nu = exp(d2b$b_lognu_Intercept + d2b[[c0]]))) }
  sm <- df |> group_by(group) |>
    summarise(m = median(nu), lo50 = quantile(nu, .25), hi50 = quantile(nu, .75),
              lo90 = quantile(nu, .05), hi90 = quantile(nu, .95), .groups = "drop")
  ggsave(file.path(out, "fig_solid_nu.png"),
    ggplot(sm, aes(y = group)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
      geom_linerange(aes(xmin = lo90, xmax = hi90)) +
      geom_linerange(aes(xmin = lo50, xmax = hi50), linewidth = 1.6, colour = "purple4") +
      geom_point(aes(x = m), size = 2.4) + coord_cartesian(xlim = c(0, quantile(df$nu, .98))) +
      labs(title = expression(paste("Solid fit -- shape ", nu)),
           subtitle = "median, 50% & 90% CrI; dashed = 1 (logistic)", x = expression(nu), y = NULL) +
      theme_minimal(base_size = 10), width = 7, height = 5, dpi = 120)
  cat("\nFigures: fig_solid_nu.png, fig_solid_ppcheck.png\n")
}, error = function(e) cat("[nu fig skipped] ", conditionMessage(e), "\n"))

rule("SOLID FIT COMPLETE")
