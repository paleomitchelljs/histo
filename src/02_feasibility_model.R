#!/usr/bin/env Rscript
# =====================================================================
# 02_feasibility_model.R
# Feasibility run for the hierarchical generalized-logistic (Richards)
# growth model of plan sec 1.1-1.5. Staged exactly as plan sec 1.7.
#
# Response : femur diaphyseal circumference (mm), modeled on the log
#            scale via family = lognormal() so variance grows with size
#            and fits stay positive.
#
# Richards (nu > 0) on the log-mean scale (numerically stable, log1p):
#   log M(t) = logA - exp(-lognu) * log1p( exp(lognu) * exp(-exp(logk)*(age - t0)) )
#   <=>  M(t) = A * (1 + nu*exp(-k*(t - t0)))^(-1/nu)
#   with A=exp(logA), k=exp(logk), nu=exp(lognu).
#
# Stages:
#   1. single best-sampled taxon (T. rex FMNH PR 2081 "Sue"): proves the
#      likelihood + whether nu is estimable from one good series.
#   2a. all 7 taxa, partial pooling on logA/logk/t0; single global nu
#       (nu borrows strength across the whole clade -- the core payoff).
#   2b. allow nu to vary by species (1|species) -- does shape shrink?
#   3. (stretch) add phylogenetic covariance on logA (plan Option A).
#
# Outputs -> outputs/ : fitted objects (.rds), feasibility_report.txt,
#            and diagnostic figures (.png).
# =====================================================================

suppressMessages({
  library(brms); library(dplyr); library(readr); library(ape); library(ggplot2)
})
# NOTE: cmdstanr's processx subprocess helper is blocked in this sandbox, so we
# use the rstan backend (compiles + samples in-process). cores=1 -> sequential
# chains, avoiding any process spawning. Switch back to cmdstanr off-sandbox.
options(mc.cores = 1, brms.backend = "rstan")
set.seed(1)

root  <- "/Users/jmitchell/Library/CloudStorage/Dropbox/Research/00_current_Histology"
clean <- file.path(root, "data_clean")
out   <- file.path(root, "outputs"); dir.create(out, showWarnings = FALSE)

con <- file(file.path(out, "feasibility_report.txt"), open = "wt")
sink(con, split = TRUE); on.exit({ sink(); close(con) }, add = TRUE)
rule <- function(s) cat("\n", strrep("=", 70), "\n", s, "\n", strrep("=", 70), "\n", sep = "")

dat <- read_csv(file.path(clean, "lag_femur_long.csv"), show_col_types = FALSE)
tip <- read_csv(file.path(clean, "species_to_tip.csv"), show_col_types = FALSE)
dat <- left_join(dat, tip, by = "species")

ITER <- 1500; WARM <- 750; CHAINS <- 4
ADAPT <- 0.97; TREEDEPTH <- 13

diag_block <- function(fit, label) {
  rule(paste("DIAGNOSTICS:", label))
  print(summary(fit))
  np <- brms::nuts_params(fit)
  ndiv <- sum(subset(np, Parameter == "divergent__")$Value)
  rh   <- brms::rhat(fit); maxrh <- max(rh, na.rm = TRUE)
  cat(sprintf("\n>> divergences: %d | max Rhat: %.3f | min bulk-ESS: %.0f\n",
              ndiv, maxrh,
              min(posterior::summarise_draws(as_draws(fit), "ess_bulk")$ess_bulk,
                  na.rm = TRUE)))
  invisible(ndiv)
}

# Richards on log-mean scale (shared across stages)
nlexpr <- circ_mm ~ logA - exp(-lognu) * log1p(exp(lognu) * exp(-exp(logk) * (age - t0)))

logAc <- max(log(dat$circ_mm))   # prior center for log-asymptote (~log max circ)

# ---------------------------------------------------------------------
# STAGE 1 : single taxon, no pooling ----------------------------------
# ---------------------------------------------------------------------
rule("STAGE 1: Tyrannosaurus rex FMNH PR 2081 ('Sue'), single series")
d1 <- filter(dat, specimen == "FMNH_PR_2081")
cat(sprintf("n = %d LAGs, circ %.0f-%.0f mm, EFS present: %s\n",
            nrow(d1), min(d1$circ_mm), max(d1$circ_mm), any(d1$is_efs)))

f1 <- bf(nlexpr, logA ~ 1, logk ~ 1, lognu ~ 1, t0 ~ 1, nl = TRUE)
p1 <- c(prior(normal(6.3, 0.5), nlpar = "logA"),
        prior(normal(-1, 1),    nlpar = "logk"),
        prior(normal(0, 0.5),   nlpar = "lognu"),
        prior(normal(0, 5),     nlpar = "t0"))

fit1 <- brm(f1, data = d1, family = lognormal(), prior = p1,
            chains = CHAINS, iter = ITER, warmup = WARM, init = 0.2,
            control = list(adapt_delta = ADAPT, max_treedepth = TREEDEPTH),
            refresh = 0, silent = 2, seed = 1)
saveRDS(fit1, file.path(out, "fit_stage1_trex.rds"))
diag_block(fit1, "Stage 1 (Sue only)")

# ---------------------------------------------------------------------
# STAGE 2a : all taxa, partial pooling, single global nu --------------
# ---------------------------------------------------------------------
rule("STAGE 2a: all 7 taxa | partial pooling on logA/logk/t0 | global nu")
f2 <- bf(nlexpr,
         logA  ~ 1 + (1 | species) + (1 | specimen),
         logk  ~ 1 + (1 | species),
         lognu ~ 1,                                   # one shape, all-clade-informed
         t0    ~ 1 + (1 | specimen),
         nl = TRUE)
p2 <- c(prior(normal(5.5, 1),  nlpar = "logA"),
        prior(normal(-1, 1),   nlpar = "logk"),
        prior(normal(0, 0.5),  nlpar = "lognu"),
        prior(normal(0, 5),    nlpar = "t0"),
        prior(student_t(3, 0, 1),   class = "sd", nlpar = "logA"),
        prior(student_t(3, 0, 1),   class = "sd", nlpar = "logk"),
        prior(student_t(3, 0, 2),   class = "sd", nlpar = "t0"))

fit2 <- brm(f2, data = dat, family = lognormal(), prior = p2,
            chains = CHAINS, iter = ITER, warmup = WARM, init = 0.2,
            control = list(adapt_delta = ADAPT, max_treedepth = TREEDEPTH),
            refresh = 0, silent = 2, seed = 1)
saveRDS(fit2, file.path(out, "fit_stage2_pooled.rds"))
diag_block(fit2, "Stage 2a (pooled, global nu)")

# ---------------------------------------------------------------------
# STAGE 2b : let nu vary by species -----------------------------------
# ---------------------------------------------------------------------
rule("STAGE 2b: same, but lognu ~ 1 + (1 | species) -- does shape shrink?")
f2b <- bf(nlexpr,
          logA  ~ 1 + (1 | species) + (1 | specimen),
          logk  ~ 1 + (1 | species),
          lognu ~ 1 + (1 | species),
          t0    ~ 1 + (1 | specimen),
          nl = TRUE)
p2b <- c(p2, prior(student_t(3, 0, 0.5), class = "sd", nlpar = "lognu"))
fit2b <- brm(f2b, data = dat, family = lognormal(), prior = p2b,
             chains = CHAINS, iter = ITER, warmup = WARM, init = 0.2,
             control = list(adapt_delta = ADAPT, max_treedepth = TREEDEPTH),
             refresh = 0, silent = 2, seed = 1)
saveRDS(fit2b, file.path(out, "fit_stage2b_nu_by_species.rds"))
diag_block(fit2b, "Stage 2b (nu by species)")

# ---------------------------------------------------------------------
# STAGE 3 : (stretch) topology-only phylogenetic covariance on logA ---
# ---------------------------------------------------------------------
ok3 <- tryCatch({
  rule("STAGE 3 (stretch): + (1 | gr(phylo_tip, cov = A_phylo)) on logA")
  tree <- read.tree(file.path(clean, "theropod_femur_tree.nwk"))
  tree_unit <- compute.brlen(tree, 1)                 # topology only
  A_phylo <- vcv(tree_unit, corr = TRUE)
  f3 <- bf(nlexpr,
           logA  ~ 1 + (1 | species) + (1 | specimen) + (1 | gr(phylo_tip, cov = A_phylo)),
           logk  ~ 1 + (1 | species),
           lognu ~ 1,
           t0    ~ 1 + (1 | specimen),
           nl = TRUE)
  fit3 <- brm(f3, data = dat, family = lognormal(), prior = p2,
              data2 = list(A_phylo = A_phylo),
              chains = CHAINS, iter = ITER, warmup = WARM, init = 0.2,
              control = list(adapt_delta = 0.99, max_treedepth = 14),
              refresh = 0, silent = 2, seed = 1)
  saveRDS(fit3, file.path(out, "fit_stage3_phylo.rds"))
  diag_block(fit3, "Stage 3 (phylo covariance)")
  TRUE
}, error = function(e) { cat("\n[Stage 3 did not complete] ", conditionMessage(e), "\n"); FALSE })

# ---------------------------------------------------------------------
# FIGURES -------------------------------------------------------------
# ---------------------------------------------------------------------
tryCatch({
  ag <- seq(min(dat$age), max(dat$age) + 2, length.out = 80)

  # Fig 1: Stage 1 fit over Sue
  nd1 <- data.frame(age = ag)
  ep1 <- posterior_epred(fit1, newdata = nd1)
  g1 <- data.frame(age = ag,
                   med = apply(ep1, 2, median),
                   lo  = apply(ep1, 2, quantile, .05),
                   hi  = apply(ep1, 2, quantile, .95))
  p_fig1 <- ggplot() +
    geom_ribbon(data = g1, aes(age, ymin = lo, ymax = hi), alpha = .25, fill = "steelblue") +
    geom_line(data = g1, aes(age, med), colour = "steelblue", linewidth = 1) +
    geom_point(data = d1, aes(age, circ_mm, shape = is_efs), size = 2) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), name = "EFS") +
    labs(title = "Stage 1: Richards fit, T. rex FMNH PR 2081 (Sue)",
         x = "LAG index (age proxy)", y = "femur circumference (mm)") +
    theme_minimal(base_size = 12)
  ggsave(file.path(out, "fig_stage1_trex_fit.png"), p_fig1, width = 6, height = 4.2, dpi = 120)

  # Fig 2: Stage 2a species-level curves over data
  sp <- sort(unique(dat$species))
  nd2 <- expand.grid(age = ag, species = sp, specimen = "._dummy_",
                     phylo_tip = "._dummy_", stringsAsFactors = FALSE)
  ep2 <- posterior_epred(fit2, newdata = nd2,
                         re_formula = ~ (1 | species), allow_new_levels = TRUE)
  nd2$med <- apply(ep2, 2, median)
  nd2$lo  <- apply(ep2, 2, quantile, .05)
  nd2$hi  <- apply(ep2, 2, quantile, .95)
  p_fig2 <- ggplot() +
    geom_ribbon(data = nd2, aes(age, ymin = lo, ymax = hi), alpha = .2, fill = "darkorange") +
    geom_line(data = nd2, aes(age, med), colour = "darkorange") +
    geom_point(data = dat, aes(age, circ_mm), size = 1) +
    facet_wrap(~ species, scales = "free_y") +
    labs(title = "Stage 2a: partially-pooled species growth curves",
         x = "LAG index (age proxy)", y = "femur circumference (mm)") +
    theme_minimal(base_size = 10)
  ggsave(file.path(out, "fig_stage2_species_curves.png"), p_fig2, width = 9, height = 6, dpi = 120)

  # Fig 3: nu posteriors -- global (2a) vs per-species (2b)
  d2a <- as_draws_df(fit2)
  nu_global <- exp(d2a$b_lognu_Intercept)
  re <- ranef(fit2b)$species[,, "lognu_Intercept"]   # posterior summary per species
  d2b <- as_draws_df(fit2b)
  sp2 <- sort(unique(dat$species))
  nu_sp <- do.call(rbind, lapply(sp2, function(s){
    col <- paste0("r_species__lognu[", s, ",Intercept]")
    if (!col %in% names(d2b)) return(NULL)
    data.frame(species = s, nu = exp(d2b$b_lognu_Intercept + d2b[[col]]))
  }))
  df_nu <- rbind(
    data.frame(group = "GLOBAL (Stage 2a)", nu = nu_global),
    transform(nu_sp, group = paste0(species, " (2b)"))[, c("group","nu")]
  )
  # interval summary (no ggdist dependency): point = median, thick = 50%, thin = 90%
  summ_nu <- df_nu %>% group_by(group) %>%
    summarise(m = median(nu),
              lo50 = quantile(nu, .25), hi50 = quantile(nu, .75),
              lo90 = quantile(nu, .05), hi90 = quantile(nu, .95), .groups = "drop")
  p_fig3 <- ggplot(summ_nu, aes(y = group)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
    geom_linerange(aes(xmin = lo90, xmax = hi90)) +
    geom_linerange(aes(xmin = lo50, xmax = hi50), linewidth = 1.6, colour = "purple4") +
    geom_point(aes(x = m), size = 2.4) +
    coord_cartesian(xlim = c(0, quantile(df_nu$nu, .98))) +
    labs(title = expression(paste("Shape parameter ", nu, ": pooled vs per-species")),
         subtitle = "median, 50% & 90% CrI; dashed = 1 (logistic); nu<1 Gompertz-ward, nu>1 super-logistic",
         x = expression(nu), y = NULL) +
    theme_minimal(base_size = 10)
  ggsave(file.path(out, "fig_nu_posteriors.png"), p_fig3, width = 7, height = 5, dpi = 120)
  cat("\nFigures written to outputs/\n")
}, error = function(e) cat("\n[figure step warning] ", conditionMessage(e), "\n"))

# ---------------------------------------------------------------------
# nu summary table ----------------------------------------------------
# ---------------------------------------------------------------------
rule("SHAPE PARAMETER nu = exp(lognu) -- posterior summaries")
nu_q <- function(x) sprintf("%.2f [%.2f, %.2f]", median(x),
                            quantile(x, .05), quantile(x, .95))
cat("Stage 1 (Sue, single series):       nu =",
    nu_q(exp(as_draws_df(fit1)$b_lognu_Intercept)), "\n")
cat("Stage 2a (global, all-clade-pooled): nu =",
    nu_q(exp(as_draws_df(fit2)$b_lognu_Intercept)), "\n")

rule("FEASIBILITY RUN COMPLETE")
cat("Fitted objects + figures in outputs/. Stage 3 phylo ran:", ok3, "\n")
