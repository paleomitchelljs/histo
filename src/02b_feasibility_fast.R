#!/usr/bin/env Rscript
# =====================================================================
# 02b_feasibility_fast.R
# Bounded-runtime rerun of the feasibility stages. The first attempt
# (02_feasibility_model.R) stalled in Stage 2a: max_treedepth=13 on a
# funnel-shaped posterior (the resorption-truncated series weakly
# identify t0/logk/lognu jointly) meant ~2^13 leapfrog steps per
# iteration -> effectively wedged.
#
# Fixes here, all aimed at GUARANTEED termination + clean enough draws:
#   * max_treedepth = 10  (hard cap on cost/iter; saturation reported)
#   * tighter priors on the weak directions (t0 location, group SDs)
#   * fewer iters (1000), progress streamed (refresh)
# Stage 1 is reused from disk (it already sampled cleanly).
# Distinct output names (*_fast) so this can run alongside the old job.
# =====================================================================

suppressMessages({ library(brms); library(dplyr); library(readr); library(ape); library(ggplot2) })
options(mc.cores = 1, brms.backend = "rstan")
set.seed(1)

root  <- "/Users/jmitchell/Library/CloudStorage/Dropbox/Research/00_current_Histology"
clean <- file.path(root, "data_clean"); out <- file.path(root, "outputs")

con <- file(file.path(out, "feasibility_report_fast.txt"), open = "wt")
sink(con, split = TRUE); on.exit({ sink(); close(con) }, add = TRUE)
rule <- function(s) cat("\n", strrep("=", 70), "\n", s, "\n", strrep("=", 70), "\n", sep = "")

dat <- read_csv(file.path(clean, "lag_femur_long.csv"), show_col_types = FALSE) |>
  left_join(read_csv(file.path(clean, "species_to_tip.csv"), show_col_types = FALSE), by = "species")

ITER <- 1000; WARM <- 500; CH <- 4; ADAPT <- 0.9; TD <- 10
nlexpr <- circ_mm ~ logA - exp(-lognu) * log1p(exp(lognu) * exp(-exp(logk) * (age - t0)))

diag_block <- function(fit, label) {
  rule(paste("DIAGNOSTICS:", label)); print(summary(fit))
  np <- brms::nuts_params(fit)
  ndiv <- sum(subset(np, Parameter == "divergent__")$Value)
  td   <- subset(np, Parameter == "treedepth__")$Value
  cat(sprintf(">> divergences: %d | iters at max treedepth(%d): %d | max Rhat: %.3f\n",
              ndiv, TD, sum(td >= TD), max(brms::rhat(fit), na.rm = TRUE)))
  invisible(ndiv)
}

# ---- Stage 1: reuse saved fit -------------------------------------
fit1 <- readRDS(file.path(out, "fit_stage1_trex.rds"))

# ---- Stage 2a: all taxa, pooled, global nu ------------------------
rule("STAGE 2a (fast): pooling on logA/logk/t0 | global nu")
f2 <- bf(nlexpr,
         logA  ~ 1 + (1 | species) + (1 | specimen),
         logk  ~ 1 + (1 | species),
         lognu ~ 1,
         t0    ~ 1 + (1 | specimen),
         nl = TRUE)
p2 <- c(prior(normal(5.5, 1),  nlpar = "logA"),
        prior(normal(-1, 1),   nlpar = "logk"),
        prior(normal(0, 0.5),  nlpar = "lognu"),
        prior(normal(0, 2),    nlpar = "t0"),                       # tightened
        prior(student_t(3, 0, 0.7), class = "sd", nlpar = "logA"),
        prior(student_t(3, 0, 0.5), class = "sd", nlpar = "logk"),
        prior(student_t(3, 0, 1),   class = "sd", nlpar = "t0"))    # tightened
fit2 <- brm(f2, data = dat, family = lognormal(), prior = p2,
            chains = CH, iter = ITER, warmup = WARM, init = 0.2, seed = 1,
            control = list(adapt_delta = ADAPT, max_treedepth = TD),
            refresh = 200)
saveRDS(fit2, file.path(out, "fit_stage2_pooled_fast.rds"))
diag_block(fit2, "Stage 2a fast")

# ---- Stage 2b: nu varies by species -------------------------------
rule("STAGE 2b (fast): lognu ~ 1 + (1 | species)")
f2b <- bf(nlexpr,
          logA  ~ 1 + (1 | species) + (1 | specimen),
          logk  ~ 1 + (1 | species),
          lognu ~ 1 + (1 | species),
          t0    ~ 1 + (1 | specimen),
          nl = TRUE)
p2b <- c(p2, prior(student_t(3, 0, 0.3), class = "sd", nlpar = "lognu"))
fit2b <- brm(f2b, data = dat, family = lognormal(), prior = p2b,
             chains = CH, iter = ITER, warmup = WARM, init = 0.2, seed = 1,
             control = list(adapt_delta = ADAPT, max_treedepth = TD),
             refresh = 200)
saveRDS(fit2b, file.path(out, "fit_stage2b_nu_by_species_fast.rds"))
diag_block(fit2b, "Stage 2b fast")

# ---- Stage 3 (stretch): topology-only phylo covariance on logA ----
ok3 <- tryCatch({
  rule("STAGE 3 (fast, stretch): + (1 | gr(phylo_tip, cov = A_phylo)) on logA")
  tree <- read.tree(file.path(clean, "theropod_femur_tree.nwk"))
  A_phylo <- vcv(compute.brlen(tree, 1), corr = TRUE)
  f3 <- bf(nlexpr,
           logA  ~ 1 + (1 | species) + (1 | specimen) + (1 | gr(phylo_tip, cov = A_phylo)),
           logk  ~ 1 + (1 | species), lognu ~ 1, t0 ~ 1 + (1 | specimen), nl = TRUE)
  fit3 <- brm(f3, data = dat, family = lognormal(), prior = p2,
              data2 = list(A_phylo = A_phylo),
              chains = CH, iter = ITER, warmup = WARM, init = 0.2, seed = 1,
              control = list(adapt_delta = 0.95, max_treedepth = 11), refresh = 200)
  saveRDS(fit3, file.path(out, "fit_stage3_phylo_fast.rds"))
  diag_block(fit3, "Stage 3 fast"); TRUE
}, error = function(e) { cat("\n[Stage 3 skipped] ", conditionMessage(e), "\n"); FALSE })

# ---- Figures ------------------------------------------------------
tryCatch({
  ag <- seq(min(dat$age), max(dat$age) + 2, length.out = 80)
  # Stage 1 fit on Sue
  d1 <- filter(dat, specimen == "FMNH_PR_2081")
  e1 <- posterior_epred(fit1, newdata = data.frame(age = ag))
  g1 <- data.frame(age = ag, med = apply(e1, 2, median),
                   lo = apply(e1, 2, quantile, .05), hi = apply(e1, 2, quantile, .95))
  ggsave(file.path(out, "fig_stage1_trex_fit.png"),
    ggplot() +
      geom_ribbon(data = g1, aes(age, ymin = lo, ymax = hi), alpha = .25, fill = "steelblue") +
      geom_line(data = g1, aes(age, med), colour = "steelblue", linewidth = 1) +
      geom_point(data = d1, aes(age, circ_mm, shape = is_efs), size = 2) +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), name = "EFS") +
      labs(title = "Stage 1: Richards fit, T. rex FMNH PR 2081 (Sue)",
           x = "LAG index (age proxy)", y = "femur circumference (mm)") +
      theme_minimal(base_size = 12),
    width = 6, height = 4.2, dpi = 120)

  # Stage 2a species curves
  sp <- sort(unique(dat$species))
  nd <- expand.grid(age = ag, species = sp, specimen = "._d", phylo_tip = "._d",
                    stringsAsFactors = FALSE)
  e2 <- posterior_epred(fit2, newdata = nd, re_formula = ~ (1 | species), allow_new_levels = TRUE)
  nd$med <- apply(e2, 2, median); nd$lo <- apply(e2, 2, quantile, .05); nd$hi <- apply(e2, 2, quantile, .95)
  ggsave(file.path(out, "fig_stage2_species_curves.png"),
    ggplot() +
      geom_ribbon(data = nd, aes(age, ymin = lo, ymax = hi), alpha = .2, fill = "darkorange") +
      geom_line(data = nd, aes(age, med), colour = "darkorange") +
      geom_point(data = dat, aes(age, circ_mm), size = 1) +
      facet_wrap(~ species, scales = "free_y") +
      labs(title = "Stage 2a: partially-pooled species growth curves",
           x = "LAG index (age proxy)", y = "femur circumference (mm)") +
      theme_minimal(base_size = 10),
    width = 9, height = 6, dpi = 120)

  # nu posteriors: Stage 1 (single) vs Stage 2a global vs Stage 2b per species
  draws1  <- as_draws_df(fit1);  draws2 <- as_draws_df(fit2); draws2b <- as_draws_df(fit2b)
  df_nu <- rbind(
    data.frame(group = "T. rex only (Stage 1)",  nu = exp(draws1$b_lognu_Intercept)),
    data.frame(group = "GLOBAL pooled (Stage 2a)", nu = exp(draws2$b_lognu_Intercept)))
  for (s in sp) {
    col <- paste0("r_species__lognu[", s, ",Intercept]")
    if (col %in% names(draws2b))
      df_nu <- rbind(df_nu, data.frame(group = paste0(s, " (2b)"),
                                       nu = exp(draws2b$b_lognu_Intercept + draws2b[[col]])))
  }
  summ <- df_nu |> group_by(group) |>
    summarise(m = median(nu), lo50 = quantile(nu, .25), hi50 = quantile(nu, .75),
              lo90 = quantile(nu, .05), hi90 = quantile(nu, .95), .groups = "drop")
  ggsave(file.path(out, "fig_nu_posteriors.png"),
    ggplot(summ, aes(y = group)) +
      geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
      geom_linerange(aes(xmin = lo90, xmax = hi90)) +
      geom_linerange(aes(xmin = lo50, xmax = hi50), linewidth = 1.6, colour = "purple4") +
      geom_point(aes(x = m), size = 2.4) +
      coord_cartesian(xlim = c(0, quantile(df_nu$nu, .98))) +
      labs(title = expression(paste("Shape ", nu, ": single series vs pooled vs per-species")),
           subtitle = "median, 50% & 90% CrI; dashed = 1 (logistic)", x = expression(nu), y = NULL) +
      theme_minimal(base_size = 10),
    width = 7, height = 5, dpi = 120)
  cat("\nFigures written.\n")
}, error = function(e) cat("\n[figure warning] ", conditionMessage(e), "\n"))

rule("nu = exp(lognu) posterior summaries")
nq <- function(x) sprintf("%.2f [%.2f, %.2f]", median(x), quantile(x, .05), quantile(x, .95))
cat("Stage 1  (Sue only):        nu =", nq(exp(as_draws_df(fit1)$b_lognu_Intercept)), "\n")
cat("Stage 2a (all-clade pooled): nu =", nq(exp(as_draws_df(fit2)$b_lognu_Intercept)), "\n")
rule(paste("FAST FEASIBILITY RUN COMPLETE | Stage 3 ran:", ok3))
