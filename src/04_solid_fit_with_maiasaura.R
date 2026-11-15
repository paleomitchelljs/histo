#!/usr/bin/env Rscript
# =====================================================================
# 04_solid_fit_with_maiasaura.R
# The fit the data extraction was for: add the Maiasaura full-curve series
# (perinate->adult tibia circumference) to the resorption-truncated theropod
# femora and ask whether the all-clade-pooled shape ν tightens.
#
# Data: data_clean/lag_circumference_long.csv (102 theropod femur + 51 Maiasaura
#       tibia; 2 clades, 2 elements). Age = years (LAG index for theropods, true
#       age incl. age-0 perinate for Maiasaura); the per-specimen t0 absorbs the
#       theropods' unknown resorption offset, while Maiasaura pins the lower curve.
#
# Compare ν here vs the theropod-only fit (outputs/fit_stage2_pooled.rds).
# rstan, treedepth=13 (slow but converges). Saves incrementally.
# =====================================================================
suppressMessages({ library(brms); library(dplyr); library(readr); library(ggplot2) })
options(mc.cores = 1, brms.backend = "rstan"); set.seed(1)
root <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
clean <- file.path(root, "data_clean"); out <- file.path(root, "outputs")

con <- file(file.path(out, "solid_maia_report.txt"), open = "wt")
sink(con, split = TRUE); on.exit({ sink(); close(con) }, add = TRUE)
rule <- function(s) cat("\n", strrep("=", 70), "\n", s, "\n", strrep("=", 70), "\n", sep = "")

dat <- read_csv(file.path(clean, "lag_circumference_long.csv"), show_col_types = FALSE)
cat(sprintf("data: %d rows | %d species | %d specimens | elements: %s\n",
            nrow(dat), n_distinct(dat$species), n_distinct(dat$specimen),
            paste(unique(dat$element), collapse = "/")))

nlexpr <- circ_mm ~ logA - exp(-lognu) * log1p(exp(lognu) * exp(-exp(logk) * (age - t0)))
f <- bf(nlexpr,
        logA  ~ 1 + (1 | species) + (1 | specimen),
        logk  ~ 1 + (1 | species),
        lognu ~ 1,                                  # global shape, now Maiasaura-informed
        t0    ~ 1 + (1 | specimen),
        nl = TRUE)
p <- c(prior(normal(5.5, 1),  nlpar = "logA"),
       prior(normal(-1, 1),   nlpar = "logk"),
       prior(normal(0, 0.5),  nlpar = "lognu"),
       prior(normal(0, 2),    nlpar = "t0"),
       prior(student_t(3, 0, 0.7), class = "sd", nlpar = "logA"),
       prior(student_t(3, 0, 0.5), class = "sd", nlpar = "logk"),
       prior(student_t(3, 0, 1),   class = "sd", nlpar = "t0"))

rule("SOLID FIT: theropod femora + Maiasaura tibiae | global nu (treedepth=13)")
fit <- brm(f, data = dat, family = lognormal(), prior = p,
           chains = 4, iter = 1500, warmup = 750, init = 0.2, seed = 1,
           control = list(adapt_delta = 0.97, max_treedepth = 13), refresh = 250)
saveRDS(fit, file.path(out, "fit_solid_maia.rds"))

print(summary(fit))
np <- nuts_params(fit)
cat(sprintf(">> divergences: %d | max Rhat: %.3f | min bulk-ESS: %.0f\n",
            sum(subset(np, Parameter == "divergent__")$Value),
            max(rhat(fit), na.rm = TRUE),
            min(posterior::summarise_draws(as_draws(fit), "ess_bulk")$ess_bulk, na.rm = TRUE)))

rule("DOES MAIASAURA SHARPEN nu?  (global nu = exp(lognu))")
nq <- function(x) sprintf("%.2f [%.2f, %.2f]", median(x), quantile(x, .05), quantile(x, .95))
nu_now <- exp(as_draws_df(fit)$b_lognu_Intercept)
cat("theropods + Maiasaura: nu =", nq(nu_now), "\n")
to <- tryCatch(exp(as_draws_df(readRDS(file.path(out, "fit_stage2_pooled.rds")))$b_lognu_Intercept),
               error = function(e) NULL)
if (!is.null(to)) cat("theropods only (2a):   nu =", nq(to),
    sprintf("  [90%% width %.2f -> %.2f]\n", diff(quantile(to, c(.05, .95))),
            diff(quantile(nu_now, c(.05, .95)))))

# species growth curves (incl. Maiasaura)
tryCatch({
  sp <- sort(unique(dat$species))
  gg <- do.call(rbind, lapply(sp, function(s){
    a <- seq(0, max(dat$age[dat$species == s]) + 1, length.out = 60)
    nd <- data.frame(age = a, species = s, specimen = "._d")
    e <- posterior_epred(fit, newdata = nd, re_formula = ~ (1 | species), allow_new_levels = TRUE)
    data.frame(species = s, age = a, med = apply(e, 2, median),
               lo = apply(e, 2, quantile, .05), hi = apply(e, 2, quantile, .95))
  }))
  ggsave(file.path(out, "fig_solid_maia_curves.png"),
    ggplot() + geom_ribbon(data = gg, aes(age, ymin = lo, ymax = hi), alpha = .2, fill = "seagreen") +
      geom_line(data = gg, aes(age, med), colour = "seagreen") +
      geom_point(data = dat, aes(age, circ_mm), size = .8) +
      facet_wrap(~ species, scales = "free") +
      labs(title = "Solid fit + Maiasaura: pooled growth curves",
           x = "age (yr)", y = "circumference (mm)") + theme_minimal(base_size = 9),
    width = 10, height = 6.5, dpi = 120)
  cat("\nfig_solid_maia_curves.png written\n")
}, error = function(e) cat("[fig skipped] ", conditionMessage(e), "\n"))
rule("DONE")
