#!/usr/bin/env Rscript
# =====================================================================
# 05_unified_nu_by_species.R
# 04 showed a single GLOBAL nu is misspecified: forced to Maiasaura's
# Gompertz-ward shape (nu=0.18) it underfits the large theropods (T. rex,
# Meraxes plateau ~500-570 mm but the curve caps ~310). So let shape vary:
# lognu ~ 1 + (1 | species), partially pooled. With Maiasaura supplying a
# full curve, every species' nu should now be identifiable (each borrows the
# lower-curve information through the shared hierarchy).
#
# Unified data (theropod femora + Maiasaura tibiae). rstan, treedepth=13.
# Compares to: 04 global-nu (fit_solid_maia.rds) and theropod-only 2b.
# =====================================================================
suppressMessages({ library(brms); library(dplyr); library(readr); library(ggplot2) })
options(mc.cores = 1, brms.backend = "rstan"); set.seed(1)
root <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
clean <- file.path(root, "data_clean"); out <- file.path(root, "outputs")

con <- file(file.path(out, "unified_nu_report.txt"), open = "wt")
sink(con, split = TRUE); on.exit({ sink(); close(con) }, add = TRUE)
rule <- function(s) cat("\n", strrep("=", 70), "\n", s, "\n", strrep("=", 70), "\n", sep = "")

dat <- read_csv(file.path(clean, "lag_circumference_long.csv"), show_col_types = FALSE)
nlexpr <- circ_mm ~ logA - exp(-lognu) * log1p(exp(lognu) * exp(-exp(logk) * (age - t0)))
f <- bf(nlexpr,
        logA  ~ 1 + (1 | species) + (1 | specimen),
        logk  ~ 1 + (1 | species),
        lognu ~ 1 + (1 | species),                  # shape varies, partially pooled
        t0    ~ 1 + (1 | specimen),
        nl = TRUE)
p <- c(prior(normal(5.5, 1),  nlpar = "logA"),
       prior(normal(-1, 1),   nlpar = "logk"),
       prior(normal(0, 0.5),  nlpar = "lognu"),
       prior(normal(0, 2),    nlpar = "t0"),
       prior(student_t(3, 0, 0.7), class = "sd", nlpar = "logA"),
       prior(student_t(3, 0, 0.5), class = "sd", nlpar = "logk"),
       prior(student_t(3, 0, 0.5), class = "sd", nlpar = "lognu"),
       prior(student_t(3, 0, 1),   class = "sd", nlpar = "t0"))

rule("UNIFIED FIT: nu ~ 1 + (1|species)  (treedepth=13)")
fit <- brm(f, data = dat, family = lognormal(), prior = p,
           chains = 4, iter = 1500, warmup = 750, init = 0.2, seed = 1,
           control = list(adapt_delta = 0.97, max_treedepth = 13), refresh = 250)
saveRDS(fit, file.path(out, "fit_unified_nu_by_species.rds"))
print(summary(fit))
np <- nuts_params(fit)
cat(sprintf(">> divergences: %d | max Rhat: %.3f | min bulk-ESS: %.0f\n",
            sum(subset(np, Parameter == "divergent__")$Value),
            max(rhat(fit), na.rm = TRUE),
            min(posterior::summarise_draws(as_draws(fit), "ess_bulk")$ess_bulk, na.rm = TRUE)))

rule("PER-SPECIES nu = exp(lognu_Intercept + r_species)")
d <- as_draws_df(fit)
for (s in sort(unique(dat$species))) {
  c0 <- paste0("r_species__lognu[", s, ",Intercept]")
  if (c0 %in% names(d)) { v <- exp(d$b_lognu_Intercept + d[[c0]])
    cat(sprintf("  %-30s %.2f [%.2f, %.2f]\n", s, median(v), quantile(v, .05), quantile(v, .95))) }
}
cat("\n(nu<1 Gompertz-ward, =1 logistic, >1 super-logistic)\n")

# LOO vs global-nu (04)
tryCatch({
  rule("LOO: nu-by-species vs global nu (both on unified data)")
  fg <- readRDS(file.path(out, "fit_solid_maia.rds"))
  print(loo_compare(loo(fit), loo(fg)))
  cat("(row 1 = better; 'fit' = nu-by-species, 'fg' = global nu)\n")
}, error = function(e) cat("[loo skipped] ", conditionMessage(e), "\n"))

# curves
tryCatch({
  sp <- sort(unique(dat$species))
  gg <- do.call(rbind, lapply(sp, function(s){
    a <- seq(0, max(dat$age[dat$species == s]) + 1, length.out = 60)
    e <- posterior_epred(fit, newdata = data.frame(age = a, species = s, specimen = "._d"),
                         re_formula = ~ (1 | species), allow_new_levels = TRUE)
    data.frame(species = s, age = a, med = apply(e, 2, median),
               lo = apply(e, 2, quantile, .05), hi = apply(e, 2, quantile, .95))
  }))
  ggsave(file.path(out, "fig_unified_nu_curves.png"),
    ggplot() + geom_ribbon(data = gg, aes(age, ymin = lo, ymax = hi), alpha = .2, fill = "firebrick") +
      geom_line(data = gg, aes(age, med), colour = "firebrick") +
      geom_point(data = dat, aes(age, circ_mm), size = .8) +
      facet_wrap(~ species, scales = "free") +
      labs(title = "Unified fit, nu by species: pooled growth curves",
           x = "age (yr)", y = "circumference (mm)") + theme_minimal(base_size = 9),
    width = 10, height = 6.5, dpi = 120)
  cat("\nfig_unified_nu_curves.png written\n")
}, error = function(e) cat("[fig skipped] ", conditionMessage(e), "\n"))
rule("DONE")
