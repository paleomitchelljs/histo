#!/usr/bin/env Rscript
# Fit the offset + categorical-histology Richards model (08_offset_model.stan)
# on the unified cross-Dinosauria dataset. cmdstanr backend (now available).
suppressMessages({ library(cmdstanr); library(dplyr); library(readr); library(posterior) })
root  <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
clean <- file.path(root, "data_clean"); out <- file.path(root, "outputs")

dat <- read_csv(file.path(clean, "lag_size_long.csv"), show_col_types = FALSE) %>%
  filter(!is.na(size), !is.na(age))
dat$species  <- factor(dat$species)
dat$specimen <- factor(dat$specimen)
dat$dimr     <- as.integer(factor(dat$dimension))

ind_tbl <- dat %>% distinct(specimen, species, resorption, efs_present) %>%
  left_join(dat %>% group_by(specimen) %>% summarise(logmax = log(max(size)), .groups = "drop"),
            by = "specimen") %>%
  arrange(specimen)
ind_tbl$sp_of_ind <- as.integer(ind_tbl$species)

standata <- list(
  N = nrow(dat), S = nlevels(dat$species), I = nlevels(dat$specimen), D = max(dat$dimr),
  size = dat$size, age0 = dat$age, ind = as.integer(dat$specimen), dimr = dat$dimr,
  sp_of_ind = ind_tbl$sp_of_ind,
  efs    = as.integer(ind_tbl$efs_present),
  resorb = as.integer(ind_tbl$resorption != "absent"),
  logmax = ind_tbl$logmax)

mod <- cmdstan_model(file.path(root, "src/08_offset_model.stan"))
fit <- mod$sample(data = standata, chains = 4, parallel_chains = 4,
                  iter_warmup = 1000, iter_sampling = 1000, adapt_delta = 0.95,
                  max_treedepth = 12, refresh = 250, seed = 1, init = 0.5)

fit$save_object(file.path(out, "fit_offset_model.rds"))
saveRDS(list(species = levels(dat$species), specimen = levels(dat$specimen),
             ind_tbl = ind_tbl), file.path(out, "offset_model_index.rds"))

cat("\n=== diagnostics ===\n"); print(fit$diagnostic_summary())
draws <- fit$draws(format = "df")
sm <- summarise_draws(fit$draws(c("mu_lognu","mu_logk","mu_t0","tau_off","sigma")))
print(as.data.frame(sm), digits = 3)

# per-species nu
cat("\n=== per-species nu = exp(lognu_sp) ===\n")
sp <- levels(dat$species)
nu_tab <- sapply(seq_along(sp), function(j){
  v <- exp(draws[[paste0("lognu_sp[", j, "]")]]); c(median(v), quantile(v, c(.05,.95)))
})
nu_df <- data.frame(species = sp, nu = round(nu_tab[1,],2),
                    lo = round(nu_tab[2,],2), hi = round(nu_tab[3,],2))
print(nu_df, row.names = FALSE)

# estimated missing rings per resorbed specimen
cat("\n=== estimated missing inner rings (offset) for resorbed specimens ===\n")
off <- sapply(seq_len(nrow(ind_tbl)), function(i) median(draws[[paste0("miss_yr[", i, "]")]]))
od <- data.frame(specimen = ind_tbl$specimen, resorption = ind_tbl$resorption,
                 efs = ind_tbl$efs_present, miss_yr = round(off,1))
print(od[od$resorption != "absent", ], row.names = FALSE)
cat("\nDONE offset model\n")
