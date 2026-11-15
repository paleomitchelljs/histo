#!/usr/bin/env Rscript
# Fit the LAG-width growth-rate model (15_lagwidth_model.stan) on the
# circumference taxa (radial mm/yr comparable), active-growth rings only.
suppressMessages({ library(cmdstanr); library(dplyr); library(readr); library(ape); library(posterior); library(ggplot2) })
root  <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
clean <- file.path(root,"data_clean"); out <- file.path(root,"outputs")

w <- read_csv(file.path(clean,"lag_widths_long.csv"), show_col_types=FALSE) %>%
  filter(dimension=="circumference_mm", !is_efs, lag_width>0)   # radial mm/yr, active growth
w$species  <- factor(w$species); w$specimen <- factor(w$specimen)
w$logw <- log(w$lag_width)

# prune tree to the circumference-taxa tips
tree   <- read.tree(file.path(clean,"dino_tree.nwk"))
tipmap <- read_csv(file.path(clean,"species_to_tip.csv"), show_col_types=FALSE)
sp_levels <- levels(w$species)
tips_needed <- unique(tipmap$phylo_tip[match(sp_levels, tipmap$species)])
tree <- keep.tip(tree, tips_needed)
R <- vcv(compute.brlen(tree, method="Grafen"), corr=TRUE); tips <- rownames(R); L <- t(chol(R))
tip_of_sp <- match(tipmap$phylo_tip[match(sp_levels, tipmap$species)], tips)
stopifnot(!any(is.na(tip_of_sp)))

standata <- list(N=nrow(w), S=nlevels(w$species), I=nlevels(w$specimen), T=length(tips),
  logw=w$logw, age=w$age, sp=as.integer(w$species), id=as.integer(w$specimen),
  tip_of_sp=tip_of_sp, Lphylo=L, ref_age=2)

cat("LAG-width rate model:", nrow(w),"active-growth widths |", nlevels(w$species),"species |",
    nlevels(w$specimen),"specimens |", length(tips),"tips\n")
mod <- cmdstan_model(file.path(root,"src/15_lagwidth_model.stan"))
fit <- mod$sample(data=standata, chains=4, parallel_chains=2, iter_warmup=1000, iter_sampling=1000,
                  adapt_delta=0.95, max_treedepth=11, refresh=500, seed=1)
fit$save_object(file.path(out,"fit_lagwidth.rds"))
saveRDS(list(species=sp_levels, tips=tips, tip_of_sp=tip_of_sp, tree=tree), file.path(out,"lagwidth_index.rds"))

cat("\n=== diagnostics ===\n"); print(fit$diagnostic_summary())
d <- fit$draws(format="df")
cat("\n=== phylo signal in growth-rate level (lambda_rate) + decline (beta) ===\n")
print(as.data.frame(summarise_draws(fit$draws(c("lambda_rate","beta","sigma_phy","sd_sp","sigma")))), digits=3)
cat("\n=== radial apposition rate at age 2 (mm/yr) per taxon ===\n")
rate <- do.call(rbind, lapply(seq_along(sp_levels), function(s){
  v <- exp(d[[paste0("lograte_ref[",s,"]")]])
  data.frame(species=sp_levels[s], rate_mm_yr=round(median(v),2), lo=round(quantile(v,.05),2), hi=round(quantile(v,.95),2))}))
print(rate[order(-rate$rate_mm_yr),], row.names=FALSE)
cat("\nDONE lag-width model\n")
