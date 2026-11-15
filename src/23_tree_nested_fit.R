#!/usr/bin/env Rscript
# Fit the tree-nested multivariate model (23) to the empirical LAG data and write
# the variance / ancestry-vs-independent partitions and the per-taxon precision
# table (where more data would help most).
# Usage: Rscript 23_tree_nested_fit.R [iter=1000] [chains=4] [tag]
#   tag e.g. "next" reads data_clean/*_next.* (staged data) and writes outputs/*_next.*
suppressMessages({ library(cmdstanr); library(dplyr); library(readr); library(ape); library(posterior) })
a <- commandArgs(trailingOnly=TRUE)
iter <- ifelse(is.na(a[1]), 1000L, as.integer(a[1])); chains <- ifelse(is.na(a[2]), 4L, as.integer(a[2]))
tag <- ifelse(is.na(a[3]), "", paste0("_", a[3]))
root  <- "/Users/jmitchell/Library/CloudStorage/Dropbox/Research/00_current_Histology"
clean <- file.path(root,"data_clean"); out <- file.path(root,"outputs")
source(file.path(root,"src/23_tree_helpers.R"))

dat  <- read_csv(file.path(clean,paste0("lag_size_long",tag,".csv")), show_col_types=FALSE) %>% filter(!is.na(size), !is.na(age))
tree <- read.tree(file.path(clean,paste0("dino_tree",tag,".nwk")))
tmap <- read_csv(file.path(clean,paste0("species_to_tip",tag,".csv")), show_col_types=FALSE)
sd0  <- make_standata(dat, tree, setNames(tmap$phylo_tip, tmap$species))

mod <- cmdstan_model(file.path(root,"src/23_tree_nested_model.stan"))
fit <- mod$sample(data=sd0$standata, chains=chains, parallel_chains=chains, iter_warmup=iter, iter_sampling=iter,
                  adapt_delta=0.95, max_treedepth=12, refresh=200, seed=1, init=0.3)
fit$save_object(file.path(out,paste0("fit_tree_nested",tag,".rds")))
saveRDS(sd0, file.path(out,paste0("tree_nested_index",tag,".rds")))

cat("\n=== diagnostics ===\n"); print(fit$diagnostic_summary())
cat("\n=== hyperparameters ===\n")
print(as.data.frame(summarise_draws(fit$draws(c("lambda","sd_phy","sd_res","sd_logA_ind","tau_off","sigma")),
      "median", ~quantile(.x, c(.05,.95)), "rhat", "ess_bulk")), digits=3)

S <- sd0$standata$S; B <- sd0$standata$B; P <- sd0$tm$P; tos <- sd0$tip_of_sp
delta <- get_arr(fit, "delta", B); dev <- get_arr(fit, "dev", S)
pp <- pair_partition(delta, dev, P, tos); ps <- partition_summary(pp, delta, dev, P, tos)
q <- function(m) t(apply(m, 2, quantile, c(.05,.5,.95)))
cat("\n=== realized tree share of among-species variance (5/50/95%) ===\n"); print(round(q(ps$phy_real),2))
cat("\n=== ancestry share of pairwise similarity (5/50/95%) ===\n");        print(round(q(ps$anc_share),2))

# pairwise table (posterior means + P(independent > 0))
pt <- do.call(rbind, lapply(seq_along(TRAITS), function(k) data.frame(trait=TRAITS[k],
  sp_i=sd0$species[pp$pairs[,1]], sp_j=sd0$species[pp$pairs[,2]],
  inherited=colMeans(pp$inherited[,,k]), independent=colMeans(pp$independent[,,k]),
  cross=colMeans(pp$cross[,,k]), p_parallel=colMeans(pp$independent[,,k] > 0))))
write_csv(pt, file.path(out,paste0("tree_nested_pairwise",tag,".csv")))

# per-species precision: posterior SD of each trait, and how much of the tip
# estimate is prior/tree-borrowed (posterior SD relative to the between-species SD)
spv <- do.call(rbind, lapply(seq_along(TRAITS), function(k){
  tot <- sqrt(median(apply(dev[,,k], 1, var)))
  data.frame(trait=TRAITS[k], species=sd0$species, post_sd=apply(dev[,,k], 2, sd),
             rel_uncert=apply(dev[,,k], 2, sd)/tot) }))
cnt <- dat %>% group_by(species) %>% summarise(n_spec=n_distinct(specimen), n_rings=n(), .groups="drop")
spv <- left_join(spv, cnt, by="species")
write_csv(spv, file.path(out,paste0("tree_nested_species_precision",tag,".csv")))
cat("\n=== least-constrained species per trait (rel_uncert = posterior SD / among-species SD) ===\n")
print(as.data.frame(spv %>% group_by(trait) %>% slice_max(rel_uncert, n=6) %>% mutate(across(where(is.numeric), ~round(.x,2)))))
cat("\nDONE tree-nested fit\n")
