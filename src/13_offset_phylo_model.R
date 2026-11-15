#!/usr/bin/env Rscript
# Fit the offset + categorical model WITH phylogenetic covariance on nu
# (13_offset_phylo_model.stan) on the expanded cross-Dinosauria dataset.
suppressMessages({ library(cmdstanr); library(dplyr); library(readr); library(ape); library(posterior) })
root  <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
clean <- file.path(root,"data_clean"); out <- file.path(root,"outputs")

dat <- read_csv(file.path(clean,"lag_size_long.csv"), show_col_types=FALSE) %>% filter(!is.na(size), !is.na(age))
dat$species <- factor(dat$species); dat$specimen <- factor(dat$specimen); dat$dimr <- as.integer(factor(dat$dimension))
ind_tbl <- dat %>% distinct(specimen, species, resorption, efs_present) %>%
  left_join(dat %>% group_by(specimen) %>% summarise(logmax=log(max(size)), .groups="drop"), by="specimen") %>%
  arrange(specimen)
ind_tbl$sp_of_ind <- as.integer(ind_tbl$species)

# phylogenetic correlation (topology-only: Grafen branch lengths -> ultrametric) + Cholesky
tree  <- read.tree(file.path(clean,"dino_tree.nwk"))
treeU <- compute.brlen(tree, method="Grafen")
R     <- vcv(treeU, corr=TRUE)                       # T x T correlation, tips in tree order
tips  <- rownames(R)
L     <- t(chol(R))                                  # lower Cholesky: L %*% t(L) == R
tipmap <- read_csv(file.path(clean,"species_to_tip.csv"), show_col_types=FALSE)
sp_levels <- levels(dat$species)
tip_of_sp <- match(tipmap$phylo_tip[match(sp_levels, tipmap$species)], tips)
stopifnot(!any(is.na(tip_of_sp)))

standata <- list(N=nrow(dat), S=nlevels(dat$species), I=nlevels(dat$specimen), D=max(dat$dimr), T=length(tips),
  size=dat$size, age0=dat$age, ind=as.integer(dat$specimen), dimr=dat$dimr,
  sp_of_ind=ind_tbl$sp_of_ind, tip_of_sp=tip_of_sp, Lphylo=L,
  efs=as.integer(ind_tbl$efs_present), resorb=as.integer(ind_tbl$resorption!="absent"), logmax=ind_tbl$logmax)

mod <- cmdstan_model(file.path(root,"src/13_offset_phylo_model.stan"))
fit <- mod$sample(data=standata, chains=4, parallel_chains=4, iter_warmup=1000, iter_sampling=1000,
                  adapt_delta=0.97, max_treedepth=13, refresh=250, seed=1, init=0.5)
fit$save_object(file.path(out,"fit_offset_phylo.rds"))
saveRDS(list(species=sp_levels, specimen=levels(dat$specimen), ind_tbl=ind_tbl,
             tips=tips, tip_of_sp=tip_of_sp), file.path(out,"offset_phylo_index.rds"))

cat("\n=== diagnostics ===\n"); print(fit$diagnostic_summary())
cat("\n=== phylogenetic signal in nu (lambda_nu) + key params ===\n")
print(as.data.frame(summarise_draws(fit$draws(c("lambda_nu","sigma_phy","sd_lognu","mu_lognu","tau_off")))), digits=3)
cat("\nDONE phylo model\n")
