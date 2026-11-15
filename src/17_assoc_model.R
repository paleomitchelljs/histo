#!/usr/bin/env Rscript
# Run the categorical-character -> shape (nu) association model (17_assoc_model.stan)
# with vascular complexity (simple/complex/unknown) coded from Cullen S1 + Bybee +
# D'Emic. Reference level = "unknown" (the uncoded taxa, which still inform the
# phylo/nu structure); the test of interest is the complex - simple contrast among
# the 4 coded theropods, phylogeny-controlled and nu-uncertainty propagated.
suppressMessages({ library(cmdstanr); library(dplyr); library(readr); library(ape); library(posterior) })
root  <- "/Users/jmitchell/Library/CloudStorage/Dropbox/Research/00_current_Histology"
clean <- file.path(root,"data_clean"); out <- file.path(root,"outputs")

dat <- read_csv(file.path(clean,"lag_size_long.csv"), show_col_types=FALSE) %>% filter(!is.na(size), !is.na(age))
dat$species <- factor(dat$species); dat$specimen <- factor(dat$specimen); dat$dimr <- as.integer(factor(dat$dimension))
it <- dat %>% distinct(specimen, species, resorption, efs_present) %>%
  left_join(dat %>% group_by(specimen) %>% summarise(logmax=log(max(size)), .groups="drop"), by="specimen") %>%
  arrange(specimen); it$sp_of_ind <- as.integer(it$species)
tree <- read.tree(file.path(clean,"dino_tree.nwk")); R <- vcv(compute.brlen(tree,method="Grafen"),corr=TRUE)
tips <- rownames(R); L <- t(chol(R)); tm <- read_csv(file.path(clean,"species_to_tip.csv"),show_col_types=FALSE)
spl <- levels(dat$species); tos <- match(tm$phylo_tip[match(spl,tm$species)], tips)

# Xnu from vasc_complexity (ref = unknown): col1 = simple, col2 = complex
feat <- read_csv(file.path(clean,"histo_features.csv"), show_col_types=FALSE)
vc <- feat$vasc_complexity[match(spl, feat$species)]; vc[is.na(vc)] <- "unknown"
Xnu <- cbind(simple = as.numeric(vc=="simple"), complex = as.numeric(vc=="complex"))

standata <- list(N=nrow(dat), S=nlevels(dat$species), I=nlevels(dat$specimen), D=max(dat$dimr), T=length(tips), K=2L,
  size=dat$size, age0=dat$age, ind=as.integer(dat$specimen), dimr=dat$dimr, sp_of_ind=it$sp_of_ind,
  tip_of_sp=tos, Lphylo=L, Xnu=Xnu, efs=as.integer(it$efs_present),
  resorb=as.integer(it$resorption!="absent"), logmax=it$logmax)
cat("coded taxa: simple =", sum(Xnu[,"simple"]), "| complex =", sum(Xnu[,"complex"]),
    "| unknown(ref) =", sum(rowSums(Xnu)==0), "\n")

mod <- cmdstan_model(file.path(root,"src/17_assoc_model.stan"))
fit <- mod$sample(data=standata, chains=3, parallel_chains=1, iter_warmup=800, iter_sampling=800,
                  adapt_delta=0.95, max_treedepth=11, refresh=500, seed=1, init=0.4)
fit$save_object(file.path(out,"fit_assoc_vascular.rds"))

cat("\n=== diagnostics ===\n"); print(fit$diagnostic_summary())
d <- fit$draws(format="df")
g1 <- d[["gamma_out[1]"]]; g2 <- d[["gamma_out[2]"]]; contr <- g2 - g1   # complex - simple
q <- function(x) sprintf("%.2f [%.2f, %.2f]", median(x), quantile(x,.05), quantile(x,.95))
cat("\n=== vascular -> shape (nu) association (log-nu scale) ===\n")
cat("gamma simple  vs unknown :", q(g1), "\n")
cat("gamma complex vs unknown :", q(g2), "\n")
cat("CONTRAST complex - simple:", q(contr), "  (excludes 0:", (quantile(contr,.05)>0)|(quantile(contr,.95)<0), ")\n")
cat("lambda_nu:", q(d[["lambda_nu"]]), "\n\nDONE assoc model\n")
