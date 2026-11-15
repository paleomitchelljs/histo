#!/usr/bin/env Rscript
# Fit the offset + phylogenetic-nu model (13) to a simulated scenario and score
# recovery against known truth: per-species shape nu, the per-specimen missing-LAG
# offset (resorption before first LAG), and NEONATE mass (size at true age 0 --
# never asserted, must be back-inferred via the offset), plus whether lambda_nu
# correctly reads the regime (iid -> low, BM -> high).
# Usage: Rscript 19_sim_fit_assess.R <iid|bm|shift>
suppressMessages({ library(cmdstanr); library(dplyr); library(readr); library(ape); library(posterior); library(ggplot2) })
scn <- commandArgs(trailingOnly=TRUE)[1]; if (is.na(scn)) scn <- "iid"
root  <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
clean <- file.path(root,"data_clean"); sim <- file.path(clean,"sim"); out <- file.path(root,"outputs")

dat <- read_csv(file.path(sim, paste0("scenario_",scn,"_data.csv")), show_col_types=FALSE)
tsp <- read_csv(file.path(sim, paste0("scenario_",scn,"_truth_sp.csv")), show_col_types=FALSE)
tind<- read_csv(file.path(sim, paste0("scenario_",scn,"_truth_ind.csv")), show_col_types=FALSE)
dat$species <- factor(dat$species); dat$specimen <- factor(dat$specimen); dat$dimr <- 1L
it <- dat %>% distinct(specimen, species, resorption, efs_present) %>%
  left_join(dat %>% group_by(specimen) %>% summarise(logmax=log(max(size)), .groups="drop"), by="specimen") %>%
  arrange(specimen); it$sp_of_ind <- as.integer(it$species)
spl <- levels(dat$species)
tree <- compute.brlen(keep.tip(read.tree(file.path(clean,"dino_tree.nwk")), spl), method="Grafen")
R <- vcv(tree, corr=TRUE); tips <- rownames(R); L <- t(chol(R)); tip_of_sp <- match(spl, tips)

standata <- list(N=nrow(dat), S=nlevels(dat$species), I=nlevels(dat$specimen), D=1L, T=length(tips),
  size=dat$size, age0=dat$age, ind=as.integer(dat$specimen), dimr=dat$dimr, sp_of_ind=it$sp_of_ind,
  tip_of_sp=tip_of_sp, Lphylo=L, efs=as.integer(it$efs_present),
  resorb=as.integer(it$resorption!="absent"), logmax=it$logmax)

mod <- cmdstan_model(file.path(root,"src/13_offset_phylo_model.stan"))
fit <- mod$sample(data=standata, chains=2, parallel_chains=1, iter_warmup=500, iter_sampling=500,
                  adapt_delta=0.95, max_treedepth=11, refresh=0, seed=1, init=0.4)
fit$save_object(file.path(out, paste0("fit_sim_",scn,".rds")))
d <- fit$draws(format="df")
gv <- function(p,i) d[[paste0(p,"[",i,"]")]]

# --- per-species nu recovery ---
nu_rec <- sapply(seq_along(spl), function(j) median(exp(gv("lognu_sp",j))))
nu_df <- data.frame(species=spl, nu_true=tsp$nu[match(spl,tsp$species)], nu_rec=nu_rec)
# --- per-specimen offset + neonate recovery ---
sp_idx <- it$sp_of_ind
neo_rec <- offset_rec <- numeric(nrow(it))
for (i in seq_len(nrow(it))) {
  s <- sp_idx[i]
  logA <- gv("logA_i",i); nu <- exp(gv("lognu_sp",s)); k <- exp(gv("logk_sp",s)); t0 <- gv("t0_sp",s)
  offset_rec[i] <- median(gv("miss_yr",i))
  neo_rec[i] <- median(exp(logA - (1/nu)*log1p(nu*exp(-k*(0 - t0)))))   # M(true age 0)
}
ind_df <- data.frame(specimen=as.character(it$specimen)) %>%
  left_join(tind %>% transmute(specimen, offset_true, neonate_true, resorption=NA), by="specimen")
ind_df$offset_rec <- offset_rec; ind_df$neonate_rec <- neo_rec
lam <- median(d[["lambda_nu"]])

# --- scoring ---
rmse <- function(a,b) sqrt(mean((a-b)^2))
cat(sprintf("\n===== SCENARIO %s =====\n", toupper(scn)))
cat(sprintf("nu:      cor(true,rec)=%.2f  rmse(log nu)=%.2f\n", cor(log(nu_df$nu_true),log(nu_df$nu_rec)),
            rmse(log(nu_df$nu_true),log(nu_df$nu_rec))))
cat(sprintf("offset:  cor=%.2f  rmse=%.2f yr  (true mean %.1f, rec mean %.1f)\n",
            cor(ind_df$offset_true,ind_df$offset_rec), rmse(ind_df$offset_true,ind_df$offset_rec),
            mean(ind_df$offset_true), mean(ind_df$offset_rec)))
cat(sprintf("neonate: cor(log)=%.2f  rmse(log)=%.2f\n",
            cor(log(ind_df$neonate_true),log(ind_df$neonate_rec)), rmse(log(ind_df$neonate_true),log(ind_df$neonate_rec))))
cat(sprintf("lambda_nu (recovered) = %.2f   [true-trait lambda printed by simulator]\n", lam))
div <- sum(fit$diagnostic_summary()$num_divergent); cat(sprintf("divergences=%d\n", div))

write_csv(nu_df,  file.path(out, paste0("sim_recovery_nu_",scn,".csv")))
write_csv(ind_df, file.path(out, paste0("sim_recovery_ind_",scn,".csv")))
cat("DONE sim fit", scn, "\n")
