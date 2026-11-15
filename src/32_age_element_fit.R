#!/usr/bin/env Rscript
# Fit model 32 (individual age + element level + count-and-size + three-state EFS)
# to the broad tables from src/31_stage_broad.R.
# Usage: Rscript 32_age_element_fit.R [iter=1000] [chains=4] [tag=broad] [adapt_delta=0.95] [max_treedepth=12] [data prefix=broad]
suppressMessages({ library(cmdstanr); library(dplyr); library(readr); library(ape); library(posterior) })
a <- commandArgs(trailingOnly=TRUE)
iter <- ifelse(is.na(a[1]), 1000L, as.integer(a[1])); chains <- ifelse(is.na(a[2]), 4L, as.integer(a[2]))
tag <- ifelse(is.na(a[3]), "broad", a[3]); ad <- ifelse(is.na(a[4]), 0.95, as.numeric(a[4])); mtd <- ifelse(is.na(a[5]), 12L, as.integer(a[5])); pre <- ifelse(is.na(a[6]), "broad", a[6])
root  <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
clean <- file.path(root,"data_clean"); out <- file.path(root,"outputs")
source(file.path(root,"src/23_tree_helpers.R")); source(file.path(root,"src/32_helpers.R")); TRAITS <- c("logA","logr","nu","t50")
rd <- function(f) as.data.frame(read_csv(file.path(clean,paste0(pre,"_",f)), show_col_types=FALSE))
tmap <- rd("species_to_tip.csv")
x <- make_standata32(rd("rings.csv"), rd("elements.csv"), rd("individuals.csv"),
                     read.tree(file.path(clean,paste0(pre,"_tree.nwk"))), setNames(tmap$phylo_tip, tmap$species))
d <- x$standata
cat(sprintf("S %d | tips %d | I %d (%d multi-element, %d with age estimate) | ring elements %d (%d eroded) | marks %d | count-and-size %d (%d sized, %d counted) | dims %d\n",
            d$S, d$T, d$I, d$IM, sum(d$age_type), d$ER, d$NERO, d$N, d$EC, sum(d$c_has_size), sum(d$c_has_count), d$D))
mod <- cmdstan_model(file.path(root,"src/32_age_element_model.stan"))
set.seed(1)
fit <- mod$sample(data=d, chains=chains, parallel_chains=chains, iter_warmup=iter, iter_sampling=iter,
                  adapt_delta=ad, max_treedepth=mtd, refresh=100, seed=1, init=make_init32(x))
fit$save_object(file.path(out,paste0("fit_age_element_",tag,".rds"))); saveRDS(x, file.path(out,paste0("age_element_index_",tag,".rds")))

options(width=200)
cat("\n=== diagnostics ===\n"); print(fit$diagnostic_summary())
sm <- function(v) as.data.frame(summarise_draws(fit$draws(v), "median", ~quantile(.x, c(.05,.95)), "rhat", "ess_bulk"))
cat("\n=== hyperparameters ===\n")
h <- sm(c("lambda","sd_phy","sd_res","mu_logr","mu_nu","mu_t50","sd_logA_ind","rho_ind","tau_off","sigma_cs","beta_el","sigma","mu_logA"))
h$label <- NA; h$label[grep("^beta_el", h$variable)] <- ELTYPES; h$label[grep("^sigma\\[", h$variable)] <- x$rdims; h$label[grep("^mu_logA", h$variable)] <- x$dims
h$label[grep("^(lambda|sd_phy|sd_res)\\[", h$variable)] <- rep(TRAITS, 3)
print(h, digits=3)
all <- fit$summary(NULL, "rhat", "ess_bulk"); cat(sprintf("\nworst R-hat %.3f | parameters with R-hat > 1.05: %d of %d | min bulk ESS %.0f\n",
    max(all$rhat, na.rm=TRUE), sum(all$rhat > 1.05, na.rm=TRUE), nrow(all), min(all$ess_bulk, na.rm=TRUE)))

# species traits
sp <- do.call(rbind, lapply(c("logr_sp","nu_sp","t50_sp"), function(v){ s <- sm(v); data.frame(trait=sub("_sp","",v), species=x$species, s[,-1], check.names=FALSE) }))
nobs <- bind_rows(x$rings %>% count(species), x$ec %>% filter(has_size) %>% count(species)) %>% group_by(species) %>% summarise(n_sized=sum(n))
sp <- left_join(sp, nobs, by="species") %>% left_join(x$inds %>% count(species, name="n_ind"), by="species")
write_csv(sp, file.path(out,paste0("age_element_species_",tag,".csv")))
cat("\n=== per-species shape nu (1 logistic, 0 Gompertz, -1/3 von Bertalanffy, -1 monomolecular) ===\n")
print(sp %>% filter(trait=="nu") %>% select(species, median, `5%`, `95%`, rhat, n_ind, n_sized), digits=2)

# individual ages and element losses
ti <- sm("Tage"); it <- cbind(x$inds[,c("species","individual","age_type","age_est","Lmax")], age_med=ti$median, age_lo=ti$`5%`, age_hi=ti$`95%`, rhat=ti$rhat)
write_csv(it, file.path(out,paste0("age_element_individual_ages_",tag,".csv")))
mi <- sm("miss"); la <- sm("logA_e")
et <- cbind(x$er[,c("species","individual","element","dimension","efs","eroded","last")], lost_inner_med=mi$median, lost_inner_lo=mi$`5%`, lost_inner_hi=mi$`95%`,
            A_med=exp(la$median), frac_attained=exp(x$er$logmax - la$median))
et <- left_join(et, it[,c("individual","age_med")], by="individual"); et$lost_outer_med <- pmax(0, et$age_med - et$last - et$lost_inner_med)
write_csv(et, file.path(out,paste0("age_element_ring_elements_",tag,".csv")))
cl <- sm("c_lost"); cr <- sm("c_rel")
ct <- cbind(x$ec[,c("species","individual","element","dimension","efs","count","size","nrep")], lost_med=cl$median, lost_lo=cl$`5%`, lost_hi=cl$`95%`, frac_attained=cr$median)
ct <- left_join(ct, it[,c("individual","age_med","age_lo","age_hi")], by="individual")
write_csv(ct, file.path(out,paste0("age_element_cs_elements_",tag,".csv")))
cat("\n=== individuals with several ring-series elements: one age, element-specific losses ===\n")
mm <- et %>% group_by(individual) %>% filter(n() > 1) %>% ungroup() %>%
  transmute(species, individual, element, efs=c("absent","present","unknown")[efs+1], eroded, marks=last, age=round(age_med,1),
            lost_inner=round(lost_inner_med,1), lost_outer=round(lost_outer_med,1), frac_A=round(frac_attained,2))
print(as.data.frame(mm))

# ancestry-vs-independent partition, as for model 26
S <- d$S; B <- d$B; P <- x$tm$P; tos <- x$tip_of_sp
delta <- get_arr(fit, "delta", B); dev <- get_arr(fit, "dev", S)
pp <- pair_partition(delta, dev, P, tos); ps <- partition_summary(pp, delta, dev, P, tos)
q <- function(m) { r <- t(apply(m, 2, quantile, c(.05,.5,.95))); rownames(r) <- TRAITS; r }
cat("\n=== realized tree share of among-species variance (5/50/95%) ===\n"); print(round(q(ps$phy_real),2))
pt <- do.call(rbind, lapply(seq_along(TRAITS), function(k) data.frame(trait=TRAITS[k],
  sp_i=x$species[pp$pairs[,1]], sp_j=x$species[pp$pairs[,2]], inherited=colMeans(pp$inherited[,,k]),
  independent=colMeans(pp$independent[,,k]), cross=colMeans(pp$cross[,,k]), p_parallel=colMeans(pp$independent[,,k] > 0))))
write_csv(pt, file.path(out,paste0("age_element_pairwise_",tag,".csv")))
cat("\nDONE model 32 fit\n")
