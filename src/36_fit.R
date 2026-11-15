#!/usr/bin/env Rscript
# Fit model 36 (src/36_age_element_model.stan: model 35 + species x element offsets + non-centred nu +
# periosteal sizes not counted as marks) to the tables from src/31_stage_broad.R; writes the M3 rate summaries.
# Usage: Rscript 36_fit.R [iter=700] [chains=4] [tag=r36] [adapt_delta=0.95] [max_treedepth=11] [data prefix=r36]
suppressMessages({ library(cmdstanr); library(dplyr); library(readr); library(ape); library(posterior) })
a <- commandArgs(trailingOnly=TRUE)
iter <- ifelse(is.na(a[1]), 700L, as.integer(a[1])); chains <- ifelse(is.na(a[2]), 4L, as.integer(a[2]))
tag <- ifelse(is.na(a[3]), "r36", a[3]); ad <- ifelse(is.na(a[4]), 0.95, as.numeric(a[4])); mtd <- ifelse(is.na(a[5]), 11L, as.integer(a[5])); pre <- ifelse(is.na(a[6]), "r36", a[6])
root  <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
clean <- file.path(root,"data_clean"); out <- file.path(root,"outputs")
source(file.path(root,"src/23_tree_helpers.R")); source(file.path(root,"src/32_helpers.R")); source(file.path(root,"src/35_helpers.R")); source(file.path(root,"src/36_helpers.R"))
TRAITS <- c("logA","logr","nu","t50")
rd <- function(f) as.data.frame(read_csv(file.path(clean,paste0(pre,"_",f)), show_col_types=FALSE))
tmap <- rd("species_to_tip.csv")
x <- make_standata36(rd("rings.csv"), rd("elements.csv"), rd("individuals.csv"),
                     read.tree(file.path(clean,paste0(pre,"_tree.nwk"))), setNames(tmap$phylo_tip, tmap$species))
d <- x$standata
options(width=200)
cat(sprintf("S %d | tips %d | I %d (%d multi-element, %d with age estimate, %d explicit timing effects) | ring elements %d | marks %d | count-and-size %d (%d sized, %d counted, %d interval counts) | dims %d | channels %d (%d free exponents) | relative intercepts %d | sex species %d\n",
            d$S, d$T, d$I, d$IM, sum(d$age_type), d$IU, d$ER, d$N, d$EC, sum(d$c_has_size), sum(d$c_cnt_type > 0), sum(d$c_cnt_type == 2), d$D, d$C, d$NB, d$NO, d$NSX))
cat("\n=== channels ===\n"); print(x$channels)
cat(sprintf("\nspecies x element offsets: %d (reference type per species: %s)\n", d$NG,
    paste(names(x$gamma_refs), x$gamma_refs, sep="=", collapse="; ")))
cat("\nsex species:", paste(x$sex_species, collapse=", "), "\nclades (lost-mark scale):\n"); print(table(x$clade_sp))
mod <- cmdstan_model(file.path(root,"src/36_age_element_model.stan"))
set.seed(1)
fit <- mod$sample(data=d, chains=chains, parallel_chains=chains, iter_warmup=iter, iter_sampling=iter,
                  adapt_delta=ad, max_treedepth=mtd, refresh=100, seed=1, init=make_init36(x))
fit$save_object(file.path(out,paste0("fit_m36_",tag,".rds"))); saveRDS(x, file.path(out,paste0("m36_index_",tag,".rds")))

cat("\n=== diagnostics ===\n"); print(fit$diagnostic_summary())
sm <- function(v) as.data.frame(summarise_draws(fit$draws(v), "median", ~quantile(.x, c(.05,.95)), "rhat", "ess_bulk"))
cat("\n=== hyperparameters ===\n")
h <- sm(c("lambda","sd_phy","sd_res","mu_logr","mu_nu","mu_t50","sd_logA_ind","rho_ind","sd_u","sd_gam","tau_off","sd_tel","sd_tcl","sigma_cs","gbar_sex","sd_sex","beta_el","sigma","mu_logA"))
h$label <- NA; h$label[grep("^beta_el", h$variable)] <- ELTYPES; h$label[grep("^sigma\\[", h$variable)] <- x$rdims; h$label[grep("^mu_logA", h$variable)] <- x$dims
h$label[grep("^(lambda|sd_phy|sd_res)\\[", h$variable)] <- rep(TRAITS, 3)
print(h, digits=3)
all <- fit$summary(NULL, "rhat", "ess_bulk"); cat(sprintf("\nworst R-hat %.3f | parameters with R-hat > 1.05: %d of %d | min bulk ESS %.0f\n",
    max(all$rhat, na.rm=TRUE), sum(all$rhat > 1.05, na.rm=TRUE), nrow(all), min(all$ess_bulk, na.rm=TRUE)))
cat("\n=== M1 channel exponents (1 = isometric with the latent curve) ===\n")
b <- sm("beta_ch"); bt <- cbind(x$channels, b[,-1]); print(bt, digits=3); write_csv(bt, file.path(out,paste0("m36_channels_",tag,".csv")))
cat("\n=== M7 lost-mark scale multipliers ===\n")
te <- sm("tau_el"); te$group <- TAU_ELTYPES; tc <- sm("tau_cl"); tc$group <- x$clades; print(rbind(te, tc), digits=3)
if (d$NG > 0) { cat("\n=== (a) species x element offsets (log asymptote, relative to the species' reference type) ===\n")
  g <- sm("gam"); g$key <- x$gamma_keys; print(g[, c("key","median","5%","95%","rhat")], digits=3) }
if (d$NSX > 0) { cat("\n=== M6 sex effect on log asymptote (male - female) ===\n"); g <- sm("g_sex"); g$species <- x$sex_species; print(g, digits=3) }

# species traits
sp <- do.call(rbind, lapply(c("logr_sp","nu_sp","t50_sp"), function(v){ s <- sm(v); data.frame(trait=sub("_sp","",v), species=x$species, s[,-1], check.names=FALSE) }))
nobs <- bind_rows(x$rings %>% count(species), x$ec %>% filter(has_size) %>% count(species)) %>% group_by(species) %>% summarise(n_sized=sum(n))
sp <- left_join(sp, nobs, by="species") %>% left_join(x$inds %>% count(species, name="n_ind"), by="species")
write_csv(sp, file.path(out,paste0("m36_species_",tag,".csv")))
cat("\n=== per-species traits ===\n"); print(sp %>% select(trait, species, median, `5%`, `95%`, rhat, n_ind, n_sized), digits=2)

# individuals
ti <- sm("Tage"); tu <- sm("u_ind")
it <- cbind(x$inds[,c("species","individual","age_type","age_est","Lmax","sex")], explicit_u=x$explicit_u, age_med=ti$median, age_lo=ti$`5%`, age_hi=ti$`95%`, age_rhat=ti$rhat,
            u_med=tu$median, u_lo=tu$`5%`, u_hi=tu$`95%`)
write_csv(it, file.path(out,paste0("m36_individual_ages_",tag,".csv")))
mi <- sm("miss"); la <- sm("logA_e")
et <- cbind(x$er[,c("species","individual","element","channel","efs","eroded","last")], lost_inner_med=mi$median, lost_inner_lo=mi$`5%`, lost_inner_hi=mi$`95%`,
            A_med=exp(la$median), frac_attained=exp(x$er$logmax - la$median))
write_csv(et, file.path(out,paste0("m36_ring_elements_",tag,".csv")))
cl <- sm("c_lost"); cr <- sm("c_rel")
ct <- cbind(x$ec[,c("species","individual","element","channel","efs","count","count_lo","count_hi","size","nrep")], lost_med=cl$median, lost_lo=cl$`5%`, lost_hi=cl$`95%`, frac_attained=cr$median)
write_csv(ct, file.path(out,paste0("m36_cs_elements_",tag,".csv")))

# ---- M3 rate summaries (generated in R from the posterior draws) ----
# logr_sp is the mid-growth rate of the latent curve = the reference channel (femur circumference, beta = 1).
# For species with an observed femur circumference: asymptotic femur circumference A_fc, body mass by
# Anderson et al. 1985 eq. 9 (bipeds, W[g] = 0.16 C_f^2.73; applied to every taxon: ROUGH), the
# maximum absolute mass growth rate, and the radial apposition rate dr/dt = (1/2pi) dC/dt at t50.
rich_log_R <- function(nu, xx) { ex <- exp(pmin(xx, 50)); u <- nu*ex; w <- 1 + u
  ifelse(abs(u) < 1e-4, -ex*(1 - 0.5*u + u^2/3), -log(0.5*(w + sqrt(w^2 + 1e-12)))/nu) }
gf <- function(f, nu) ifelse(abs(nu*log(f)) < 1e-8, log(-log(f)), log(expm1(-nu*log(f))/nu))
dr <- as_draws_matrix(fit$draws(c("theta","mu_logA","logr_sp","nu_sp","t50_sp")))
dr <- dr[round(seq(1, nrow(dr), length.out=min(400, nrow(dr)))), ]
fc_dim <- match("circumference_mm", x$dims)
fc_sp <- unique(c(x$er$species[x$er$channel == REF_CHANNEL], x$ec$species[x$ec$channel == REF_CHANNEL & x$ec$has_size]))
tg <- seq(0, 60, by=0.05)
m3 <- do.call(rbind, lapply(seq_along(x$species), function(s){
  nu <- dr[, sprintf("nu_sp[%d]", s)]; lr <- dr[, sprintf("logr_sp[%d]", s)]; t50 <- dr[, sprintf("t50_sp[%d]", s)]
  k <- 2*exp(lr)*(gf(0.25, nu) - gf(0.75, nu)); g50 <- gf(0.5, nu)
  out <- data.frame(species=x$species[s], logr_med=median(lr), r_med=median(exp(lr)), has_femur_circ=x$species[s] %in% fc_sp,
                    A_fc_mm=NA, mass_kg=NA, Gmax_kg_yr=NA, Gmax_lo=NA, Gmax_hi=NA, app_um_day=NA, app_lo=NA, app_hi=NA)
  if (out$has_femur_circ) {
    A <- exp(dr[, sprintf("mu_logA[%d]", fc_dim)] + dr[, sprintf("theta[%d,1]", s)])
    G <- app <- numeric(nrow(dr))
    for (j in seq_len(nrow(dr))) { rel <- exp(rich_log_R(nu[j], -k[j]*(tg - t50[j]) + g50[j]))
      M <- 0.16 * (A[j]*rel)^2.73 / 1000; G[j] <- max(diff(M))/0.05
      i50 <- which.min(abs(tg - t50[j])); i50 <- max(2, min(i50, length(tg) - 1))
      app[j] <- A[j] * (rel[i50+1] - rel[i50-1]) / 0.1 / (2*pi) * 1000 / 365.25 }
    q <- function(v) quantile(v, c(.5,.05,.95), na.rm=TRUE)
    out$A_fc_mm <- median(A); out$mass_kg <- median(0.16*A^2.73/1000)
    out[, c("Gmax_kg_yr","Gmax_lo","Gmax_hi")] <- q(G); out[, c("app_um_day","app_lo","app_hi")] <- q(app) }
  out }))
write_csv(m3, file.path(out,paste0("m36_rates_",tag,".csv")))
cat("\n=== M3 rates (femur-circumference species: mass by Anderson 1985 eq. 9, rough) ===\n"); print(m3, digits=3)

# ancestry-vs-independent partition, as for models 26 / 32
S <- d$S; B <- d$B; P <- x$tm$P; tos <- x$tip_of_sp
delta <- get_arr(fit, "delta", B); dev <- get_arr(fit, "dev", S)
pp <- pair_partition(delta, dev, P, tos); ps <- partition_summary(pp, delta, dev, P, tos)
q <- function(m) { r <- t(apply(m, 2, quantile, c(.05,.5,.95))); rownames(r) <- TRAITS; r }
cat("\n=== realized tree share of among-species variance (5/50/95%) ===\n"); print(round(q(ps$phy_real),2))
cat("\nDONE model 36 fit\n")
