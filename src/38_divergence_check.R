#!/usr/bin/env Rscript
# Where are the divergent transitions of a model 36/37 fit? For every sampled quantity: shift of its mean in the
# divergent draws relative to the rest, in posterior SDs, and the same for the tree-scale block alone.
# Written 2026-11-16 (question: do poorly constrained phylogenetic quantities drive the divergences?).
# Usage: Rscript 38_divergence_check.R [fit tag=r38] [model=37]
suppressMessages({ library(dplyr); library(posterior) })
a <- commandArgs(trailingOnly=TRUE); tag <- ifelse(is.na(a[1]), "r38", a[1]); model <- ifelse(is.na(a[2]), "37", a[2])
root  <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
out <- file.path(root, "outputs"); options(width=200)
fit <- readRDS(file.path(out, sprintf("fit_m%s_%s.rds", model, tag))); x <- readRDS(file.path(out, sprintf("m%s_index_%s.rds", model, tag)))
sd_ <- as_draws_df(fit$sampler_diagnostics()); div <- sd_$divergent__ == 1
cat(sprintf("divergent: %d of %d; by chain: %s\n", sum(div), length(div), paste(tapply(div, sd_$.chain, sum), collapse=" ")))
pars <- c("mu_logA","mu_logr","mu_nu","mu_t50","sd_phy","sd_res","L_Om","kappa_raw","theta_o","z_nu","nu_sp","sd_gam","beta_el",
          "sd_logA_ind","rho_ind","log_tau0","sd_tel","sd_tcl","sigma","sigma_cs","log_beta","sd_u","gbar_sex","sd_sex","lambda",
          "Tage","miss","logA_e","t_raw","f_ero","o_rel","z_u","z_gam","z_sex","z_ind")
dm <- as_draws_matrix(fit$draws(intersect(pars, fit$metadata()$stan_variables)))
dm <- dm[, apply(dm, 2, sd) > 0, drop=FALSE]
shift <- (colMeans(dm[div, , drop=FALSE]) - colMeans(dm[!div, , drop=FALSE])) / apply(dm, 2, sd)
ix <- function(v) as.integer(sub("^[^[]*\\[(\\d+).*$", "\\1", v))
lab <- names(shift)
for (v in c("nu_sp","z_nu")) { j <- grep(paste0("^", v, "\\["), lab); lab[j] <- paste(lab[j], x$species[ix(lab[j])]) }
j <- grep("^theta_o\\[", lab); lab[j] <- paste(lab[j], x$species[ix(lab[j])], c("logA","logr","t50")[as.integer(sub(".*,(\\d+)\\]", "\\1", lab[j]))])
for (v in c("Tage","t_raw")) { j <- grep(paste0("^", v, "\\["), lab); lab[j] <- paste(lab[j], x$inds$individual[ix(lab[j])], x$inds$species[ix(lab[j])]) }
for (v in c("miss","logA_e")) { j <- grep(paste0("^", v, "\\["), lab); lab[j] <- paste(lab[j], x$er$individual[ix(lab[j])], x$er$element[ix(lab[j])]) }
j <- grep("^(lambda|sd_phy|sd_res)\\[", lab); lab[j] <- paste(lab[j], c("logA","logr","nu","t50")[ix(lab[j])])
res <- data.frame(quantity=lab, shift_sd=round(shift, 2), row.names=NULL)
cat("\n== largest shifts, all quantities (divergent minus rest, posterior SDs) ==\n"); print(head(res[order(-abs(res$shift_sd)), ], 25), row.names=FALSE)
tree <- grepl("^(sd_phy|sd_res|lambda|L_Om|kappa_raw|mu_)", res$quantity)
cat("\n== tree-scale block ==\n"); print(res[tree, ][order(-abs(res$shift_sd[tree])), ], row.names=FALSE)
grp <- sub("\\[.*", "", names(shift))
cat("\n== mean |shift| by parameter group ==\n"); print(round(sort(tapply(abs(shift), grp, mean), decreasing=TRUE), 2))
# where do divergences sit in the funnel? species-trait spread vs its scale, divergent vs rest
for (k in 1:4) { s1 <- sprintf("sd_phy[%d]", k); s2 <- sprintf("sd_res[%d]", k)
  if (all(c(s1, s2) %in% colnames(dm))) cat(sprintf("%s: sd_phy median div %.3f vs rest %.3f | sd_res div %.3f vs rest %.3f\n",
    c("logA","logr","nu","t50")[k], median(dm[div, s1]), median(dm[!div, s1]), median(dm[div, s2]), median(dm[!div, s2]))) }
