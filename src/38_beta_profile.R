#!/usr/bin/env Rscript
# Which observations want a high channel exponent? For every element that uses a free exponent (r_bi / c_bi > 0) in
# the chosen channels: log likelihood of its sizes over a grid of beta, every other quantity held at each posterior
# draw (conditional, not marginal), averaged over draws. Reports per element the beta that maximises the mean
# conditional log likelihood, and the drop in log likelihood at beta = 1 (how hard the element pulls away from
# isometry, given the rest of the fit). EFS terms are ignored. Written 2026-11-16 (docs/allometry_prior_plan.md).
# Usage: Rscript 38_beta_profile.R [fit tag=r38] [model=37] [channel regex=circumference]
suppressMessages({ library(dplyr); library(posterior) })
a <- commandArgs(trailingOnly=TRUE); tag <- ifelse(is.na(a[1]), "r38", a[1]); model <- ifelse(is.na(a[2]), "37", a[2])
rx <- ifelse(is.na(a[3]), "(tibia|radius|ulna|humerus)\\|circumference", a[3])
root  <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
out <- file.path(root, "outputs"); options(width=200)
fit <- readRDS(file.path(out, sprintf("fit_m%s_%s.rds", model, tag))); x <- readRDS(file.path(out, sprintf("m%s_index_%s.rds", model, tag)))
d <- x$standata; ch <- x$channels$channel
rich_log_R <- function(nu, xx) { ex <- exp(pmin(xx, 50)); u <- nu*ex; w <- 1 + u
  ifelse(abs(u) < 1e-4, -ex*(1 - 0.5*u + u^2/3), -log(0.5*(w + sqrt(w^2 + 1e-12)))/nu) }
rich_dlog_dx_R <- function(nu, xx) { ex <- exp(pmin(xx, 50)); -ex / pmax(1 + nu*ex, 1e-3) }
dm <- as_draws_matrix(fit$draws(c("logA_e","miss","t50_i","nu_sp","k_sp","g50_sp","beta_ch","sigma","sigma_cs","theta","mu_logA",
                                  "beta_el","gam","o_rel","a_ind","sd_el","Tage","sd_u")))
dm <- dm[round(seq(1, nrow(dm), length.out=400)), ]
g <- function(v, i=NULL) { n <- if (is.null(i)) v else sprintf("%s[%s]", v, i); if (n %in% colnames(dm)) dm[, n] else 0 }
grid <- seq(0.4, 3.5, by=0.05)
prof <- function(llfun) { m <- vapply(grid, function(b) mean(llfun(b)), 0); c(best=grid[which.max(m)], drop_at_1=max(m) - m[grid == 1]) }
rows <- list()
for (e in which(d$r_bi > 0 & grepl(rx, ch[d$r_ch]))) {
  i <- d$r_ind[e]; s <- d$sp_of_ind[i]; n <- which(d$el == e)
  lr <- sapply(n, function(j) rich_log_R(g("nu_sp", s), -g("k_sp", s)*(d$age0[j] + g("miss", e) - g("t50_i", i)) + g("g50_sp", s)))
  ll <- function(b) rowSums(sapply(seq_along(n), function(q) dlnorm(d$size[n[q]], g("logA_e", e) + b*lr[, q], g("sigma", d$r_sdim[e]), log=TRUE)))
  rows[[length(rows)+1]] <- data.frame(kind="ring", channel=ch[d$r_ch[e]], species=x$er$species[e], individual=x$er$individual[e],
    n_sizes=length(n), beta_post=median(g("beta_ch", d$r_ch[e])), t(prof(ll)))
}
for (e in which(d$c_bi > 0 & d$c_has_size == 1 & grepl(rx, ch[d$c_ch]))) {
  i <- d$c_ind[e]; s <- d$sp_of_ind[i]
  xx <- -g("k_sp", s)*(g("Tage", i) + d$c_tadd[e] - g("t50_i", i)) + g("g50_sp", s); lr <- rich_log_R(g("nu_sp", s), xx)
  bel <- (if (d$c_elt[e] > 0) g("beta_el", d$c_elt[e]) else 0) + (if (d$c_g[e] > 0) g("gam", d$c_g[e]) else 0)
  base <- if (d$c_relch[e] > 0) g("o_rel", d$c_relch[e]) else g("mu_logA", d$c_dim[e]) + dm[, sprintf("theta[%d,1]", s)] + bel
  ll <- function(b) { vu <- if (d$u_of_ind[i] == 0) (b*g("k_sp", s)*rich_dlog_dx_R(g("nu_sp", s), xx)*g("sd_u"))^2 / d$c_nrep[e] else 0
    dnorm(d$c_logsize[e], base + g("a_ind", i) + b*lr, sqrt(g("sd_el", i)^2/d$c_nrep[e] + g("sigma_cs")^2 + d$c_logsd[e]^2 + vu), log=TRUE) }
  rows[[length(rows)+1]] <- data.frame(kind="count+size", channel=ch[d$c_ch[e]], species=x$ec$species[e], individual=x$ec$individual[e],
    n_sizes=1, beta_post=median(g("beta_ch", d$c_ch[e])), t(prof(ll)))
}
res <- bind_rows(rows) %>% arrange(channel, desc(drop_at_1))
res$beta_post <- round(res$beta_post, 2); res$drop_at_1 <- round(res$drop_at_1, 2)
print(as.data.frame(res), row.names=FALSE)
