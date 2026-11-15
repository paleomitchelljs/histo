#!/usr/bin/env Rscript
# What identifies the channel exponents (beta_ch, M1) of a model 36/37 fit? Two checks, written 2026-11-15.
# (1) Direct evidence: within an individual, the slope of log(other circumference) on log(femur circumference) at
#     matched mark numbers (ring series), and across individuals for count-and-size elements, estimates beta
#     directly (femur circumference has beta = 1). Lost-mark offsets differ between elements, so this is rough.
# (2) Indirect evidence: posterior correlation of each circumference exponent with every other quantity (species
#     traits, lost marks, element offsets). Large correlations with the curve shape of species that have only that
#     channel mean the exponent is being set by trading off against that species' shape.
# (PSIS leave-one-species-out was tried and is unreliable here: Pareto k > 0.7 for nearly every species.)
# Usage: Rscript 37_channel_check.R [fit tag=r36] [model=36|37]
suppressMessages({ library(dplyr); library(posterior) })
a <- commandArgs(trailingOnly=TRUE); tag <- ifelse(is.na(a[1]), "r36", a[1]); model <- ifelse(is.na(a[2]), "36", a[2])
root  <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
out <- file.path(root, "outputs"); options(width=200)
source(file.path(root, "src/32_helpers.R"))   # ELTYPES
x <- readRDS(file.path(out, sprintf("m%s_index_%s.rds", model, tag)))
er <- x$er
rg <- x$rings %>% select(-any_of(c("species","individual","channel","eltype"))) %>%
  left_join(er %>% select(el_id, species, channel, individual), by="el_id")
fc <- rg %>% filter(channel=="femur|circumference_mm") %>% select(individual, age0, fsz=size)
oth <- rg %>% filter(grepl("circumference", channel), channel!="femur|circumference_mm") %>% select(species, individual, channel, age0, osz=size)
m <- inner_join(oth, fc, by=c("individual","age0"), relationship="many-to-many")
cat("ring marks matched by individual and mark number (slope of log other on log femur, individual intercepts):\n")
print(as.data.frame(m %>% group_by(channel, species) %>% summarise(n=n(), n_ind=n_distinct(individual),
  slope=if (n() > n_distinct(individual) + 1) coef(if (n_distinct(individual) > 1) lm(log(osz) ~ log(fsz) + factor(individual)) else lm(log(osz) ~ log(fsz)))[2] else NA, .groups="drop")), digits=3)
ec <- x$ec %>% filter(has_size)
w <- ec %>% filter(grepl("circumference", channel)) %>% group_by(species, individual, channel) %>% summarise(size=mean(size), .groups="drop")
f2 <- w %>% filter(channel=="femur|circumference_mm") %>% select(individual, fsz=size)
m2 <- inner_join(w %>% filter(channel!="femur|circumference_mm"), f2, by="individual")
cat("\ncount-and-size elements in individuals that also have a femur circumference (across individuals):\n")
print(as.data.frame(m2 %>% group_by(channel) %>% summarise(n=n(), species=paste(unique(species), collapse=","),
  slope=if (n()>2) coef(lm(log(size) ~ log(fsz)))[2] else NA, .groups="drop")), digits=3)
fit <- readRDS(file.path(out, sprintf("fit_m%s_%s.rds", model, tag)))
ch <- x$channels$channel
dm <- as_draws_matrix(fit$draws(c("beta_ch","logr_sp","nu_sp","t50_sp","miss","Tage","beta_el","mu_logA","theta","gam")))
dm <- dm[, !grepl("^theta\\[\\d+,[2-4]\\]", colnames(dm))]       # duplicates of logr_sp, nu_sp, t50_sp
ix <- function(v) as.integer(sub("^[^[]*\\[(\\d+).*$", "\\1", v))
nm <- colnames(dm)
pretty <- nm
for (v in c("logr_sp","nu_sp","t50_sp")) { j <- grep(paste0("^", v, "\\["), nm); pretty[j] <- paste(v, x$species[ix(nm[j])]) }
j <- grep("^miss\\[", nm); pretty[j] <- paste("miss", x$er$individual[ix(nm[j])], x$er$element[ix(nm[j])])
j <- grep("^theta\\[\\d+,1\\]", nm); pretty[j] <- paste("logA", x$species[ix(nm[j])])
j <- grep("^gam\\[", nm); pretty[j] <- paste("gam", x$gamma_keys[ix(nm[j])])
j <- grep("^Tage\\[", nm); pretty[j] <- paste("Tage", x$inds$individual[ix(nm[j])], x$inds$species[ix(nm[j])])
j <- grep("^beta_el\\[", nm); pretty[j] <- paste("beta_el", ELTYPES[ix(nm[j])])
for (c in grep("(tibia|ulna|radius|humerus)\\|circumference", ch)) {
  b <- dm[, sprintf("beta_ch[%d]", c)]
  r <- suppressWarnings(cor(b, dm)); r <- r[1, ]; r <- r[!grepl("^beta_ch", nm)]
  o <- order(-abs(r))[1:12]
  cat("\n==", ch[c], " median", round(median(b), 2), "\n"); print(data.frame(r=round(r[o], 2), quantity=pretty[!grepl("^beta_ch", nm)][o]), row.names=FALSE)
}
