#!/usr/bin/env Rscript
# Evaluate the branch-shift model (14) on the SHIFT simulation scenario:
# does the regularized-horseshoe per-branch shift recover the TWO true shift
# edges (recorded in scenario_shift_truth_edges.csv)?
suppressMessages({ library(cmdstanr); library(dplyr); library(readr); library(ape); library(posterior) })
root <- "/Users/jmitchell/Library/CloudStorage/Dropbox/Research/00_current_Histology"
clean <- file.path(root,"data_clean"); sim <- file.path(clean,"sim"); out <- file.path(root,"outputs")

dat <- read_csv(file.path(sim,"scenario_shift_data.csv"), show_col_types=FALSE)
truth_edges <- read_csv(file.path(sim,"scenario_shift_truth_edges.csv"), show_col_types=FALSE)
dat$species <- factor(dat$species); dat$specimen <- factor(dat$specimen); dat$dimr <- 1L
it <- dat %>% distinct(specimen, species, resorption, efs_present) %>%
  left_join(dat %>% group_by(specimen) %>% summarise(logmax=log(max(size)), .groups="drop"), by="specimen") %>%
  arrange(specimen); it$sp_of_ind <- as.integer(it$species)
spl <- levels(dat$species)
tree <- compute.brlen(keep.tip(read.tree(file.path(clean,"dino_tree.nwk")), spl), method="Grafen")
tips <- tree$tip.label; T <- length(tips); edges <- tree$edge; B <- nrow(edges); rootnode <- T+1
# path matrix P[t,b] = 1 if edge b on root->tip t path
P <- matrix(0, T, B)
for (t in 1:T) { pth <- nodepath(tree, from=rootnode, to=t)
  for (k in 1:(length(pth)-1)) { e <- which(edges[,1]==pth[k] & edges[,2]==pth[k+1]); P[t,e] <- 1 } }
tip_of_sp <- match(spl, tips)

standata <- list(N=nrow(dat), S=nlevels(dat$species), I=nlevels(dat$specimen), D=1L, T=T, B=B,
  size=dat$size, age0=dat$age, ind=as.integer(dat$specimen), dimr=dat$dimr, sp_of_ind=it$sp_of_ind,
  tip_of_sp=tip_of_sp, P=P, efs=as.integer(it$efs_present), resorb=as.integer(it$resorption!="absent"),
  logmax=it$logmax, slab_scale=2.0, tau0=0.30)   # much weaker shrinkage: can real shifts emerge?

mod <- cmdstan_model(file.path(root,"src/14_branch_shift_model.stan"))
csvdir <- file.path(out,"bs_csv"); dir.create(csvdir, showWarnings=FALSE)
fit <- mod$sample(data=standata, chains=1, iter_warmup=400, iter_sampling=400,
                  adapt_delta=0.99, max_treedepth=9, refresh=200, seed=1, init=0.1, output_dir=csvdir)
fit$save_object(file.path(out,"fit_sim_branchshift.rds"))
d <- fit$draws(format="df")

# per-branch shift summary
delta <- data.frame(edge=1:B,
  med=sapply(1:B,function(b) median(d[[paste0("delta_out[",b,"]")]])),
  absmed=sapply(1:B,function(b) median(abs(d[[paste0("delta_out[",b,"]")]]))),
  p_big=sapply(1:B,function(b) mean(abs(d[[paste0("delta_out[",b,"]")]]) > 0.5)),
  child_node=edges[,2])
delta$is_true_shift <- delta$edge %in% truth_edges$edge
delta <- delta %>% left_join(truth_edges %>% transmute(edge, true_jump=jump), by="edge")
delta_rank <- delta %>% arrange(desc(absmed))
write_csv(delta, file.path(out,"sim_branchshift_delta.csv"))

cat("\n===== BRANCH-SHIFT RECOVERY (shift scenario) =====\n")
cat("true shift edges:", paste(truth_edges$edge,collapse=","), "| true jumps:", paste(round(truth_edges$jump,2),collapse=","), "\n")
cat("divergences:", sum(fit$diagnostic_summary()$num_divergent), "\n\n")
cat("--- top 6 branches by |median shift| ---\n")
print(head(delta_rank %>% transmute(edge, child_node, median_shift=round(med,2), abs_med=round(absmed,2),
  p_excl=round(p_big,2), TRUE_SHIFT=ifelse(is_true_shift,"<== TRUE","")), 6), row.names=FALSE)
cat("\n--- the true shift edges ---\n")
print(delta %>% filter(is_true_shift) %>% transmute(edge, true_jump=round(true_jump,2),
  median_shift=round(med,2), p_abs_gt_0.5=round(p_big,2),
  rank_by_absmed=match(edge, delta_rank$edge)), row.names=FALSE)
cat("\nDONE branch-shift sim\n")
