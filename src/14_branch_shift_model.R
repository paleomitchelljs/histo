#!/usr/bin/env Rscript
# Fit the branch-shift model (14_branch_shift_model.stan): per-branch horseshoe
# increments to nu on the expanded cross-Dinosauria tree. Reports which branches
# carry shape shifts.
suppressMessages({ library(cmdstanr); library(dplyr); library(readr); library(ape); library(posterior) })
root  <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
clean <- file.path(root,"data_clean"); out <- file.path(root,"outputs")

dat <- read_csv(file.path(clean,"lag_size_long.csv"), show_col_types=FALSE) %>% filter(!is.na(size), !is.na(age))
dat$species <- factor(dat$species); dat$specimen <- factor(dat$specimen); dat$dimr <- as.integer(factor(dat$dimension))
ind_tbl <- dat %>% distinct(specimen, species, resorption, efs_present) %>%
  left_join(dat %>% group_by(specimen) %>% summarise(logmax=log(max(size)), .groups="drop"), by="specimen") %>%
  arrange(specimen); ind_tbl$sp_of_ind <- as.integer(ind_tbl$species)

tree <- read.tree(file.path(clean,"dino_tree.nwk"))
edge <- tree$edge; B <- nrow(edge); T <- length(tree$tip.label)
P <- matrix(0, T, B)                                  # root->tip path matrix
for (i in seq_len(T)) { node <- i
  repeat { e <- which(edge[,2]==node); if (length(e)==0) break; P[i,e] <- 1; node <- edge[e,1] } }
tips <- tree$tip.label
tipmap <- read_csv(file.path(clean,"species_to_tip.csv"), show_col_types=FALSE)
sp_levels <- levels(dat$species)
tip_of_sp <- match(tipmap$phylo_tip[match(sp_levels, tipmap$species)], tips)
stopifnot(!any(is.na(tip_of_sp)))

standata <- list(N=nrow(dat), S=nlevels(dat$species), I=nlevels(dat$specimen), D=max(dat$dimr),
  T=T, B=B, size=dat$size, age0=dat$age, ind=as.integer(dat$specimen), dimr=dat$dimr,
  sp_of_ind=ind_tbl$sp_of_ind, tip_of_sp=tip_of_sp, P=P,
  efs=as.integer(ind_tbl$efs_present), resorb=as.integer(ind_tbl$resorption!="absent"),
  logmax=ind_tbl$logmax, slab_scale=0.8, tau0=0.2)

mod <- cmdstan_model(file.path(root,"src/14_branch_shift_model.stan"))
fit <- mod$sample(data=standata, chains=4, parallel_chains=1, iter_warmup=1000, iter_sampling=1000,
                  adapt_delta=0.95, max_treedepth=12, refresh=250, seed=1, init=0.2)
fit$save_object(file.path(out,"fit_branch_shift.rds"))
# branch -> (parent,child) + descendant tips, for interpreting shifts
edge_df <- data.frame(branch=seq_len(B), parent=edge[,1], child=edge[,2])
saveRDS(list(species=sp_levels, specimen=levels(dat$specimen), ind_tbl=ind_tbl,
             tips=tips, tip_of_sp=tip_of_sp, edge=edge, P=P, tree=tree), file.path(out,"branch_shift_index.rds"))

cat("\n=== diagnostics ===\n"); print(fit$diagnostic_summary())
cat("\n=== branches with shape shifts (delta posterior excluding 0) ===\n")
d <- fit$draws("delta_out", format="df")
shifts <- do.call(rbind, lapply(seq_len(B), function(j){
  v <- d[[paste0("delta_out[",j,"]")]]
  data.frame(branch=j, child=edge[j,2], med=median(v), lo=quantile(v,.05), hi=quantile(v,.95),
             p_pos=mean(v>0))}))
# descendant tip labels per branch (child node)
desc <- function(node){ if(node<=T) return(tips[node]); ch<-edge[edge[,1]==node,2]; paste(unlist(lapply(ch,desc)),collapse=",") }
shifts$clade <- sapply(shifts$child, function(c) substr(desc(c),1,40))
sig <- shifts[shifts$lo>0 | shifts$hi<0,]
sig <- sig[order(-abs(sig$med)),]
if (nrow(sig)) print(sig, row.names=FALSE) else cat("(no branch's 90% CrI excludes 0 -- shape shifts are diffuse, consistent with high phylo signal)\n")
cat("\nDONE branch-shift model\n")
