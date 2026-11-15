#!/usr/bin/env Rscript
# =====================================================================
# 18_simulate.R  -- parameter-recovery / SBC simulator for the growth model.
# Generates LAG-by-LAG data under three generative regimes for the SHAPE
# parameter lognu (the thing the phylo/shift machinery targets):
#   (1) iid    -- curves vary randomly, NO phylogenetic structure
#   (2) BM     -- curves vary along the tree by Brownian motion
#   (3) shift  -- a few large step-changes on internal branches
# In ALL regimes: realistic Richards curves, RANDOM neonate masses, and a KNOWN
# per-specimen resorption offset (inner rings absorbed before the first LAG).
# Truth is saved so we can score recovery of curve params, neonate mass, and the
# missing-LAG offset, and check whether the model correctly reads the regime.
#
# Output: data_clean/sim/scenario_{iid,bm,shift}_{data,truth_sp,truth_ind}.csv
# =====================================================================
suppressMessages({ library(ape); library(phytools); library(dplyr); library(readr) })
set.seed(20260609)
root  <- "/Users/jmitchell/Library/CloudStorage/Dropbox/Research/00_current_Histology"
clean <- file.path(root,"data_clean"); simdir <- file.path(clean,"sim"); dir.create(simdir, showWarnings=FALSE)

tree  <- compute.brlen(read.tree(file.path(clean,"dino_tree.nwk")), method="Grafen")
tips  <- tree$tip.label; T <- length(tips)
N_PER <- 2                      # specimens per species
SIGMA <- 0.05                   # lognormal obs noise (matches real fits)
# realistic centres / scales (from the real fits)
MU_LOGNU <- -0.4; SD_NU <- 1.0  # shape: nu ~ 0.1-2.5
MU_LOGK  <- -1.1; SD_LOGK <- 0.35
MU_T0    <- 2.0;  SD_T0 <- 1.2
MU_LOGA  <- 5.6;  SD_LOGA_SP <- 0.45; SD_LOGA_IND <- 0.20

richards_logM <- function(t, logA, nu, k, t0) logA - (1/nu)*log1p(nu*exp(-k*(t - t0)))

# scenario-specific lognu per tip (the regime lives here)
gen_lognu <- function(scn){
  if (scn=="iid")  return(setNames(rnorm(T, MU_LOGNU, SD_NU), tips))
  if (scn=="bm") {                                        # BM, rescaled to SD_NU + clear signal
    repeat {
      raw <- fastBM(tree, sig2=1, a=0)
      x <- MU_LOGNU + SD_NU * (raw - mean(raw)) / sd(raw)
      if (phylosig(tree, setNames(x, tips), method="lambda")$lambda > 0.7) return(setNames(x, tips))
    }
  }
  if (scn=="shift"){
    v <- setNames(rep(MU_LOGNU, T), tips)
    # two large clade shifts on internal nodes (descendants get a big jump)
    inodes <- (T+1):(T+tree$Nnode); set.seed(7); shifts <- list()
    for (jmp in c(+1.8, -1.6)) {                         # big steps in lognu
      nd <- sample(inodes[sapply(inodes,function(n) length(getDescendants(tree,n))>=3 & length(getDescendants(tree,n))<=8)],1)
      desc_tips <- intersect(getDescendants(tree,nd), 1:T)
      v[tips[desc_tips]] <- v[tips[desc_tips]] + jmp
      shifts[[length(shifts)+1]] <- data.frame(node=nd, edge=which(tree$edge[,2]==nd),
        jump=jmp, n_desc_tips=length(desc_tips))
    }
    attr(v,"shifts") <- bind_rows(shifts)
    return(v)
  }
}

simulate_one <- function(scn){
  lognu_sp <- gen_lognu(scn)
  logk_sp  <- setNames(rnorm(T, MU_LOGK, SD_LOGK), tips)   # rate: iid realistic in all regimes
  t0_sp    <- setNames(rnorm(T, MU_T0, SD_T0), tips)
  logA_sp  <- setNames(rnorm(T, MU_LOGA, SD_LOGA_SP), tips)
  rows <- list(); tr_ind <- list(); ii <- 0
  for (s in tips) {
    nu <- exp(lognu_sp[s]); k <- exp(logk_sp[s]); t0 <- t0_sp[s]
    for (r in 1:N_PER) {
      ii <- ii+1; spec <- paste0(s,"_",r)
      logA_i <- logA_sp[s] + rnorm(1,0,SD_LOGA_IND)
      death  <- sample(7:22,1)                              # true age at death
      resorbed <- rbinom(1,1,0.6)                           # 60% have inner-LAG resorption
      offset <- if (resorbed) min(sample(1:7,1), death-3) else 0   # KNOWN missing inner rings
      efs    <- as.integer(death >= 15 & runif(1)<0.6)      # some reach EFS (growth ceased)
      neonate <- exp(richards_logM(0, logA_i, nu, k, t0))   # TRUE neonate mass (size at age 0)
      # observed rings: recorded index j=1.., true age = j+offset, up to death.
      # perinate-anchored (offset 0) also records age 0.
      true_ages <- if (offset==0) c(0, 1:death) else (offset+1):death
      size <- exp(sapply(true_ages, richards_logM, logA=logA_i, nu=nu, k=k, t0=t0) + rnorm(length(true_ages),0,SIGMA))
      is_efs <- (true_ages >= death-2) & efs==1
      rows[[ii]] <- data.frame(scenario=scn, species=s, specimen=spec,
        lag_index=seq_along(true_ages)-as.integer(offset==0),  # recorded ring index
        age=true_ages-offset,                                  # recorded age (= true_age - offset)
        true_age=true_ages, size=size, is_efs=is_efs,
        resorption=ifelse(resorbed,"present","absent"), efs_present=efs==1, dimension="sim_mm")
      tr_ind[[ii]] <- data.frame(scenario=scn, species=s, specimen=spec,
        logA=logA_i, offset_true=offset, neonate_true=neonate, death=death, efs=efs)
    }
  }
  list(data=bind_rows(rows),
       truth_sp=data.frame(scenario=scn, species=tips, lognu=as.numeric(lognu_sp[tips]), logk=logk_sp[tips],
                           t0=t0_sp[tips], logA_sp=logA_sp[tips], nu=exp(as.numeric(lognu_sp[tips]))),
       truth_ind=bind_rows(tr_ind), shifts=attr(lognu_sp,"shifts"))
}

for (scn in c("iid","bm","shift")) {
  s <- simulate_one(scn)
  write_csv(s$data,      file.path(simdir, paste0("scenario_",scn,"_data.csv")))
  write_csv(s$truth_sp,  file.path(simdir, paste0("scenario_",scn,"_truth_sp.csv")))
  write_csv(s$truth_ind, file.path(simdir, paste0("scenario_",scn,"_truth_ind.csv")))
  if (!is.null(s$shifts)) { write_csv(s$shifts, file.path(simdir,"scenario_shift_truth_edges.csv"))
    cat("  true shift edges:", paste(s$shifts$edge, collapse=","), "| nodes:", paste(s$shifts$node,collapse=","), "\n") }
  ps <- phylosig(tree, setNames(s$truth_sp$lognu, s$truth_sp$species), method="lambda")
  cat(sprintf("[%s] %d obs | %d specimens | true nu range %.2f-%.2f | true lognu Pagel-lambda=%.2f\n",
      scn, nrow(s$data), n_distinct(s$data$specimen), min(s$truth_sp$nu), max(s$truth_sp$nu), ps$lambda))
}
cat("\nwrote data_clean/sim/scenario_*_{data,truth_sp,truth_ind}.csv\n")
