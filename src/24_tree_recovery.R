#!/usr/bin/env Rscript
# Simulation / recovery study for the tree-nested model (23).
# Regimes are applied (independently) to logk and lognu on the real 22-tip
# topology with unit branch lengths; logA and t0 are iid across species.
#   iid   -- species differ, no tree structure          (all change is lineage-private)
#   bm    -- change accumulates on every branch         (diffuse shared ancestry)
#   shift -- one jump on an internal branch (4-9 tips)  (similarity by shared ancestry)
#   conv  -- the SAME jump on two distant small clades  (similarity by independent acquisition)
# Scored against known truth: lambda, realized tree share of variance, ancestry
# share of pairwise similarity, pairwise inherited / independent terms, offsets.
# Usage: Rscript 24_tree_recovery.R <nreps=10> <chunk=0> <n_per=2> [regimes=iid,bm,shift,conv]
suppressMessages({ library(cmdstanr); library(ape); library(dplyr); library(readr); library(posterior) })
a <- commandArgs(trailingOnly=TRUE)
nreps <- ifelse(is.na(a[1]), 10L, as.integer(a[1])); chunk <- ifelse(is.na(a[2]), 0L, as.integer(a[2]))
N_PER <- ifelse(is.na(a[3]), 2L, as.integer(a[3]))
regimes <- if (is.na(a[4])) c("iid","bm","shift","conv") else strsplit(a[4], ",")[[1]]
set.seed(20261113 + 1000*chunk + 17*N_PER)
root  <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
clean <- file.path(root,"data_clean"); out <- file.path(root,"outputs")
source(file.path(root,"src/23_tree_helpers.R"))
resfile <- file.path(out,"tree_recovery_results.csv")

tree <- read.tree(file.path(clean,"dino_tree.nwk")); tm <- tree_mats(tree)
P <- tm$P; B <- tm$B; T <- tm$T; tips <- tm$tips; edge <- tm$edge
ndesc <- colSums(P)                                         # tips subtended by each branch
SIGMA<-0.05; MU <- c(logA=5.6, logk=-1.1, lognu=-0.4, t0=2.0); SDV <- c(logA=0.45, logk=0.35, lognu=1.0, t0=1.2)
SD_LOGA_IND <- 0.20
rich <- function(t,logA,nu,k,t0) logA-(1/nu)*log1p(nu*exp(-k*(t-t0)))

gen_trait <- function(scn, sdv){                            # -> list(delta[B], res[T], target edges)
  delta <- numeric(B); tg <- integer(0)
  if (scn=="iid") res <- rnorm(T, 0, sdv)
  if (scn=="bm") { d <- rnorm(B); delta <- d * sdv / sd(as.numeric(P %*% d)); res <- rnorm(T, 0, 0.15*sdv) }
  if (scn=="shift") { tg <- sample(which(ndesc>=4 & ndesc<=9), 1)
    delta[tg] <- sample(c(-1,1),1) * 1.6*sdv; res <- rnorm(T, 0, 0.3*sdv) }
  if (scn=="conv") { cand <- which(ndesc>=2 & ndesc<=4)
    repeat { tg <- sample(cand, 2)                          # disjoint clades whose lineages share <= 2 branches
      t1 <- which(P[,tg[1]]==1); t2 <- which(P[,tg[2]]==1)
      if (!length(intersect(t1,t2)) && sum(P[t1[1],]*P[t2[1],]) <= 2) break }
    delta[tg] <- sample(c(-1,1),1) * 1.6*sdv; res <- rnorm(T, 0, 0.3*sdv) }
  list(delta=delta, res=res, tg=tg)
}

gen_data <- function(scn){
  tr <- list(logA=gen_trait("iid", SDV["logA"]), logk=gen_trait(scn, SDV["logk"]),
             lognu=gen_trait(scn, SDV["lognu"]), t0=gen_trait("iid", SDV["t0"]))
  delta <- sapply(tr, `[[`, "delta"); dev <- P %*% delta + sapply(tr, `[[`, "res")   # B x 4, T x 4
  rows <- list(); ti <- list(); ii <- 0
  for (s in seq_len(T)) { nu <- exp(MU["lognu"]+dev[s,3]); k <- exp(MU["logk"]+dev[s,2]); tt <- MU["t0"]+dev[s,4]
    for (r in seq_len(N_PER)) { ii <- ii+1; spec <- paste0(tips[s],"_",r)
      logAi <- MU["logA"]+dev[s,1]+rnorm(1,0,SD_LOGA_IND)
      death <- sample(7:22,1); resorbed <- rbinom(1,1,0.6); offset <- if (resorbed) min(sample(1:7,1), death-3) else 0
      efs <- as.integer(death>=15 & runif(1)<0.6)
      ta <- if (offset==0) c(0,1:death) else (offset+1):death
      sz <- exp(sapply(ta, rich, logA=logAi, nu=nu, k=k, t0=tt) + rnorm(length(ta),0,SIGMA))
      rows[[ii]] <- data.frame(species=tips[s], specimen=spec, age=ta-offset, size=sz, dimension="sim",
                               resorption=ifelse(resorbed,"present","absent"), efs_present=efs==1)
      ti[[ii]] <- data.frame(specimen=spec, offset_true=offset) } }
  list(data=bind_rows(rows), delta=delta, dev=dev, tr=tr, truth_ind=bind_rows(ti))
}

mod <- cmdstan_model(file.path(root,"src/23_tree_nested_model.stan"))
fit_rep <- function(g, scn){
  dat <- g$data; dat$species <- factor(dat$species, levels=tips)
  sd0 <- make_standata(dat, tree, setNames(tips, tips))
  f <- mod$sample(data=sd0$standata, chains=2, parallel_chains=2, iter_warmup=300, iter_sampling=300,
                  adapt_delta=0.9, max_treedepth=10, refresh=0, seed=sample.int(1e6,1), init=0.3,
                  show_messages=FALSE, show_exceptions=FALSE)
  tos <- sd0$tip_of_sp                                       # = 1:T
  delta <- get_arr(f, "delta", B); dev <- get_arr(f, "dev", T)
  pp <- pair_partition(delta, dev, P, tos); ps <- partition_summary(pp, delta, dev, P, tos)
  dT <- array(g$delta, c(1,B,4)); vT <- array(g$dev, c(1,T,4))
  ppT <- pair_partition(dT, vT, P, tos); psT <- partition_summary(ppT, dT, vT, P, tos)
  lam <- apply(as_draws_matrix(f$draws("lambda")), 2, median)
  off <- apply(as_draws_matrix(f$draws("miss_yr")), 2, median)
  ti <- g$truth_ind[match(as.character(sd0$ind_tbl$specimen), g$truth_ind$specimen),]
  sc <- function(x,y) if (sd(x)>0 && sd(y)>0) cor(x,y) else NA_real_
  bind_rows(lapply(1:4, function(k){
    tg <- g$tr[[k]]$tg; pr <- pp$pairs; tgt <- rep(FALSE, nrow(pr))
    if (scn=="shift" && length(tg)) { m <- which(P[,tg]==1); tgt <- pr[,1] %in% m & pr[,2] %in% m }
    if (scn=="conv"  && length(tg)) { m1 <- which(P[,tg[1]]==1); m2 <- which(P[,tg[2]]==1)
      tgt <- (pr[,1] %in% m1 & pr[,2] %in% m2) | (pr[,1] %in% m2 & pr[,2] %in% m1) }
    ppar <- colMeans(pp$independent[,,k] > 0); inh <- colMeans(pp$inherited[,,k]); ind <- colMeans(pp$independent[,,k])
    data.frame(trait=TRAITS[k], lambda=lam[k],
      phy_real=median(ps$phy_real[,k]), phy_real_true=psT$phy_real[1,k],
      anc_share=median(ps$anc_share[,k]), anc_share_true=psT$anc_share[1,k],
      tip_cor=sc(colMeans(dev[,,k]), g$dev[,k]),
      inh_cor=sc(inh, ppT$inherited[1,,k]), ind_cor=sc(ind, ppT$independent[1,,k]),
      n_target=sum(tgt),
      tgt_ppar=if (any(tgt)) mean(ppar[tgt]) else NA, oth_ppar=mean(ppar[!tgt]),
      tgt_inh=if (any(tgt)) mean(inh[tgt]) else NA, tgt_inh_true=if (any(tgt)) mean(ppT$inherited[1,tgt,k]) else NA,
      tgt_ind=if (any(tgt)) mean(ind[tgt]) else NA, tgt_ind_true=if (any(tgt)) mean(ppT$independent[1,tgt,k]) else NA,
      offset_cor=cor(ti$offset_true, off), ndiv=sum(f$diagnostic_summary()$num_divergent))
  }))
}

t0all <- Sys.time()
for (rep in seq_len(nreps)) for (scn in regimes) {
  r <- tryCatch(fit_rep(gen_data(scn), scn), error=function(e){ cat("rep err", scn, rep, conditionMessage(e), "\n"); NULL })
  if (is.null(r)) next
  r$regime <- scn; r$rep <- rep; r$chunk <- chunk; r$n_per <- N_PER
  write_csv(r, resfile, append=file.exists(resfile))
  cat(sprintf("[%s] rep %d/%d | %.1f min\n", scn, rep, nreps, as.numeric(Sys.time()-t0all, units="mins")))
}
cat("DONE tree recovery\n")
