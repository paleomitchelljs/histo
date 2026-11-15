# Helpers for the tree-nested model (23): tree matrices, Stan data, and the
# ancestry-vs-independent partition computed from posterior draws.
TRAITS <- c("logA","logk","lognu","t0")

# root->tip path matrix P (T x B), plus edge table
tree_mats <- function(tree){
  edge <- tree$edge; B <- nrow(edge); T <- length(tree$tip.label)
  P <- matrix(0, T, B)
  for (i in seq_len(T)) { node <- i
    repeat { e <- which(edge[,2]==node); if (!length(e)) break; P[i,e] <- 1; node <- edge[e,1] } }
  rownames(P) <- tree$tip.label
  list(P=P, edge=edge, B=B, T=T, tips=tree$tip.label)
}

# draws -> arrays [draw, B, K] and [draw, S, K]
get_arr <- function(fit, var, n1, K=4){
  m <- posterior::as_draws_matrix(fit$draws(var))
  array(as.numeric(m), c(nrow(m), n1, K))
}

# Partition for every species pair (i,j), per trait, per draw.
#   x_i = species deviation from the root state
#   a   = summed change on the branches i and j SHARE        (common ancestry)
#   d_i = x_i - a = change private to i's lineage (private branches + species residual)
#   x_i * x_j = a^2 [inherited] + d_i*d_j [independent] + a*(d_i+d_j) [cross]
# inherited >= 0 always; independent > 0 = parallel change on separate lineages,
# < 0 = divergence. delta: [n,B,K], dev: [n,S,K].
pair_partition <- function(delta, dev, P, tip_of_sp){
  n <- dim(delta)[1]; S <- dim(dev)[2]; K <- dim(delta)[3]
  pr <- t(combn(S, 2)); np <- nrow(pr)
  inh <- ind <- crs <- array(0, c(n, np, K))
  for (p in seq_len(np)) {
    i <- pr[p,1]; j <- pr[p,2]
    m <- P[tip_of_sp[i],] * P[tip_of_sp[j],]
    for (k in seq_len(K)) {
      a  <- as.numeric(matrix(delta[,,k], n) %*% m)
      di <- dev[,i,k] - a; dj <- dev[,j,k] - a
      inh[,p,k] <- a^2; ind[,p,k] <- di*dj; crs[,p,k] <- a*(di+dj)
    }
  }
  list(pairs=pr, inherited=inh, independent=ind, cross=crs)
}

# Per-draw summaries per trait:
#  anc_share  = sum(inherited) / (sum(inherited) + sum(|independent|))  over pairs
#  phy_real   = realized share of among-species variance carried by the tree
#               component: cov(phy, dev) / var(dev)   (finite-population)
partition_summary <- function(pp, delta, dev, P, tip_of_sp){
  n <- dim(delta)[1]; K <- dim(delta)[3]
  anc <- phy_real <- matrix(NA, n, K)
  for (k in seq_len(K)) {
    anc[,k] <- rowSums(matrix(pp$inherited[,,k], n)) /
      (rowSums(matrix(pp$inherited[,,k], n)) + rowSums(abs(matrix(pp$independent[,,k], n))))
    ph <- matrix(delta[,,k], n) %*% t(P[tip_of_sp,,drop=FALSE])       # n x S
    dv <- matrix(dev[,,k], n)
    phc <- ph - rowMeans(ph); dvc <- dv - rowMeans(dv)
    phy_real[,k] <- rowSums(phc*dvc) / rowSums(dvc^2)
  }
  colnames(anc) <- colnames(phy_real) <- TRAITS[seq_len(K)]
  list(anc_share=anc, phy_real=phy_real)
}

# Stan data from a long LAG table (columns species, specimen, size, age,
# dimension, resorption, efs_present) + tree + species->tip map.
make_standata <- function(dat, tree, tip_of_species, est_kappa=0L, brlen=NULL){
  dat$species <- factor(dat$species); dat$specimen <- factor(dat$specimen)
  dat$dimr <- as.integer(factor(dat$dimension))
  it <- dat |> dplyr::distinct(specimen, species, resorption, efs_present, dimr) |>
    dplyr::left_join(dat |> dplyr::group_by(specimen) |>
      dplyr::summarise(logmax=log(max(size)), .groups="drop"), by="specimen") |> dplyr::arrange(specimen)
  stopifnot(nrow(it) == nlevels(dat$specimen))
  tm <- tree_mats(tree); spl <- levels(dat$species)
  tip_of_sp <- match(tip_of_species[spl], tm$tips); stopifnot(!anyNA(tip_of_sp))
  if (is.null(brlen)) brlen <- rep(1, tm$B)
  list(standata=list(N=nrow(dat), S=length(spl), I=nrow(it), D=max(dat$dimr), T=tm$T, B=tm$B,
      size=dat$size, age0=dat$age, ind=as.integer(dat$specimen), dimr=dat$dimr,
      dim_of_ind=it$dimr, sp_of_ind=as.integer(it$species), tip_of_sp=tip_of_sp, P=tm$P,
      brlen=brlen, est_kappa=est_kappa, efs=as.integer(it$efs_present),
      resorb=as.integer(it$resorption!="absent"), logmax=it$logmax),
    species=spl, ind_tbl=it, tm=tm, tip_of_sp=tip_of_sp)
}
