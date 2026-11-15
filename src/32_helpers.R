# Stan data for model 32 from the broad tables written by src/31_stage_broad.R.
# Uses tree_mats() from src/23_tree_helpers.R.
ELTYPES <- c("humerus","ulna","radius","tibia","fibula")     # femur and body measures are the reference (0)
make_standata32 <- function(rings, els, inds, tree, tip_of_species, est_kappa=0L, brlen=NULL,
                            efs_thr=0.90, efs_slope=0.04, tau_fresh=0.25){
  spl <- sort(unique(inds$species)); inds <- inds[order(inds$species, inds$individual),]
  dims <- sort(unique(els$dimension))
  er <- els[els$kind=="ring",]; ec <- els[els$kind=="cs",]
  er <- er[order(match(er$individual, inds$individual), er$el_id),]; ec <- ec[order(match(ec$individual, inds$individual), ec$el_id),]
  rings <- rings[order(match(rings$el_id, er$el_id), rings$age0),]
  rs <- do.call(rbind, lapply(split(rings, rings$el_id), function(d) data.frame(el_id=d$el_id[1], last=max(d$age0), logmax=log(max(d$size)))))
  er <- merge(er, rs, by="el_id", sort=FALSE); er <- er[order(match(er$individual, inds$individual), er$el_id),]
  er$last <- pmax(er$last, er$count, na.rm=TRUE)
  ec$has_count <- !is.na(ec$count); ec$has_size <- !is.na(ec$size)
  # one EFS term per individual among its count-and-size elements (they share one point on the curve)
  for (i in unique(ec$individual)) { j <- which(ec$individual==i); st <- ec$efs[j]
    keep <- if (any(st==1)) j[st==1][1] else if (any(st==0)) j[st==0][1] else NA
    ec$efs[setdiff(j, keep)] <- 2L }
  # individual lower bound for the age, softness, multi-element index
  cnt <- c(tapply(er$last, er$individual, max), tapply(ec$count[ec$has_count], ec$individual[ec$has_count], max))
  inds$Lmax <- vapply(inds$individual, function(i){ v <- cnt[names(cnt)==i]; if (length(v)) max(v) else 0 }, 0)
  res_any <- c(tapply(er$resorb, er$individual, max), tapply(ec$resorb | ec$count_min, ec$individual, max))
  inds$csoft <- ifelse(vapply(inds$individual, function(i) any(res_any[names(res_any)==i] > 0), TRUE), 4, 16)
  nsized <- table(c(er$individual, ec$individual[ec$has_size]))
  multi <- inds$individual %in% names(nsized)[nsized > 1]
  inds$im <- ifelse(multi, cumsum(multi), 0L)
  inds$age_type <- as.integer(!is.na(inds$age_est))
  # age parametrization: 0 = tied to t50 (count-based animals with a resorbed ring series), 1 = direct
  # Only RESORBED ring series are tied to t50. Tying the non-resorbed ones as well was tried
  # (fit `core2`, 2026-11-14) and was far worse: R-hat up to 3.8 against 1.05.
  ring_res <- inds$individual %in% er$individual[er$resorb==1]
  inds$tpar <- as.integer(!(ring_res & inds$age_type==0))
  inds$t_loc <- ifelse(inds$age_type==1, pmax(inds$age_est - inds$Lmax, 0), 0)
  inds$t_scale <- ifelse(inds$age_type==1, inds$age_sd, 1)
  rdims <- sort(unique(er$dimension))
  er$ero <- ifelse(er$eroded==1, cumsum(er$eroded==1), 0L)
  tm <- tree_mats(tree); tip_of_sp <- match(tip_of_species[spl], tm$tips); stopifnot(!anyNA(tip_of_sp))
  if (is.null(brlen)) brlen <- rep(1, tm$B)
  ob <- rbind(data.frame(species=rings$species, dimension=er$dimension[match(rings$el_id, er$el_id)]),
              data.frame(species=ec$species[ec$has_size], dimension=ec$dimension[ec$has_size]))
  main_dim <- vapply(spl, function(sp){ t <- table(ob$dimension[ob$species==sp]); match(names(t)[which.max(t)], dims) }, 1L)
  ix <- function(x, lev) { m <- match(x, lev); m[is.na(m)] <- 0L; as.integer(m) }
  sd <- list(S=length(spl), D=length(dims), T=tm$T, B=tm$B, L=length(ELTYPES), tip_of_sp=tip_of_sp, P=tm$P, brlen=brlen, est_kappa=est_kappa,
    I=nrow(inds), sp_of_ind=match(inds$species, spl), age_type=inds$age_type, Lmax=inds$Lmax,
    age_est=ifelse(is.na(inds$age_est), 0, inds$age_est), age_sd=ifelse(is.na(inds$age_sd), 1, inds$age_sd), csoft=inds$csoft,
    IM=sum(multi), im_of_ind=as.integer(inds$im), tpar=inds$tpar, t_loc=inds$t_loc, t_scale=inds$t_scale,
    DR=length(rdims), r_sdim=match(er$dimension, rdims),
    ER=nrow(er), r_ind=match(er$individual, inds$individual), r_dim=match(er$dimension, dims), r_elt=ix(er$eltype, ELTYPES),
    r_efs=as.integer(er$efs), r_resorb=as.integer(er$resorb), NERO=sum(er$eroded==1), r_ero=as.integer(er$ero),
    r_last=er$last, r_logmax=er$logmax,
    N=nrow(rings), size=rings$size, age0=rings$age0, el=match(rings$el_id, er$el_id),
    EC=nrow(ec), c_ind=match(ec$individual, inds$individual), c_dim=match(ec$dimension, dims), c_elt=ix(ec$eltype, ELTYPES),
    c_efs=as.integer(ec$efs), c_resorb=as.integer(ec$resorb | ec$count_min), c_has_count=as.integer(ec$has_count),
    c_count=ifelse(ec$has_count, ec$count, 0), c_has_size=as.integer(ec$has_size), c_logsize=ifelse(ec$has_size, log(ec$size), 0),
    c_logsd=ec$size_logsd, c_nrep=ec$nrep, c_tadd=ec$t_add, main_dim=as.array(unname(main_dim)), efs_thr=efs_thr, efs_slope=efs_slope, tau_fresh=tau_fresh)
  for (n in c("r_ind","r_dim","r_elt","r_efs","r_resorb","r_ero","c_ind","c_dim","c_elt","c_efs","c_resorb","c_has_count","c_has_size","el","sp_of_ind","age_type","im_of_ind","tip_of_sp","tpar","r_sdim"))
    sd[[n]] <- as.array(sd[[n]])
  list(standata=sd, species=spl, dims=dims, rdims=rdims, inds=inds, er=er, ec=ec, rings=rings, tm=tm, tip_of_sp=tip_of_sp)
}
# starting values near the data (the default random inits put asymptotes at ~1)
make_init32 <- function(x) function(){ d <- x$standata; j <- function(n, s=0.1) rnorm(n, 0, s)
  mla <- vapply(seq_len(d$D), function(k){ v <- c(d$r_logmax[d$r_dim==k], d$c_logsize[d$c_dim==k & d$c_has_size==1]); if (length(v)) max(v) else 4 }, 0)
  t50 <- 5 + j(d$S, 0.3); want <- ifelse(d$age_type==1, pmax(d$age_est - d$Lmax, 0.3), 0.5)
  list(mu_logA=mla + j(d$D), mu_logr=-2.5 + j(1), mu_nu=j(1), mu_t50=5 + j(1), sd_phy=rep(0.3, 4)*exp(j(4)), sd_res=rep(0.3, 4)*exp(j(4)),
       theta_raw=cbind(mla[d$main_dim] + j(d$S), -2.5 + j(d$S), j(d$S), t50), beta_el=as.array(j(d$L)), sd_logA_ind=0.3*exp(j(1)), rho_ind=0.5,
       z_ind=as.array(j(d$IM)), logA_e=as.array(d$r_logmax + 0.3 + j(d$ER)), kappa_raw=0.5, tau_off=2*exp(j(1)),
       t_raw=as.array(ifelse(d$tpar==0, t50[d$sp_of_ind] - want, (want - d$t_loc)/d$t_scale) + j(d$I)),
       f_ero=as.array(rep(0.5, d$NERO)), sigma=as.array(0.1*exp(j(d$DR))), sigma_cs=0.1*exp(j(1))) }
