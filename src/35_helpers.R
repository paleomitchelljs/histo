# Stan data for model 35 (src/35_revision_model.stan) from the tables written by src/31_stage_broad.R.
# Builds on make_standata32 (src/32_helpers.R) and adds: size channels and which exponents are free
# (M1), relative-channel intercepts, explicit vs integrated individual timing effects (M2),
# interval counts (M5), sex (M6), and lost-mark scale groups by element type and clade (M7).
# Requires src/23_tree_helpers.R and src/32_helpers.R to be sourced first.
REF_CHANNEL <- "femur|circumference_mm"
REL_DIMS <- c("pct_max_size", "pct_bone_length")          # scaled to the largest specimen of a sample
BODY_DIMS <- c("svl_mm", "total_length_mm", "bsl_mm")
TAU_ELTYPES <- c("femur","humerus","ulna","radius","tibia","fibula","other")
CLADE_TIPS <- list(                                        # tips by major clade, for the lost-mark scale (M7)
  Lissamphibia = c("Andrias","Oreophrynella","Pristimantis"),
  Lepidosauria = c("Ophisops","Phymaturus","Sauromalus","Phrynocephalus"),
  Crocodylomorpha = c("Orthosuchus","Crocodylus","Diplocynodon","Alligator","Caiman"),
  Dinosauriformes = c("Asilisaurus","Jeholosaurus","Orodromeus","Psittacosaurus","Tenontosaurus","Dryosaurus","Maiasaura","Hypacrosaurus",
    "Plateosaurus","Massospondylus","Mamenchisaurus","Northampton","Janenschia","Apatosaurus","Diplodocus","Camarasaurus","Alamosaurus",
    "Syntarsus","Coelophysis","Ceratosaurus","Allosaurus","Acrocanthosaurus","Meraxes","Albertosaurus","Gorgosaurus","Tyrannosaurus",
    "Rativates","Saurornitholestes"))                       # every other tip: Synapsida
channel_of <- function(eltype, dimension) ifelse(dimension %in% BODY_DIMS | eltype == "other", dimension, paste(eltype, dimension, sep="|"))
clade_of_tip <- function(tip) { cl <- rep("Synapsida", length(tip))
  for (n in names(CLADE_TIPS)) cl[tip %in% CLADE_TIPS[[n]]] <- n; cl }

make_standata35 <- function(rings, els, inds, tree, tip_of_species, ...){
  # interval counts: the age floor (Lmax) and the 32 builder use the LOWER bound of each count
  els$count_mid <- els$count
  int <- !is.na(els$count) & (els$count_lo != els$count | els$count_hi != els$count)
  els$count[int] <- els$count_lo[int]
  x <- make_standata32(rings, els, inds, tree, tip_of_species, ...)
  d <- x$standata; er <- x$er; ec <- x$ec; inds <- x$inds; spl <- x$species
  # ---- channels (M1) ----
  er$channel <- channel_of(er$eltype, er$dimension); ec$channel <- channel_of(ec$eltype, ec$dimension)
  obs <- rbind(data.frame(species=er$species, channel=er$channel, n=1), data.frame(species=ec$species[ec$has_size], channel=ec$channel[ec$has_size], n=1))
  chans <- sort(unique(c(er$channel, ec$channel)))
  # co-occurrence components (channels observed in a common species)
  comp <- setNames(seq_along(chans), chans)
  for (sp in unique(obs$species)) { cs <- unique(obs$channel[obs$species == sp]); if (length(cs) > 1) { k <- comp[cs]; comp[comp %in% k] <- min(k) } }
  nobs <- table(factor(obs$channel, chans))
  fixed <- rep(FALSE, length(chans)); names(fixed) <- chans
  for (k in unique(comp)) { m <- chans[comp == k]
    if (REF_CHANNEL %in% m) fixed[REF_CHANNEL] <- TRUE else fixed[m[which.max(nobs[m])]] <- TRUE }
  beta_idx <- ifelse(fixed, 0L, cumsum(!fixed))
  chtab <- data.frame(channel=chans, component=match(comp, unique(comp)), n_sized=as.integer(nobs), beta=ifelse(fixed, "fixed 1", "free"),
                      n_species=vapply(chans, function(c) length(unique(obs$species[obs$channel == c])), 0L), row.names=NULL)
  # relative channels: one free intercept per species x channel
  rel_key <- function(sp, ch, dim) ifelse(dim %in% REL_DIMS, paste(sp, ch), NA)
  rk <- unique(na.omit(c(rel_key(er$species, er$channel, er$dimension), rel_key(ec$species, ec$channel, ec$dimension))))
  ix0 <- function(v, lev) { m <- match(v, lev); m[is.na(m)] <- 0L; as.integer(m) }
  # main dimension per species: best-sampled NON-relative dimension (fallback: circumference)
  ob2 <- rbind(data.frame(species=rings$species, dimension=er$dimension[match(rings$el_id, er$el_id)]),
               data.frame(species=ec$species[ec$has_size], dimension=ec$dimension[ec$has_size]))
  ob2 <- ob2[!ob2$dimension %in% REL_DIMS, ]
  main_dim <- vapply(spl, function(sp){ t <- table(ob2$dimension[ob2$species == sp])
    if (length(t)) match(names(t)[which.max(t)], x$dims) else match("circumference_mm", x$dims) }, 1L)
  # ---- individual timing (M2): explicit where an individual has a ring series or > 1 sized element ----
  nsz <- table(c(er$individual, ec$individual[ec$has_size]))
  expl <- inds$individual %in% er$individual | inds$individual %in% names(nsz)[nsz > 1]
  # ---- sex (M6) ----
  sexc <- ifelse(is.na(inds$sex), 0, ifelse(inds$sex == "male", 0.5, -0.5))
  both <- names(which(tapply(inds$sex, inds$species, function(s) all(c("female","male") %in% s))))
  # ---- lost-mark groups (M7) ----
  tip_names <- x$tm$tips[x$tip_of_sp]; clade_sp <- setNames(clade_of_tip(tip_names), spl)
  clades <- sort(unique(clade_sp))
  # ---- count types (M5) ----
  ec$cnt_type <- ifelse(!ec$has_count, 0L, ifelse(ec$count_lo == ec$count_hi, 1L, 2L))
  d <- c(d, list(
    IU=sum(expl), u_of_ind=as.array(ifelse(expl, cumsum(expl), 0L)),
    sexc=as.array(sexc), NSX=length(both), sx_of_sp=as.array(ix0(spl, both)),
    C=length(chans), NB=sum(!fixed), beta_idx=as.array(unname(beta_idx)), NO=length(rk),
    NTE=length(TAU_ELTYPES), NTC=length(clades),
    r_ch=as.array(match(er$channel, chans)), r_relch=as.array(ix0(rel_key(er$species, er$channel, er$dimension), rk)),
    r_te=as.array(match(er$eltype, TAU_ELTYPES)), r_tc=as.array(match(clade_sp[er$species], clades)),
    c_ch=as.array(match(ec$channel, chans)), c_relch=as.array(ix0(rel_key(ec$species, ec$channel, ec$dimension), rk)),
    c_te=as.array(match(ec$eltype, TAU_ELTYPES)), c_tc=as.array(match(clade_sp[ec$species], clades)),
    c_cnt_type=as.array(ec$cnt_type), c_lo=ifelse(ec$has_count, ec$count_lo, 0),
    c_hi=ifelse(ec$has_count & is.finite(ec$count_hi), ec$count_hi, 0), c_open=as.array(as.integer(ec$has_count & !is.finite(ec$count_hi)))))
  d$c_count <- ifelse(ec$has_count, ec$count_mid, 0)     # exact counts: the printed value
  d$main_dim <- as.array(unname(main_dim))
  d$c_has_count <- NULL
  stopifnot(!anyNA(d$r_ch), !anyNA(d$c_ch), !anyNA(d$r_te), !anyNA(d$c_te))
  x$standata <- d; x$er <- er; x$ec <- ec; x$channels <- chtab; x$rel_keys <- rk; x$clades <- clades; x$clade_sp <- clade_sp
  x$sex_species <- both; x$explicit_u <- expl
  x
}
make_init35 <- function(x) { f32 <- make_init32(x); function(){ d <- x$standata; v <- f32(); j <- function(n, s=0.1) rnorm(n, 0, s)
  v$beta_el <- as.array(j(d$L))
  v$log_tau0 <- log(2) + j(1); v$sd_tel <- 0.2; v$sd_tcl <- 0.2; v$z_tel <- as.array(j(d$NTE)); v$z_tcl <- as.array(j(d$NTC)); v$tau_off <- NULL
  v$log_beta <- as.array(j(d$NB, 0.05)); v$o_rel <- as.array(log(100) + j(d$NO))
  v$sd_u <- 0.3; v$z_u <- as.array(j(d$IU)); v$gbar_sex <- 0; v$sd_sex <- 0.05; v$z_sex <- as.array(j(d$NSX))
  v } }
