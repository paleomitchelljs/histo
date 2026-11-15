# Stan data for model 36 (src/36_age_element_model.stan). Requires src/23_tree_helpers.R, 32_helpers.R and
# 35_helpers.R to be sourced first. Adds to make_standata35:
#  (a) species x element-type offset indices (r_g, c_g, NG): for a species with sized elements of two or more
#      long-bone types, every type except its reference type (femur if present, else its best-sampled type)
#      gets an offset; relative channels (pct_*) get none;
#  (c) marks preserved (r_last) counted from mark rows only: a ring row with zone "periosteal" is a size at death
#      and does not count as a mark. Lmax and t_loc are recomputed accordingly.
LONG_BONES <- c("femur", "humerus", "ulna", "radius", "tibia", "fibula")

make_standata36 <- function(rings, els, inds, tree, tip_of_species, ...){
  if (!"zone" %in% names(rings)) rings$zone <- "mark"
  rings$zone[is.na(rings$zone)] <- "mark"
  x <- make_standata35(rings, els, inds, tree, tip_of_species, ...)
  d <- x$standata; er <- x$er; ec <- x$ec; inds <- x$inds
  # ---- (c) marks preserved from mark rows only ----
  mk <- x$rings[is.na(x$rings$zone) | x$rings$zone != "periosteal", ]
  lm <- tapply(mk$age0, mk$el_id, max)
  er$last <- pmax(as.numeric(lm[er$el_id]), er$count, na.rm=TRUE)
  er$last[is.na(er$last)] <- 0
  cnt <- c(tapply(er$last, er$individual, max), tapply(ec$count[ec$has_count], ec$individual[ec$has_count], max))
  inds$Lmax <- vapply(inds$individual, function(i){ v <- cnt[names(cnt) == i]; if (length(v)) max(v) else 0 }, 0)
  inds$t_loc <- ifelse(inds$age_type == 1, pmax(inds$age_est - inds$Lmax, 0), 0)
  d$r_last <- er$last; d$Lmax <- inds$Lmax; d$t_loc <- inds$t_loc
  # ---- (a) species x element-type offsets ----
  sized <- rbind(data.frame(species=er$species, eltype=er$eltype, dimension=er$dimension),
                 data.frame(species=ec$species[ec$has_size], eltype=ec$eltype[ec$has_size], dimension=ec$dimension[ec$has_size]))
  sized <- sized[sized$eltype %in% LONG_BONES & !sized$dimension %in% REL_DIMS, ]
  keys <- character(0); refs <- list()
  for (sp in unique(sized$species)) {
    tb <- table(sized$eltype[sized$species == sp]); if (length(tb) < 2) next
    ref <- if ("femur" %in% names(tb)) "femur" else names(tb)[which.max(tb)]
    refs[[sp]] <- ref; keys <- c(keys, paste(sp, setdiff(names(tb), ref), sep=" | ")) }
  gkey <- function(sp, et, dim) ifelse(et %in% LONG_BONES & !dim %in% REL_DIMS, paste(sp, et, sep=" | "), NA)
  ix0 <- function(v, lev) { m <- match(v, lev); m[is.na(m)] <- 0L; as.integer(m) }
  d$NG <- length(keys)
  d$r_g <- as.array(ix0(gkey(er$species, er$eltype, er$dimension), keys))
  d$c_g <- as.array(ix0(gkey(ec$species, ec$eltype, ec$dimension), keys))
  x$standata <- d; x$er <- er; x$inds <- inds; x$gamma_keys <- keys; x$gamma_refs <- unlist(refs)
  x
}
make_init36 <- function(x) { f35 <- make_init35(x); function(){ d <- x$standata; v <- f35(); j <- function(n, s=0.1) rnorm(n, 0, s)
  v$theta_o <- v$theta_raw[, c(1, 2, 4)]; v$z_nu <- as.array(j(d$S, 0.3)); v$theta_raw <- NULL
  v$sd_gam <- 0.2; v$z_gam <- as.array(j(d$NG))
  v } }
