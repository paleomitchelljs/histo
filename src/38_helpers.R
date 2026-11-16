# Stan data for model 38 (src/38_age_element_model.stan). Requires src/23_tree_helpers.R, 32_helpers.R, 35_helpers.R,
# 36_helpers.R and 37_helpers.R to be sourced first. Adds to make_standata37: cl_of_sp (clade of each species, the M7
# clades) for the clade-pooled nu (g). Elements and counts are unchanged; the inner/outer split (f) is model-side.
TREE_TRAITS <- c("logA", "logr", "t50")                 # traits on the tree in model 38 (nu is pooled by clade)
make_standata38 <- function(rings, els, inds, tree, tip_of_species, ...){
  x <- make_standata37(rings, els, inds, tree, tip_of_species, ...)
  d <- x$standata
  d$cl_of_sp <- as.array(match(x$clade_sp[x$species], x$clades))
  stopifnot(!anyNA(d$cl_of_sp))
  x$standata <- d
  x
}
make_init38 <- function(x) { f37 <- make_init37(x); function(){ d <- x$standata; v <- f37(); j <- function(n, s=0.1) rnorm(n, 0, s)
  v$sd_phy <- v$sd_phy[c(1, 2, 4)]; v$sd_res <- v$sd_res[c(1, 2, 4)]
  v$z_nu <- NULL; v$kappa_raw <- NULL; v$f_ero <- NULL
  v$sd_nu_cl <- 0.2; v$sd_nu_sp <- 0.1; v$z_nu_cl <- as.array(j(d$NTC)); v$z_nu_sp <- as.array(j(d$S))
  v$f_el <- as.array(rep(0.5, d$ER)); v$tau_out <- 0.5
  v } }
