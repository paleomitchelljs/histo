# Stan data for model 37 (src/37_age_element_model.stan). Requires src/23_tree_helpers.R, 32_helpers.R,
# 35_helpers.R and 36_helpers.R to be sourced first. Adds to make_standata36:
#  (d) s_open: scale (marks) of the unknown extra preserved marks behind an "at least lo" count;
#  (e) r_bi / c_bi: per-element exponent index. An element uses its channel's free exponent only if its species has
#      sized data in the reference (fixed) channel of that channel's co-occurrence component; otherwise 0 (beta = 1).
make_standata37 <- function(rings, els, inds, tree, tip_of_species, s_open=2, ...){
  x <- make_standata36(rings, els, inds, tree, tip_of_species, ...)
  d <- x$standata; ch <- x$channels; er <- x$er; ec <- x$ec
  ref_of_comp <- setNames(ch$channel[ch$beta == "fixed 1"], ch$component[ch$beta == "fixed 1"])
  sized <- unique(rbind(data.frame(species=er$species, channel=er$channel),
                        data.frame(species=ec$species[ec$has_size], channel=ec$channel[ec$has_size])))
  has_ref <- function(sp, chan) { r <- ref_of_comp[as.character(ch$component[match(chan, ch$channel)])]
    mapply(function(s, rr) any(sized$species == s & sized$channel == rr), sp, r) }
  bi <- function(sp, chan) { b <- d$beta_idx[match(chan, ch$channel)]; as.array(as.integer(ifelse(b > 0 & has_ref(sp, chan), b, 0L))) }
  d$r_bi <- bi(er$species, er$channel); d$c_bi <- bi(ec$species, ec$channel)
  d$s_open <- s_open
  # species that inform each free exponent (sized data in the channel AND in its reference channel)
  inf <- rbind(data.frame(species=er$species, channel=er$channel, used=d$r_bi > 0),
               data.frame(species=ec$species[ec$has_size], channel=ec$channel[ec$has_size], used=d$c_bi[ec$has_size] > 0))
  ch$n_species_informing <- vapply(ch$channel, function(c) length(unique(inf$species[inf$channel == c & inf$used])), 0L)
  ch$species_informing <- vapply(ch$channel, function(c) paste(sort(unique(inf$species[inf$channel == c & inf$used])), collapse=", "), "")
  ch$n_sized_informing <- vapply(ch$channel, function(c) sum(inf$channel == c & inf$used), 0L)
  x$standata <- d; x$channels <- ch
  x
}
make_init37 <- function(x) make_init36(x)
