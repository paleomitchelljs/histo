#!/usr/bin/env Rscript
# =====================================================================
# 06_phylo_growth_figure.R
# Composite figure: phylogeny (branches colored by the most-variable
# parameter = shape nu) + silhouette at each tip + per-taxon posterior
# growth-curve ribbon (normalized shape) to the right.
#
# Silhouettes: tries rphylopic (needs internet). In a network-blocked
# environment it falls back to built-in clade glyphs so the figure still
# renders complete; run with internet to get real PhyloPic art.
# =====================================================================
suppressMessages({
  library(ape); library(ggplot2); library(dplyr); library(readr)
  library(viridisLite); library(grid); library(brms)
})
root <- "/Users/jmitchell/Library/CloudStorage/Dropbox/Research/00_current_Histology"
out  <- file.path(root, "outputs")
have_net <- requireNamespace("rphylopic", quietly = TRUE) && {
  old <- options(timeout = 4); on.exit(options(old), add = TRUE)
  !inherits(try(suppressWarnings(readLines(url("https://api.phylopic.org"), n = 1, warn = FALSE)),
                silent = TRUE), "try-error")
}

# ---- tree (8 tips; topology from clade columns) --------------------
nwk <- paste0("(Maiasaura_peeblesorum,((Allosaurus_fragilis,",
              "(Acrocanthosaurus_atokensis,Campanas_Carcharodontosaurid)),",
              "(Rativates_evadens,(Gorgosaurus_libratus,",
              "(Tyrannosaurus_rex,cf.Tyrannosaurus_rex)))));")
tree <- read.tree(text = nwk)
treeU <- compute.brlen(tree, method = "Grafen")     # ultrametric -> aligned tips

disp <- c(Tyrannosaurus_rex = "Tyrannosaurus rex", cf.Tyrannosaurus_rex = "cf. Tyrannosaurus (Jane)",
          Gorgosaurus_libratus = "Gorgosaurus", Rativates_evadens = "Rativates",
          Campanas_Carcharodontosaurid = "Meraxes gigas", Acrocanthosaurus_atokensis = "Acrocanthosaurus",
          Allosaurus_fragilis = "Allosaurus", Maiasaura_peeblesorum = "Maiasaura")
clade <- c(Tyrannosaurus_rex = "theropod", cf.Tyrannosaurus_rex = "theropod",
           Gorgosaurus_libratus = "theropod", Rativates_evadens = "ornithomimid",
           Campanas_Carcharodontosaurid = "theropod", Acrocanthosaurus_atokensis = "theropod",
           Allosaurus_fragilis = "theropod", Maiasaura_peeblesorum = "hadrosaur")

# ---- posterior per-species shape, and normalized growth ribbons ----
fit <- readRDS(file.path(out, "fit_unified_nu_by_species.rds"))
dat <- read_csv(file.path(root, "data_clean/lag_circumference_long.csv"), show_col_types = FALSE)
d <- as_draws_df(fit)
getp <- function(par, s) { c0 <- paste0("r_species__", par, "[", s, ",Intercept]")
  b <- paste0("b_", par, "_Intercept"); if (c0 %in% names(d)) d[[b]] + d[[c0]] else d[[b]] }
lognu_tip <- sapply(tree$tip.label, function(s) median(getp("lognu", s)))

# representative (best-sampled) specimen per species for the fitted curve
rep_spec <- dat %>% count(species, specimen) %>% group_by(species) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>% ungroup()
curve_df <- do.call(rbind, lapply(tree$tip.label, function(s) {
  sp_id <- rep_spec$specimen[rep_spec$species == s]
  amax <- max(dat$age[dat$species == s])
  ag <- seq(0, amax, length.out = 60)
  e <- posterior_epred(fit, newdata = data.frame(age = c(ag, 1e4), species = s, specimen = sp_id),
                       re_formula = NULL)
  A <- e[, ncol(e)]; rel <- e[, seq_along(ag)] / A          # normalize by asymptote
  data.frame(species = s, rel_age = ag / amax,
             med = apply(rel, 2, median), lo = apply(rel, 2, quantile, .05),
             hi = apply(rel, 2, quantile, .95))
}))

# ---- ancestral reconstruction of lognu for branch coloring ---------
av <- lognu_tip[treeU$tip.label]
ace_fit <- ape::ace(av, treeU, type = "continuous", method = "REML")
nodeval <- c(av, ace_fit$ace)                                # index = node number

# ---- tree layout coordinates (ape, headless) -----------------------
pdf(NULL); plot(treeU, plot = FALSE); pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv); invisible(dev.off())
xx <- pp$xx; yy <- pp$yy; ntip <- length(tree$tip.label); span <- max(xx[1:ntip]); tipx <- span
edge <- treeU$edge
seg_h <- data.frame(x = xx[edge[, 1]], xend = xx[edge[, 2]], y = yy[edge[, 2]], yend = yy[edge[, 2]],
                    val = (nodeval[edge[, 1]] + nodeval[edge[, 2]]) / 2)
# vertical connectors at each internal node
seg_v <- do.call(rbind, lapply(unique(edge[, 1]), function(p) {
  ch <- edge[edge[, 1] == p, 2]
  data.frame(x = xx[p], xend = xx[p], y = min(yy[ch]), yend = max(yy[ch]))
}))

# ---- column geometry to the right of the tips ----------------------
# layout: tip -> silhouette -> taxon name -> growth-curve panel
x_sil <- tipx + 0.17 * span; sw <- 0.30 * span   # silhouette centre + width
x_lab <- tipx + 0.36 * span                       # taxon name
xc0   <- tipx + 0.92 * span; W <- 0.52 * span     # growth-curve panel
H <- 0.78                                          # curve vertical extent ~ row
xlim_hi <- xc0 + W + 0.08 * span

tipdf <- data.frame(species = tree$tip.label, x = tipx, y = yy[1:ntip],
                    lab = disp[tree$tip.label], clade = clade[tree$tip.label],
                    val = lognu_tip[tree$tip.label])

# place curves into tip rows
cdf <- curve_df %>% left_join(tipdf[, c("species", "y", "val")], by = "species") %>%
  mutate(px = xc0 + rel_age * W,
         pmed = y + (med - .5) * H, plo = y + (lo - .5) * H, phi = y + (hi - .5) * H)
# light backing panel + baseline/top guides for each curve
panel <- data.frame(xmin = xc0, xmax = xc0 + W, ymin = tipdf$y - H/2, ymax = tipdf$y + H/2)

# ---- built-in fallback silhouettes (normalized polygons) -----------
theropod <- data.frame(
  x = c(.00,.08,.15,.20,.27,.45,.58,.72,1.00,.72,.62,.60,.63,.69,.60,.56,.50,.42,.40,.45,.37,.35,.28,.18,.08),
  y = c(.52,.60,.66,.62,.66,.74,.76,.70,.50,.56,.52,.30,.12,.08,.10,.30,.48,.46,.16,.10,.12,.34,.50,.50,.50))
hadrosaur <- data.frame(
  x = c(.00,.06,.13,.18,.25,.42,.60,.76,1.00,.76,.64,.62,.67,.58,.56,.44,.34,.32,.37,.28,.26,.20,.10),
  y = c(.46,.50,.56,.54,.58,.68,.66,.60,.46,.50,.44,.16,.08,.10,.40,.42,.40,.16,.10,.12,.40,.44,.44))
glyph <- list(theropod = theropod, hadrosaur = hadrosaur, ornithomimid = theropod)
# resolve real PhyloPic silhouettes when online; else fall back to clade glyphs
phylopic_name <- c(Tyrannosaurus_rex = "Tyrannosaurus", cf.Tyrannosaurus_rex = "Tyrannosaurus",
  Gorgosaurus_libratus = "Gorgosaurus", Rativates_evadens = "Ornithomimidae",
  Campanas_Carcharodontosaurid = "Giganotosaurus", Acrocanthosaurus_atokensis = "Acrocanthosaurus",
  Allosaurus_fragilis = "Allosaurus", Maiasaura_peeblesorum = "Maiasaura")
tipdf$uuid <- NA_character_
if (have_net) for (i in seq_len(nrow(tipdf)))
  tipdf$uuid[i] <- tryCatch(rphylopic::get_uuid(name = phylopic_name[tipdf$species[i]]),
                            error = function(e) NA_character_)
glyph_tips <- which(is.na(tipdf$uuid))
sil_poly <- if (length(glyph_tips)) do.call(rbind, lapply(glyph_tips, function(i) {
  g <- glyph[[tipdf$clade[i]]]; sh <- 0.92
  data.frame(species = tipdf$species[i],
             x = x_sil - sw/2 + g$x * sw, y = tipdf$y[i] - sh/2 + g$y * sh)
})) else data.frame(species = character(), x = numeric(), y = numeric())

# ---- assemble plot -------------------------------------------------
pal <- viridis(256, option = "turbo")
nu_breaks <- log(c(0.05, 0.2, 1, 5))
p <- ggplot() +
  geom_segment(data = seg_v, aes(x, y, xend = xend, yend = yend), colour = "grey55", linewidth = .5) +
  geom_segment(data = seg_h, aes(x, y, xend = xend, yend = yend, colour = val), linewidth = 2.4, lineend = "round") +
  # silhouettes (fallback polygons; rphylopic overrides below if online)
  geom_polygon(data = sil_poly, aes(x, y, group = species), fill = "grey20") +
  # growth-curve panels
  geom_rect(data = panel, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey96", colour = "grey85", linewidth = .3) +
  geom_segment(data = tipdf, aes(x = xc0, xend = xc0 + W, y = y - H/2, yend = y - H/2),
               colour = "grey80", linewidth = .3) +
  geom_ribbon(data = cdf, aes(x = px, ymin = plo, ymax = phi, group = species), fill = "grey70", alpha = .55) +
  geom_line(data = cdf, aes(x = px, y = pmed, group = species, colour = val), linewidth = .9) +
  geom_text(data = tipdf, aes(x = x_lab, y = y, label = lab), hjust = 0, fontface = "italic", size = 3.4) +
  annotate("text", x = xc0 + W/2, y = max(yy) + .75, label = "growth shape (size / asymptote  vs  relative age)",
           size = 3, colour = "grey35") +
  scale_colour_gradientn(colours = pal, name = expression(paste("shape  ", nu)),
                         breaks = nu_breaks, labels = c("0.05", "0.2", "1", "5"),
                         guide = guide_colourbar(barheight = 8, barwidth = .8)) +
  coord_cartesian(xlim = c(0, xlim_hi), ylim = c(0.4, max(yy) + 1.1), clip = "off") +
  labs(title = "Generalized-logistic growth across a dinosaur phylogeny",
       subtitle = "branch colour = shape parameter ν (ancestral REML reconstruction);  ribbon = 90% posterior of the fitted growth shape") +
  theme_void(base_size = 12) +
  theme(legend.position = c(.06, .5), plot.title = element_text(face = "bold"),
        plot.margin = margin(10, 16, 10, 10))

# real PhyloPic silhouettes on top (online only; glyph polygons cover the rest)
if (have_net) for (i in which(!is.na(tipdf$uuid)))
  p <- tryCatch(p + rphylopic::geom_phylopic(x = x_sil, y = tipdf$y[i], uuid = tipdf$uuid[i],
                                             height = 0.95, fill = "grey20"),
                error = function(e) p)

ggsave(file.path(out, "fig_phylo_growth.png"), p, width = 11, height = 6.2, dpi = 200, bg = "white")
cat("wrote outputs/fig_phylo_growth.png | silhouettes:",
    if (have_net) "PhyloPic (online)" else "built-in glyphs (offline)", "\n")
