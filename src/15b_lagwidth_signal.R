#!/usr/bin/env Rscript
# Lightweight (no-Stan) LAG-width growth-rate analysis -- robust under memory
# pressure. Per-specimen log(radial LAG width) ~ age regressions give a growth-
# rate LEVEL (rate at a reference age) and a DECLINE slope; phylogenetic signal
# in the growth-rate level is estimated with Pagel's lambda + Blomberg's K.
# Circumference taxa only (radial mm/yr comparable), active-growth rings only.
suppressMessages({ library(dplyr); library(readr); library(ape); library(phytools); library(ggplot2) })
root  <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
clean <- file.path(root,"data_clean"); out <- file.path(root,"outputs")
REF <- 2  # reference age (yr) for the reported rate level

w <- read_csv(file.path(clean,"lag_widths_long.csv"), show_col_types=FALSE) %>%
  filter(dimension=="circumference_mm", !is_efs, lag_width>0) %>%
  mutate(logw=log(lag_width))

# per-specimen regression -> rate@REF (log) and decline slope
spec <- w %>% group_by(species, specimen) %>% filter(n()>=2) %>%
  summarise(b=coef(lm(logw~age))[2], a=coef(lm(logw~age))[1], n=n(), .groups="drop") %>%
  mutate(lograte_ref=a+b*REF)
# species-level summaries
sp <- spec %>% group_by(species) %>%
  summarise(lograte_ref=mean(lograte_ref), decline=mean(b), n_spec=n(), .groups="drop") %>%
  mutate(rate_mm_yr=exp(lograte_ref))

# map to tree tips (mean per tip), prune tree, phylo signal
tipmap <- read_csv(file.path(clean,"species_to_tip.csv"), show_col_types=FALSE)
sp$tip <- tipmap$phylo_tip[match(sp$species, tipmap$species)]
tip_val <- tapply(sp$lograte_ref, sp$tip, mean)
tree <- compute.brlen(keep.tip(read.tree(file.path(clean,"dino_tree.nwk")), names(tip_val)), method="Grafen")
x <- tip_val[tree$tip.label]
lam <- phylosig(tree, x, method="lambda", test=TRUE)
kk  <- phylosig(tree, x, method="K", test=TRUE)

cat("=== LAG-width growth-rate level: phylogenetic signal ===\n")
cat(sprintf("Pagel's lambda = %.2f (p = %.3f) | Blomberg's K = %.2f (p = %.3f)\n",
            lam$lambda, lam$P, kk$K, kk$P))
cat("\n=== radial apposition rate @ age", REF, "(mm/yr) + decline (per yr, log scale) ===\n")
print(sp %>% arrange(desc(rate_mm_yr)) %>%
        transmute(species, rate_mm_yr=round(rate_mm_yr,2), decline=round(decline,2), n_spec) %>% as.data.frame(),
      row.names=FALSE)
write_csv(sp, file.path(clean,"lagwidth_rates_by_species.csv"))

# ---- figures (house style: no titles, white panels) ----
DIVERGING <- c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac')
disp <- c(Tyrannosaurus_rex="T. rex", cf.Tyrannosaurus_rex="cf. T. rex (Jane)", Gorgosaurus_libratus="Gorgosaurus",
  Rativates_evadens="Rativates", Campanas_Carcharodontosaurid="Meraxes", Acrocanthosaurus_atokensis="Acrocanthosaurus",
  Allosaurus_fragilis="Allosaurus", Maiasaura_peeblesorum="Maiasaura", Alligator_mississippiensis="Alligator")
w$lab <- ifelse(w$species %in% names(disp), disp[w$species], w$species)
th <- theme_minimal(base_size=10)+theme(panel.grid.minor=element_blank(), panel.background=element_blank())

p1 <- ggplot(w, aes(age, lag_width, group=specimen)) +
  geom_point(size=.8, colour="grey30", alpha=.6) +
  geom_smooth(aes(group=lab), method="lm", se=FALSE, colour="#b2182b", linewidth=.7, formula=y~x) +
  scale_y_log10() + facet_wrap(~lab, scales="free") +
  labs(x="age (yr)", y="radial LAG width (mm/yr, log scale)") + th
ggsave(file.path(out,"fig_lagwidth_trajectories.png"), p1, width=10, height=6, dpi=140, bg="white")

spp <- sp %>% mutate(lab=ifelse(species %in% names(disp), disp[species], species)) %>% arrange(rate_mm_yr) %>%
  mutate(lab=factor(lab, levels=lab))
p2 <- ggplot(spp, aes(rate_mm_yr, lab, fill=lograte_ref)) +
  geom_col() + scale_fill_gradientn(colours=rev(DIVERGING), guide="none") +
  labs(x=paste0("radial apposition rate @ age ",REF," (mm/yr)"), y=NULL) + th
ggsave(file.path(out,"fig_lagwidth_rates.png"), p2, width=6.5, height=4, dpi=140, bg="white")
cat("\nwrote fig_lagwidth_trajectories.png, fig_lagwidth_rates.png, lagwidth_rates_by_species.csv\n")
