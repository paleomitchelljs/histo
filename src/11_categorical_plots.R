#!/usr/bin/env Rscript
# Summary plots for the categorical histology characters and skeletal elements:
#   (1) data coverage / structure  -> fig_categorical_coverage.png
#   (2) how the model USES them     -> fig_categorical_effects.png
# House style: categorical palette, no titles, white panels.
suppressMessages({ library(cmdstanr); library(dplyr); library(readr); library(ggplot2); library(tidyr)
  library(patchwork) })
root <- "/Users/jmitchell/Library/CloudStorage/Dropbox/Research/00_current_Histology"
clean <- file.path(root,"data_clean"); out <- file.path(root,"outputs")
CATEG <- c('#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#bc80bd','#ccebc5','#d9d9d9','#ffed6f')
th <- theme_minimal(base_size=11) +
  theme(panel.grid.minor=element_blank(), panel.background=element_blank(),
        plot.background=element_blank())

dat <- read_csv(file.path(clean,"lag_size_long.csv"), show_col_types=FALSE)
spec <- dat %>% group_by(specimen) %>%
  summarise(clade=first(clade), element=first(element), dimension=first(dimension),
            resorption=first(resorption), efs_present=first(efs_present),
            n_points=n(), age_min=min(age), age_max=max(age), maxsize=max(size), .groups="drop")

# ---- (1) coverage --------------------------------------------------
# 1a: element x clade -- the confound (block diagonal)
ec <- spec %>% count(clade, element)
p1a <- ggplot(ec, aes(element, clade, fill=element)) +
  geom_tile(colour="white", linewidth=1) +
  geom_text(aes(label=n), fontface="bold") +
  scale_fill_manual(values=CATEG, guide="none") +
  labs(x="element", y=NULL, tag="A") +
  th + theme(axis.text.x=element_text(angle=20, hjust=1))

# 1b: categorical metadata actually used by the model (resorption x EFS)
re <- spec %>% count(resorption, efs_present) %>%
  mutate(efs_present=factor(efs_present, c(FALSE,TRUE), c("no EFS","EFS present")))
p1b <- ggplot(re, aes(resorption, n, fill=efs_present)) +
  geom_col(position=position_dodge(.7), width=.6) +
  geom_text(aes(label=n), position=position_dodge(.7), vjust=-.3, size=3.3) +
  scale_fill_manual(values=CATEG[c(3,4)], name=NULL) +
  labs(x="resorption category", y="specimens", tag="B") + th

# 1c: sampling effort (LAG points) by element, colored by clade
p1c <- ggplot(spec, aes(reorder(element, n_points, sum), n_points, fill=clade)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values=CATEG, name="clade") +
  labs(x="element", y="LAG / annulus observations", tag="C") + th

cov <- p1a + p1b + p1c + plot_layout(widths=c(1.1,1,1.2))
ggsave(file.path(out,"fig_categorical_coverage.png"), cov, width=13, height=4.2, dpi=150, bg="white")
cat("wrote fig_categorical_coverage.png\n")

# ---- (2) how the model uses the categoricals -----------------------
fit <- readRDS(file.path(out,"fit_offset_phylo.rds"))
idx <- readRDS(file.path(out,"offset_phylo_index.rds"))
d   <- fit$draws(format="df")
it  <- idx$ind_tbl                                   # specimen order in the model
gv  <- function(par,i) d[[paste0(par,"[",i,"]")]]

eff <- lapply(seq_len(nrow(it)), function(i){
  miss <- gv("miss_yr",i); A <- exp(gv("logA_i",i))
  ms <- spec$maxsize[match(as.character(it$specimen[i]), spec$specimen)]
  data.frame(specimen=as.character(it$specimen[i]),
             clade=spec$clade[match(as.character(it$specimen[i]),spec$specimen)],
             resorption=it$resorption[i], efs=it$efs_present[i],
             miss=median(miss), miss_lo=quantile(miss,.05), miss_hi=quantile(miss,.95),
             ratio=median(A/ms), ratio_lo=quantile(A/ms,.05), ratio_hi=quantile(A/ms,.95))
}) %>% bind_rows()

# 2a: resorption category -> estimated missing inner rings
p2a <- ggplot(eff, aes(miss, resorption, colour=clade)) +
  geom_linerange(aes(xmin=miss_lo, xmax=miss_hi), position=position_jitter(height=.18,seed=1), alpha=.6) +
  geom_point(position=position_jitter(height=.18,seed=1), size=2.4) +
  scale_colour_manual(values=CATEG, name="clade") +
  labs(x="estimated missing inner rings (yr)", y="resorption", tag="A",
       subtitle=NULL) + th

# 2b: EFS status -> asymptote anchored to largest ring (A / max observed size)
eff$efslab <- factor(eff$efs, c(FALSE,TRUE), c("no EFS\n(still growing)","EFS present\n(growth ceased)"))
p2b <- ggplot(eff, aes(efslab, ratio, colour=clade)) +
  geom_hline(yintercept=1, linetype=2, colour="grey50") +
  geom_linerange(aes(ymin=ratio_lo, ymax=ratio_hi), position=position_jitter(width=.12,seed=2), alpha=.6) +
  geom_point(position=position_jitter(width=.12,seed=2), size=2.4) +
  scale_colour_manual(values=CATEG, guide="none") +
  labs(x=NULL, y="asymptote A / largest observed ring", tag="B") + th

ggsave(file.path(out,"fig_categorical_effects.png"), p2a + p2b + plot_layout(widths=c(1.2,1)),
       width=11, height=4.4, dpi=150, bg="white")
cat("wrote fig_categorical_effects.png\n")

# small console summary
cat("\nEFS present: A/maxsize =", round(median(eff$ratio[eff$efs]),2),
    "| no EFS: A/maxsize =", round(median(eff$ratio[!eff$efs]),2), "\n")
cat("resorption absent: miss =", round(median(eff$miss[eff$resorption=='absent']),2),
    "| present: miss =", round(median(eff$miss[eff$resorption!='absent']),2), "\n")
