#!/usr/bin/env Rscript
# =====================================================================
# 10_figures.R  -- figures for the offset + categorical model on the full
# cross-Dinosauria dataset. House style:
#   * diverging palette for continuous params (nu), centered at nu=1 (logistic)
#   * categorical palette for groups
#   * NO titles/subtitles; no grey panel backgrounds
#   * ribbons color-matched: dark line + low-opacity same-color ribbon
#   * real PhyloPic silhouettes (cached in data_clean/phylopic/)
# =====================================================================
suppressMessages({ library(cmdstanr); library(ape); library(ggplot2); library(dplyr)
  library(readr); library(grid); library(scales); library(rphylopic); library(magick) })
root <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
clean <- file.path(root,"data_clean"); out <- file.path(root,"outputs")

DIVERGING <- c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac')
CATEG     <- c('#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#bc80bd','#ccebc5','#d9d9d9','#ffed6f')
nu_cols   <- rev(DIVERGING)          # low nu -> blue, nu=1 -> white, high nu -> red

dat  <- read_csv(file.path(clean,"lag_size_long.csv"), show_col_types=FALSE)
idx  <- readRDS(file.path(out,"offset_phylo_index.rds"))
fit  <- readRDS(file.path(out,"fit_offset_phylo.rds"))
d    <- fit$draws(format="df")
sp   <- idx$species                  # species in model order
tree <- read.tree(file.path(clean,"dino_tree.nwk"))
tipmap <- read_csv(file.path(clean,"species_to_tip.csv"), show_col_types=FALSE)

gv <- function(par,j) d[[paste0(par,"[",j,"]")]]      # draws of par[j]
lognu_med <- sapply(seq_along(sp), function(j) median(gv("lognu_sp",j)))
names(lognu_med) <- sp
disp <- c(Tyrannosaurus_rex="Tyrannosaurus rex", cf.Tyrannosaurus_rex="cf. Tyrannosaurus (Jane)",
  Gorgosaurus_libratus="Gorgosaurus", Rativates_evadens="Rativates",
  Campanas_Carcharodontosaurid="Meraxes gigas", Acrocanthosaurus_atokensis="Acrocanthosaurus",
  Allosaurus_fragilis="Allosaurus", Maiasaura_peeblesorum="Maiasaura",
  Plateosaurus_engelhardti="Plateosaurus", Mamenchisauridae_indet="Mamenchisaurid",
  Apatosaurus_sp="Apatosaurus", Diplodocidae_indet="Diplodocid", Camarasaurus_sp="Camarasaurus")
dnm <- function(s) ifelse(s %in% names(disp), disp[s], gsub("_"," ",s))   # generic display name

# tip-level lognu (map species->tip; for tips with 2 species e.g. Tyrannosaurus, average)
tip_lognu <- sapply(tree$tip.label, function(tp){
  spp <- tipmap$species[tipmap$phylo_tip==tp]; mean(lognu_med[spp], na.rm=TRUE) })
# diverging colour centred on the MEDIAN taxon with robust half-range, so the
# (mostly sub-logistic) taxa span the full palette instead of all reading blue.
allv <- c(tip_lognu, lognu_med); midv <- median(allv)
devv <- max(as.numeric(quantile(abs(allv - midv), 0.9)), 0.3)
clo <- midv - devv; chi <- midv + devv
nu_breaks <- c(clo, midv, chi); nu_labs <- sprintf("%.2f", exp(nu_breaks))

# ---- ancestral reconstruction of lognu on the tree ----------------
treeU <- compute.brlen(tree, method="Grafen")
ace_v <- ape::ace(tip_lognu[treeU$tip.label], treeU, type="continuous", method="REML")$ace
nodeval <- c(tip_lognu[treeU$tip.label], ace_v)
pdf(NULL); plot(treeU, plot=FALSE); pp<-get("last_plot.phylo",envir=ape::.PlotPhyloEnv); invisible(dev.off())
xx<-pp$xx; yy<-pp$yy; ntip<-length(tree$tip.label); span<-max(xx[1:ntip]); tipx<-span
edge<-treeU$edge
seg_h<-data.frame(x=xx[edge[,1]],xend=xx[edge[,2]],y=yy[edge[,2]],yend=yy[edge[,2]],
                  val=(nodeval[edge[,1]]+nodeval[edge[,2]])/2)
seg_v<-do.call(rbind,lapply(unique(edge[,1]),function(p){ch<-edge[edge[,1]==p,2]
  data.frame(x=xx[p],xend=xx[p],y=min(yy[ch]),yend=max(yy[ch]))}))

# ---- per-tip growth-shape ribbons (normalized; from species params) -
ns <- 400; ridx <- sample(nrow(d), ns)
shape_for_tip <- function(tp){
  spp <- tipmap$species[tipmap$phylo_tip==tp]
  j <- which(sp %in% spp)[1]                      # one species per tip (Tyrannosaurus: use T. rex)
  amax <- max(dat$age[dat$species %in% spp]) + 1
  ag <- seq(0, amax, length.out=50)
  nu <- exp(gv("lognu_sp",j)[ridx]); k<-exp(gv("logk_sp",j)[ridx]); t0<-gv("t0_sp",j)[ridx]
  M <- sapply(ag, function(a) (1+nu*exp(-k*(a-t0)))^(-1/nu))   # rel size, ns x len
  data.frame(tip=tp, rel_age=ag/amax,
             med=apply(M,2,median), lo=apply(M,2,quantile,.05), hi=apply(M,2,quantile,.95),
             val=tip_lognu[tp])
}
curve_df <- do.call(rbind, lapply(tree$tip.label, shape_for_tip))

# layout columns
x_sil<-tipx+0.18*span; sil_W<-0.30*span; sil_H<-0.82; x_lab<-tipx+0.36*span
xc0<-tipx+0.95*span; W<-0.55*span; H<-0.80; xlim_hi<-xc0+W+0.06*span
tipdf<-data.frame(tip=tree$tip.label, x=tipx, y=yy[1:ntip],
                  lab=disp[ tipmap$species[match(tree$tip.label,tipmap$phylo_tip)] ],
                  val=tip_lognu[tree$tip.label])
# labels = genus tip names (clean), with a couple of overrides
tipdf$lab <- ifelse(tipdf$tip=="Northampton", "Northampton sauropod",
             ifelse(tipdf$tip=="Meraxes", "Meraxes gigas", tipdf$tip))

cdf <- curve_df %>% left_join(tipdf[,c("tip","y")], by="tip") %>%
  mutate(px=xc0+rel_age*W, pmed=y+(med-.5)*H, plo=y+(lo-.5)*H, phi=y+(hi-.5)*H)

# ---- assemble composite -------------------------------------------
p <- ggplot() +
  geom_segment(data=seg_v, aes(x,y,xend=xend,yend=yend), colour="grey55", linewidth=.5) +
  geom_segment(data=seg_h, aes(x,y,xend=xend,yend=yend,colour=val), linewidth=2.6, lineend="round")
# baseline guide per curve (thin), then color-matched ribbon + line
p <- p +
  geom_segment(data=tipdf, aes(x=xc0,xend=xc0+W,y=y-H/2,yend=y-H/2), colour="grey88", linewidth=.3) +
  geom_ribbon(data=cdf, aes(x=px,ymin=plo,ymax=phi,group=tip,fill=val), alpha=.35) +
  geom_line(data=cdf, aes(x=px,y=pmed,group=tip,colour=val), linewidth=1.2) +
  geom_text(data=tipdf, aes(x=x_lab,y=y,label=lab), hjust=0, fontface="italic", size=3.5) +
  scale_colour_gradientn(colours=nu_cols, limits=c(clo,chi), oob=scales::squish,
    breaks=nu_breaks, labels=nu_labs,
    name=expression(nu), guide=guide_colourbar(barheight=9,barwidth=.8)) +
  scale_fill_gradientn(colours=nu_cols, limits=c(clo,chi), oob=scales::squish, guide="none") +
  coord_cartesian(xlim=c(0,xlim_hi), ylim=c(0.3,ntip+0.9), clip="off") +
  theme_void(base_size=12) + theme(legend.position=c(.05,.45))

# silhouettes via add_phylopic (correct aspect/coords), uniform height for
# comparability; flip_tips mirrors taxa (via flip_phylopic) for consistent facing.
# theropod silhouettes face right; flip them so all taxa face left (toward tree)
flip_tips <- c("Tyrannosaurus","Gorgosaurus","Albertosaurus","Saurornitholestes",
               "Rativates","Meraxes","Acrocanthosaurus","Allosaurus","Syntarsus")
for (i in seq_len(nrow(tipdf))) {
  f <- file.path(clean,"phylopic",paste0(tipdf$tip[i],".rds"))
  if (!file.exists(f)) next
  img <- readRDS(f)
  if (tipdf$tip[i] %in% flip_tips) img <- tryCatch(flip_phylopic(img, horizontal=TRUE), error=function(e) img)
  p <- tryCatch(p + add_phylopic(img=img, x=x_sil, y=tipdf$y[i], height=0.95, fill="grey15"),
                error=function(e) p)
}
ggsave(file.path(out,"fig_phylo_growth.png"), p, width=11.5, height=max(7, 0.52*ntip+1.5), dpi=200, bg="white", limitsize=FALSE)
cat("wrote fig_phylo_growth.png\n")

# ---- figure: estimated missing inner rings (the offset feature) ----
miss_df <- do.call(rbind, lapply(seq_len(nrow(idx$ind_tbl)), function(i){
  v<-gv("miss_yr",i); data.frame(specimen=idx$ind_tbl$specimen[i],
    resorption=idx$ind_tbl$resorption[i],
    m=median(v),lo=quantile(v,.05),hi=quantile(v,.95))})) %>%
  left_join(dat %>% distinct(specimen,clade), by="specimen") %>%
  arrange(resorption,m) %>% mutate(specimen=factor(specimen,levels=specimen))
pm <- ggplot(miss_df, aes(m, specimen, colour=clade)) +
  geom_linerange(aes(xmin=lo,xmax=hi)) + geom_point(size=2) +
  scale_colour_manual(values=CATEG, name="clade") +
  facet_grid(resorption~., scales="free_y", space="free_y") +
  labs(x="estimated missing inner rings (years)", y=NULL) +
  theme_minimal(base_size=10) + theme(panel.grid.minor=element_blank())
ggsave(file.path(out,"fig_offset_estimates.png"), pm, width=7.5, height=7, dpi=150, bg="white")
cat("wrote fig_offset_estimates.png\n")

# ---- figure: per-specimen fitted curves over data (all clades) -----
# fitted curve per specimen using its logA_i, miss, and its species k/nu/t0
ind_tbl <- idx$ind_tbl
fitcurve <- do.call(rbind, lapply(seq_len(nrow(ind_tbl)), function(i){
  s <- ind_tbl$sp_of_ind[i]; spc <- as.character(ind_tbl$specimen[i])
  dd <- dat[dat$specimen==spc,]; amax<-max(dd$age)
  ag <- seq(0, amax, length.out=40)
  logA<-gv("logA_i",i)[ridx]; nu<-exp(gv("lognu_sp",s)[ridx]); k<-exp(gv("logk_sp",s)[ridx])
  t0<-gv("t0_sp",s)[ridx]; mi<-gv("miss_yr",i)[ridx]
  M <- sapply(ag, function(a){ ta<-a+mi; exp(logA - (1/nu)*log1p(nu*exp(-k*(ta-t0)))) })
  data.frame(species=as.character(ind_tbl$species[i]), specimen=spc, age=ag,
             med=apply(M,2,median)) }))
pc <- ggplot() +
  geom_line(data=fitcurve, aes(age,med,group=specimen), colour="#b2182b", alpha=.8, linewidth=.6) +
  geom_point(data=dat, aes(age,size), size=.7, colour="grey25") +
  facet_wrap(~species, scales="free") +
  labs(x="age (years; recorded + estimated missing rings)", y="skeletal size (mm circ. or % bone length)") +
  theme_minimal(base_size=9) + theme(panel.grid.minor=element_blank())
ggsave(file.path(out,"fig_growth_curves_all.png"), pc, width=11, height=7.5, dpi=150, bg="white")
cat("wrote fig_growth_curves_all.png\n")

# ---- figure: per-species nu forest (diverging colored) -------------
nuf <- do.call(rbind, lapply(seq_along(sp), function(j){
  v<-exp(gv("lognu_sp",j)); data.frame(species=dnm(sp[j]), lognu_med=median(log(v)),
    m=median(v),lo=quantile(v,.05),hi=quantile(v,.95))})) %>% arrange(m) %>%
  mutate(species=factor(species,levels=species))
pf <- ggplot(nuf, aes(m, species, colour=lognu_med)) +
  geom_vline(xintercept=1, linetype=2, colour="grey50") +
  geom_linerange(aes(xmin=lo,xmax=hi), linewidth=1) + geom_point(size=2.6) +
  scale_colour_gradientn(colours=nu_cols, limits=c(clo,chi), oob=scales::squish, guide="none") +
  scale_x_log10() + labs(x=expression(paste("shape ",nu," (log scale; 1 = logistic)")), y=NULL) +
  theme_minimal(base_size=10) + theme(panel.grid.minor=element_blank())
ggsave(file.path(out,"fig_nu_forest.png"), pf, width=6.5, height=5, dpi=150, bg="white")
cat("wrote fig_nu_forest.png\nALL FIGURES DONE\n")
