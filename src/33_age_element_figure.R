#!/usr/bin/env Rscript
# =====================================================================
# 33_age_element_figure.R -- the main figure for model 32 (copy of 29 reading model 32 output).
#   left   : phylogeny, branches coloured by estimated growth rate. Node values
#            are the MODEL's node estimates (root + summed branch changes), not a
#            post-hoc reconstruction.
#   middle : silhouette, taxon, and the estimated growth curve at each tip
#            (relative size vs true age; median + 90% ribbon)
#   right  : ridgelines of the posterior for the key parameters per taxon:
#            shape nu, mid-growth rate r, age at half adult size t50
# Growth rate = the model's rate trait: proportion of adult size gained per
# year through mid-growth, 0.5 / (years from 25% to 75% of the asymptote). It
# stays comparable across curve shapes, which the Richards k does not once nu
# varies in sign.
# Usage: Rscript 33_age_element_figure.R [fit tag=broad] [data prefix=broad] [note printed on the figure]
# Colour: ONE colour per taxon = its mid-growth rate, used for its tip branch,
# silhouette, growth curve and all three ridgelines. Internal branches take the
# model's node estimate. Dark background; diverging palette centred on the median.
# =====================================================================
suppressMessages({ library(cmdstanr); library(posterior); library(ape); library(ggplot2)
  library(dplyr); library(readr); library(scales); library(rphylopic) })
a <- commandArgs(trailingOnly=TRUE)
tag  <- paste0("_", ifelse(is.na(a[1]), "broad", a[1])); pre <- ifelse(is.na(a[2]), "broad", a[2]); note <- ifelse(is.na(a[3]), "", a[3])
root <- "/Users/jmitchell/Library/CloudStorage/Dropbox/Research/00_current_Histology"
clean <- file.path(root,"data_clean"); out <- file.path(root,"outputs")
DIVERGING <- c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac')
rate_cols <- rev(DIVERGING)                       # slow -> blue, fast -> red
BG <- '#14161b'; FG <- '#e8e8e8'; FG2 <- '#a9adb5'; RULE <- '#3a3e47'   # dark mode

fit <- readRDS(file.path(out, paste0("fit_age_element", tag, ".rds")))
idx <- readRDS(file.path(out, paste0("age_element_index", tag, ".rds")))
sp <- idx$species; S <- length(sp); tm <- idx$tm; P <- tm$P; edge <- tm$edge; B <- tm$B; T <- tm$T
tos <- idx$tip_of_sp
dm <- function(v) as_draws_matrix(fit$draws(v))
delta <- array(as.numeric(dm("delta")), c(nrow(dm("mu_logr")), B, 4)); nd <- dim(delta)[1]
mu_r <- as.numeric(dm("mu_logr"))
logr_sp <- dm("logr_sp"); nu_sp <- dm("nu_sp"); t50_sp <- dm("t50_sp")

# ---- curve helpers ---------------------------------------------------
gfun <- function(f, nu) ifelse(abs(nu) < 1e-6, log(-log(f)), log((f^(-nu) - 1) / nu))
rate_lab <- "mid-growth rate\n(proportion of adult size per yr)"

# ---- model node values: root + summed changes on the path to each node --
nn <- T + length(unique(edge[,1])); PN <- matrix(0, nn, B)
for (n in seq_len(nn)) { node <- n
  repeat { e <- which(edge[,2]==node); if (!length(e)) break; PN[n,e] <- 1; node <- edge[e,1] } }
node_rate <- apply(mu_r + matrix(delta[,,2], nd) %*% t(PN), 2, median)   # log mid-growth rate, draws x nodes
sp_rate   <- apply(logr_sp, 2, median)
desc_tips <- function(node){ if (node<=T) return(node); unlist(lapply(edge[edge[,1]==node,2], desc_tips)) }
model_sets <- lapply(seq_len(nn), function(n) sort(desc_tips(n)))

# ---- display tree: one tip per SPECIES (tips holding 2 species become a cherry)
nwk <- write.tree(read.tree(file.path(clean, paste0(pre, "_tree.nwk"))))
for (tp in tm$tips) { spp <- sp[tos == match(tp, tm$tips)]
  rep <- if (length(spp) > 1) paste0("(", paste(spp, collapse=","), ")") else spp
  nwk <- gsub(paste0("([(,])", tp, "([,)])"), paste0("\\1", rep, "\\2"), nwk) }
dtree <- compute.brlen(read.tree(text=nwk), method="Grafen", power=0.55)   # power < 1 spreads the shallow nodes
pdf(NULL); plot(dtree, plot=FALSE); pp <- get("last_plot.phylo", envir=ape::.PlotPhyloEnv); invisible(dev.off())
xx <- pp$xx; yy <- pp$yy; nt <- length(dtree$tip.label); de <- dtree$edge; span <- max(xx[1:nt])
ddesc <- function(node){ if (node<=nt) return(node); unlist(lapply(de[de[,1]==node,2], ddesc)) }
dval <- sapply(seq_len(nt + dtree$Nnode), function(n){
  if (n <= nt) return(sp_rate[match(dtree$tip.label[n], sp)])
  tipset <- sort(unique(tos[match(dtree$tip.label[ddesc(n)], sp)]))
  node_rate[which(sapply(model_sets, identical, tipset))[1]] })
seg_h <- data.frame(x=xx[de[,1]], xend=xx[de[,2]], y=yy[de[,2]], val=dval[de[,2]])
seg_v <- do.call(rbind, lapply(unique(de[,1]), function(p){ ch <- de[de[,1]==p,2]
  data.frame(x=xx[p], y=min(yy[ch]), yend=max(yy[ch]), val=dval[p]) }))
midv <- median(dval); half <- max(quantile(abs(dval - midv), 0.95), 0.15)
colfun <- function(v) gradient_n_pal(rate_cols)(pmin(pmax((v - midv)/(2*half) + 0.5, 0), 1))

tipdf <- data.frame(species=dtree$tip.label, y=yy[1:nt]) %>%
  mutate(j=match(species, sp), val=sp_rate[j], col=colfun(val),
         lab=gsub("_"," ", sub("^cf\\.","cf. ", species)), genus=tm$tips[tos[j]])

# ---- column geometry -------------------------------------------------
x_sil <- span*1.24; x_lab <- span*1.52
xc0 <- span*2.40; Wc <- span*0.34; Hc <- 0.86             # growth-curve panel
gap <- span*0.10; Wr <- span*0.50                          # ridgeline columns
xr0 <- xc0 + Wc + gap*1.6 + (0:2)*(Wr + gap)
AGE_MAX <- 30; ages <- seq(0, AGE_MAX, length.out=80)

rich_rel <- function(t, r, nu, t50){ k <- 2*r*(gfun(.25,nu) - gfun(.75,nu)); x <- pmin(-k*(t - t50) + gfun(.5,nu), 50)
  u <- nu*exp(x); ifelse(abs(u) < 1e-8, exp(-exp(x)), ifelse(1+u <= 0, 0, (1+u)^(-1/nu))) }
curves <- do.call(rbind, lapply(seq_len(nt), function(i){ j <- tipdf$j[i]
  rr <- exp(as.numeric(logr_sp[,j])); nu <- as.numeric(nu_sp[,j]); t50 <- as.numeric(t50_sp[,j])
  r <- sapply(ages, function(t) rich_rel(t, rr, nu, t50))             # draws x ages
  q <- apply(r, 2, quantile, c(.05,.5,.95))
  data.frame(species=tipdf$species[i], col=tipdf$col[i], px=xc0 + ages/AGE_MAX*Wc,
             lo=tipdf$y[i] + (q[1,]-.5)*Hc, med=tipdf$y[i] + (q[2,]-.5)*Hc, hi=tipdf$y[i] + (q[3,]-.5)*Hc) }))

# ridgelines: posterior density per species, scaled to a common peak height
ridge_spec <- list(
  list(name="shape ν", draws=function(j) as.numeric(nu_sp[,j]), lim=c(-1.5, 3), brk=c(-1, 0, 1, 2, 3), lab=identity,
       ref=c(-1, -1/3, 0, 1), reflab=c("M","vB","G","L")),
  list(name="mid-growth rate (per yr)", draws=function(j) as.numeric(logr_sp[,j]), lim=log(c(0.01, 1)), brk=log(c(0.01,0.03,0.1,0.3,1)),
       lab=function(b) signif(exp(b),1), ref=NULL, reflab=NULL),
  list(name="age at half adult size (yr)", draws=function(j) as.numeric(t50_sp[,j]), lim=c(-5, 20), brk=c(-5,0,5,10,15,20), lab=identity,
       ref=0, reflab=NULL))
ridges <- list(); axes <- list(); refs <- list(); reftxt <- list()
for (c in seq_along(ridge_spec)) { rs <- ridge_spec[[c]]; sx <- function(v) xr0[c] + (v - rs$lim[1])/diff(rs$lim)*Wr
  for (i in seq_len(nt)) { v <- rs$draws(tipdf$j[i]); v <- v[v >= rs$lim[1] & v <= rs$lim[2]]
    if (length(v) < 20) next
    d <- density(v, from=rs$lim[1], to=rs$lim[2], n=160)
    ridges[[length(ridges)+1]] <- data.frame(grp=paste(c,i), col=tipdf$col[i], x=sx(d$x),
      ymin=tipdf$y[i] - 0.38, ymax=tipdf$y[i] - 0.38 + 1.25*d$y/max(d$y), ord=-tipdf$y[i]) }
  axes[[c]] <- data.frame(x=sx(rs$brk), lab=as.character(rs$lab(rs$brk)), x0=xr0[c], x1=xr0[c]+Wr, title=rs$name, xm=xr0[c]+Wr/2)
  if (!is.null(rs$ref)) refs[[c]] <- data.frame(x=sx(rs$ref))
  if (!is.null(rs$reflab)) reftxt[[c]] <- data.frame(x=sx(rs$ref), lab=rs$reflab) }
ridges <- bind_rows(ridges) %>% arrange(ord); axes <- bind_rows(axes); refs <- bind_rows(refs); reftxt <- bind_rows(reftxt)
ridges$grp <- factor(ridges$grp, levels=unique(ridges$grp))
ybot <- 0.35; ytop <- nt + 1.0
cax <- data.frame(x=xc0 + c(0,10,20,30)/AGE_MAX*Wc, lab=c(0,10,20,30))

# colour legend (built by hand so everything shares one coordinate system)
lg <- data.frame(v=seq(midv-half, midv+half, length.out=60)) %>% mutate(x=span*0.02 + (seq_along(v)-1)/59*span*0.42)
lgt <- data.frame(v=signif(exp(midv + c(-0.8, 0, 0.8)*half), 2)) %>% filter(log(v) >= midv-half, log(v) <= midv+half) %>%
  mutate(x=span*0.02 + (log(v)-(midv-half))/(2*half)*span*0.42)

p <- ggplot() +
  geom_segment(data=refs, aes(x=x, xend=x, y=ybot, yend=ytop-0.3), colour=RULE, linewidth=0.3) +
  geom_segment(data=seg_v, aes(x=x, xend=x, y=y, yend=yend), colour=colfun(seg_v$val), linewidth=1.5, lineend="round") +
  geom_segment(data=seg_h, aes(x=x, xend=xend, y=y, yend=y), colour=colfun(seg_h$val), linewidth=1.5, lineend="round") +
  geom_text(data=tipdf, aes(x=x_lab, y=y, label=lab), hjust=0, size=3.1, fontface="italic", colour=FG) +
  # growth curves
  geom_segment(data=tipdf, aes(x=xc0, xend=xc0+Wc, y=y-Hc/2, yend=y-Hc/2), colour=RULE, linewidth=0.25) +
  geom_ribbon(data=curves, aes(x=px, ymin=lo, ymax=hi, group=species, fill=col), alpha=0.35) +
  geom_line(data=curves, aes(x=px, y=med, group=species, colour=col), linewidth=0.7) +
  # ridgelines
  geom_ribbon(data=ridges, aes(x=x, ymin=ymin, ymax=ymax, group=grp, fill=col), colour=BG, linewidth=0.25, alpha=0.9) +
  scale_fill_identity() + scale_colour_identity() +
  # axes
  geom_segment(data=axes %>% distinct(x0, x1), aes(x=x0, xend=x1, y=ybot, yend=ybot), colour=FG2, linewidth=0.3) +
  geom_text(data=axes, aes(x=x, y=ybot-0.35, label=lab), size=2.5, colour=FG2) +
  geom_text(data=axes %>% distinct(xm, title), aes(x=xm, y=ybot-1.05, label=title), size=3, colour=FG) +
  geom_text(data=reftxt, aes(x=x, y=ytop-0.05, label=lab), size=2.4, colour=FG2) +
  annotate("segment", x=xc0, xend=xc0+Wc, y=ybot, yend=ybot, colour=FG2, linewidth=0.3) +
  geom_text(data=cax, aes(x=x, y=ybot-0.35, label=lab), size=2.5, colour=FG2) +
  annotate("text", x=xc0+Wc/2, y=ybot-1.05, label="age (yr)", size=3, colour=FG) +
  # legend
  geom_tile(data=lg, aes(x=x, y=ybot-0.75), fill=colfun(lg$v), height=0.35, width=span*0.42/59*1.05) +
  geom_text(data=lgt, aes(x=x, y=ybot-1.25, label=v), size=2.5, colour=FG2) +
  annotate("text", x=span*0.23, y=ybot-2.05, label=rate_lab, size=3, colour=FG, lineheight=0.9) +
  annotate("text", x=span*0.02, y=ybot-2.95, label="one colour per taxon: branch, silhouette, curve, ridgelines", size=2.6, colour=FG2, hjust=0) +
  annotate("text", x=max(xr0)+Wr, y=ybot-2.95, label=note, size=3.2, colour="#f4a582", hjust=1) +
  coord_cartesian(xlim=c(-span*0.02, max(xr0)+Wr+span*0.05), ylim=c(ybot-3.2, ytop), clip="off") +
  theme_void(base_size=12) + theme(plot.background=element_rect(fill=BG, colour=NA))
for (i in seq_len(nt)) {
  f <- file.path(clean, "phylopic", paste0(tipdf$genus[i], ".rds"))
  if (!file.exists(f) && tipdf$genus[i]=="Diplocynodon") f <- file.path(clean, "phylopic", "Alligator.rds")
  if (!file.exists(f)) next
  img <- readRDS(f); asp <- dim(img)[2]/dim(img)[1]            # long, low animals are capped in width; small squat ones enlarged
  p <- tryCatch(p + add_phylopic(img=img, x=x_sil, y=tipdf$y[i], height=min(0.80, 0.62*3.2/max(asp, 1.6)), fill=tipdf$col[i]),
                                    error=function(e) p)
}
fn <- file.path(out, paste0("fig_age_element", tag, ".png"))
ggsave(fn, p, width=15, height=max(7, 0.44*nt + 2), dpi=200, bg=BG, limitsize=FALSE)
cat("wrote", fn, "\n")
