#!/usr/bin/env Rscript
# Salvage branch-shift recovery from the persisted cmdstan CSV (outputs/bs_csv),
# in case the R post-processing was OOM-killed after sampling started.
suppressMessages({ library(dplyr); library(readr); library(ape) })
root <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
clean <- file.path(root,"data_clean"); sim <- file.path(clean,"sim"); out <- file.path(root,"outputs")
csv <- list.files(file.path(out,"bs_csv"), pattern="\\.csv$", full.names=TRUE)
csv <- csv[which.max(file.info(csv)$mtime)]
cat("reading", basename(csv), "\n")
d <- read.csv(csv, comment.char="#")
truth_edges <- read_csv(file.path(sim,"scenario_shift_truth_edges.csv"), show_col_types=FALSE)
spl <- sort(unique(read_csv(file.path(sim,"scenario_shift_data.csv"),show_col_types=FALSE)$species))
tree <- compute.brlen(keep.tip(read.tree(file.path(clean,"dino_tree.nwk")), spl), method="Grafen")
B <- nrow(tree$edge)
dc <- d[, grep("^delta_out\\.", names(d)), drop=FALSE]
cat("draws salvaged:", nrow(dc), "| branches:", ncol(dc), "\n")
delta <- data.frame(edge=1:B,
  med=sapply(1:B,function(b) median(dc[[b]])), absmed=sapply(1:B,function(b) median(abs(dc[[b]]))),
  p_big=sapply(1:B,function(b) mean(abs(dc[[b]])>0.5)), child_node=tree$edge[,2])
delta$is_true_shift <- delta$edge %in% truth_edges$edge
delta <- delta %>% left_join(truth_edges %>% transmute(edge, true_jump=jump), by="edge")
dr <- delta %>% arrange(desc(absmed)); write_csv(delta, file.path(out,"sim_branchshift_delta.csv"))
cat("\ntrue shift edges:", paste(truth_edges$edge,collapse=","), "jumps:", paste(round(truth_edges$jump,2),collapse=","), "\n")
cat("--- top 6 by |median shift| ---\n")
print(head(dr %>% transmute(edge, child_node, median_shift=round(med,2), abs_med=round(absmed,2), p_excl=round(p_big,2),
  TRUE_SHIFT=ifelse(is_true_shift,"<==TRUE","")),6), row.names=FALSE)
cat("--- true edges ---\n")
print(delta %>% filter(is_true_shift) %>% transmute(edge, true_jump=round(true_jump,2), median_shift=round(med,2),
  p_abs_gt_0.5=round(p_big,2), rank=match(edge,dr$edge)), row.names=FALSE)
