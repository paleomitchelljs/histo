#!/usr/bin/env Rscript
# Clean summary plots for the Monte Carlo recovery study (outputs/mc_recovery_results.csv).
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(ggplot2) })
out <- file.path(local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) }), "outputs")
r <- read_csv(file.path(out,"mc_recovery_results.csv"), show_col_types=FALSE)
r$regime <- factor(r$regime, levels=c("iid","bm","shift"))
CAT <- c(iid='#80b1d3', bm='#fdb462', shift='#fb8072')
cat("reps per regime:\n"); print(table(r$regime))

# Panel A: recovered phylo signal lambda_nu by regime, with true-trait lambda
truA <- r %>% group_by(regime) %>% summarise(lam_true=median(lam_true), .groups="drop")
pA <- ggplot(r, aes(regime, lambda_rec, fill=regime)) +
  geom_violin(alpha=.5, colour=NA) + geom_boxplot(width=.15, outlier.size=.5, alpha=.8) +
  geom_point(data=truA, aes(regime, lam_true), shape=23, size=3, fill="black", inherit.aes=FALSE) +
  scale_fill_manual(values=CAT) + ylim(0,1) +
  labs(x=NULL, y=expression("recovered "*lambda[nu])) +
  theme_classic(base_size=12) + theme(legend.position="none")

# Panel B: recovery quality (correlation true vs recovered) per quantity x regime
corr <- r %>% select(regime, nu_cor, offset_cor, neonate_cor) %>%
  pivot_longer(-regime, names_to="quantity", values_to="cor") %>%
  mutate(quantity=recode(quantity, nu_cor="shape ν", offset_cor="missing-LAG offset", neonate_cor="neonate mass"))
corr$quantity <- factor(corr$quantity, levels=c("missing-LAG offset","neonate mass","shape ν"))
pB <- ggplot(corr, aes(quantity, cor, fill=regime)) +
  geom_boxplot(outlier.size=.4, alpha=.8) + scale_fill_manual(values=CAT) +
  geom_hline(yintercept=c(0,1), colour="grey80", linetype=3) +
  labs(x=NULL, y="recovery correlation (true vs recovered)") +
  theme_classic(base_size=12) + theme(legend.position=c(.12,.2), legend.title=element_blank())

# Panel C: missing-LAG offset bias (recovered - true), should be ~0
pC <- ggplot(r, aes(regime, offset_bias, fill=regime)) +
  geom_hline(yintercept=0, colour="grey60", linetype=2) +
  geom_violin(alpha=.5, colour=NA) + geom_boxplot(width=.15, outlier.size=.5, alpha=.8) +
  scale_fill_manual(values=CAT) + labs(x=NULL, y="offset bias (yr; recovered − true)") +
  theme_classic(base_size=12) + theme(legend.position="none")

g <- cowplot::plot_grid(pA, pB, pC, nrow=1, rel_widths=c(1,1.3,1),
                        labels=c("phylo-signal recovery","recovery accuracy","offset bias"), label_size=10, hjust=0)
ggsave(file.path(out,"fig_mc_recovery.png"), g, width=13, height=4.2, dpi=200, bg="white")
cat("wrote fig_mc_recovery.png\n")

# numeric summary
s <- r %>% group_by(regime) %>% summarise(across(c(offset_cor,neonate_cor,nu_cor,lambda_rec,offset_bias),
       list(med=~median(.x), lo=~quantile(.x,.1), hi=~quantile(.x,.9)), .names="{.col}_{.fn}"), n=n(), .groups="drop")
write_csv(s, file.path(out,"mc_recovery_summary.csv"))
print(as.data.frame(s %>% transmute(regime, n, offset_cor_med=round(offset_cor_med,2),
  neonate_cor_med=round(neonate_cor_med,2), nu_cor_med=round(nu_cor_med,2),
  lambda_rec_med=round(lambda_rec_med,2))), row.names=FALSE)
