#!/usr/bin/env Rscript
# Monte Carlo parameter-recovery study: repeat each scenario n_reps times, fit the
# ACTUAL model (13, minimal MCMC ~28s/fit), and record recovery of curve shape nu,
# the missing-LAG offset, neonate mass, and the recovered phylo signal lambda_nu.
# Results appended incrementally to outputs/mc_recovery_results.csv.
# Usage: Rscript 20_recovery_montecarlo.R <n_reps=100>
suppressMessages({ library(cmdstanr); library(ape); library(dplyr); library(readr) })  # no phytools (background-safe)
.args <- commandArgs(trailingOnly=TRUE)
nreps <- as.integer(.args[1]); if (is.na(nreps)) nreps <- 100
chunk <- as.integer(.args[2]); if (is.na(chunk)) chunk <- 0
set.seed(20260609 + 1000*chunk)            # distinct reps per chunk (for foreground chunking)
root  <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)); Sys.getenv("HISTO_ROOT", if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else getwd()) })
clean <- file.path(root,"data_clean"); out <- file.path(root,"outputs")
resfile <- file.path(out,"mc_recovery_results.csv")

tree <- compute.brlen(read.tree(file.path(clean,"dino_tree.nwk")), method="Grafen")
tips <- tree$tip.label; T <- length(tips)
Rc <- vcv(tree, corr=TRUE); L <- t(chol(Rc)); tip_of_sp <- match(tips, rownames(Rc))
LC <- t(chol(vcv(tree)))                                  # BM covariance Cholesky (replaces fastBM)
desc_tips <- function(node){                              # tips descending from an internal node (replaces getDescendants)
  kids <- tree$edge[tree$edge[,1]==node,2]; out <- integer(0)
  while (length(kids)>0){ k <- kids[1]; kids <- kids[-1]
    if (k<=T) out <- c(out,k) else kids <- c(kids, tree$edge[tree$edge[,1]==k,2]) }
  out }
Ctree <- vcv(tree); dC <- diag(Ctree)                    # Pagel's lambda MLE (replaces phylosig)
pagel_lambda <- function(x){
  ll <- function(lam){ Cl <- Ctree*lam; diag(Cl) <- dC; Ci <- solve(Cl); n <- length(x)
    mu <- sum(Ci %*% x)/sum(Ci); r <- x-mu; s2 <- as.numeric(t(r)%*%Ci%*%r)/n
    -0.5*(n*log(2*pi*s2) + as.numeric(determinant(Cl,log=TRUE)$modulus) + n) }
  optimize(ll, c(0,1), maximum=TRUE)$maximum }
N_PER<-2; SIGMA<-0.05; MU_LOGNU<- -0.4; SD_NU<-1.0; MU_LOGK<- -1.1; SD_LOGK<-0.35
MU_T0<-2.0; SD_T0<-1.2; MU_LOGA<-5.6; SD_LOGA_SP<-0.45; SD_LOGA_IND<-0.20
rich <- function(t,logA,nu,k,t0) logA-(1/nu)*log1p(nu*exp(-k*(t-t0)))

gen_lognu <- function(scn){
  if (scn=="iid")  return(setNames(rnorm(T,MU_LOGNU,SD_NU),tips))
  if (scn=="bm")   { raw<-as.numeric(LC %*% rnorm(T)); x<-MU_LOGNU+SD_NU*(raw-mean(raw))/sd(raw); return(setNames(x,tips)) }
  if (scn=="shift"){ v<-setNames(rep(MU_LOGNU,T),tips); inodes<-(T+1):(T+tree$Nnode)
    for (jmp in c(+1.8,-1.6)){ nd<-sample(inodes[sapply(inodes,function(n){dt<-length(desc_tips(n)); dt>=3&dt<=8})],1)
      dt<-desc_tips(nd); v[tips[dt]]<-v[tips[dt]]+jmp }
    return(v) } }

gen_data <- function(scn){
  lognu<-gen_lognu(scn); logk<-setNames(rnorm(T,MU_LOGK,SD_LOGK),tips)
  t0<-setNames(rnorm(T,MU_T0,SD_T0),tips); logAs<-setNames(rnorm(T,MU_LOGA,SD_LOGA_SP),tips)
  rows<-list(); tr<-list(); ii<-0
  for (s in tips){ nu<-exp(lognu[s]);k<-exp(logk[s]);tt<-t0[s]
    for (r in 1:N_PER){ ii<-ii+1; spec<-paste0(s,"_",r); logAi<-logAs[s]+rnorm(1,0,SD_LOGA_IND)
      death<-sample(7:22,1); resorbed<-rbinom(1,1,0.6); offset<-if(resorbed) min(sample(1:7,1),death-3) else 0
      efs<-as.integer(death>=15 & runif(1)<0.6); neo<-exp(rich(0,logAi,nu,k,tt))
      ta<-if(offset==0) c(0,1:death) else (offset+1):death
      sz<-exp(sapply(ta,rich,logA=logAi,nu=nu,k=k,t0=tt)+rnorm(length(ta),0,SIGMA))
      rows[[ii]]<-data.frame(species=s,specimen=spec,age=ta-offset,size=sz,
        resorption=ifelse(resorbed,"present","absent"),efs_present=efs==1)
      tr[[ii]]<-data.frame(specimen=spec,offset_true=offset,neonate_true=neo) } }
  list(data=bind_rows(rows), truth_sp=data.frame(species=tips,nu=exp(as.numeric(lognu[tips]))),
       truth_ind=bind_rows(tr), lam_true=pagel_lambda(as.numeric(lognu[tips])))
}

cat("MARK: libs+tree ready, compiling model\n", file=stderr())
mod <- cmdstan_model(file.path(root,"src/13_offset_phylo_model.stan"))
cat("MARK: model compiled\n", file=stderr())
fit_rep <- function(g){
  dat<-g$data; dat$species<-factor(dat$species,levels=tips); dat$specimen<-factor(dat$specimen)
  it<-dat%>%distinct(specimen,species,resorption,efs_present)%>%
    left_join(dat%>%group_by(specimen)%>%summarise(logmax=log(max(size)),.groups="drop"),by="specimen")%>%arrange(specimen)
  it$sp_of_ind<-as.integer(it$species)
  sd<-list(N=nrow(dat),S=length(tips),I=nlevels(dat$specimen),D=1L,T=T,size=dat$size,age0=dat$age,
    ind=as.integer(dat$specimen),dimr=rep(1L,nrow(dat)),sp_of_ind=it$sp_of_ind,tip_of_sp=tip_of_sp,
    Lphylo=L,efs=as.integer(it$efs_present),resorb=as.integer(it$resorption!="absent"),logmax=it$logmax)
  f<-mod$sample(data=sd,chains=1,iter_warmup=200,iter_sampling=100,adapt_delta=0.9,max_treedepth=10,
                refresh=0,seed=sample.int(1e6,1),init=0.4)
  d<-f$draws(format="df"); gm<-function(p,i) median(d[[paste0(p,"[",i,"]")]])
  nu_rec<-sapply(seq_along(tips),gm,p="lognu_sp")
  off_rec<-sapply(seq_len(nrow(it)),gm,p="miss_yr"); sp_idx<-it$sp_of_ind
  neo_rec<-sapply(seq_len(nrow(it)),function(i){s<-sp_idx[i];la<-d[[paste0("logA_i[",i,"]")]];nu<-exp(d[[paste0("lognu_sp[",s,"]")]]);k<-exp(d[[paste0("logk_sp[",s,"]")]]);t0<-d[[paste0("t0_sp[",s,"]")]];median(exp(la-(1/nu)*log1p(nu*exp(-k*(0-t0)))))})
  tn<-log(g$truth_sp$nu); ti<-g$truth_ind[match(as.character(it$specimen),g$truth_ind$specimen),]
  rmse<-function(a,b) sqrt(mean((a-b)^2))
  data.frame(nu_cor=cor(tn,nu_rec), nu_rmse=rmse(tn,nu_rec),
    offset_cor=cor(ti$offset_true,off_rec), offset_rmse=rmse(ti$offset_true,off_rec),
    offset_bias=mean(off_rec-ti$offset_true), neonate_cor=cor(log(ti$neonate_true),log(neo_rec)),
    neonate_rmse=rmse(log(ti$neonate_true),log(neo_rec)), lambda_rec=median(d[["lambda_nu"]]),
    lam_true=g$lam_true, ndiv=sum(f$diagnostic_summary()$num_divergent))
}

t0all<-Sys.time()
for (scn in c("iid","bm","shift")) for (rep in 1:nreps) {
  cat(sprintf("MARK: start %s rep %d\n",scn,rep), file=stderr())
  r <- tryCatch(fit_rep(gen_data(scn)), error=function(e){cat("rep err",scn,rep,conditionMessage(e),"\n");NULL})
  if (is.null(r)) next
  r$regime<-scn; r$rep<-rep
  write_csv(r, resfile, append=file.exists(resfile))
  if (rep %% 10 == 0) cat(sprintf("[%s] rep %d/%d done | elapsed %.0f min\n",scn,rep,nreps,as.numeric(Sys.time()-t0all,units="mins")))
}
cat("DONE MC recovery (", nreps, "reps x 3 regimes)\n")
