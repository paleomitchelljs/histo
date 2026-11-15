// Model 13 (offset + phylogenetic-nu) + CATEGORICAL CHARACTER -> SHAPE test.
// Adds a per-species design matrix Xnu (K predictors) as fixed effects on lognu:
//
//   lognu_sp = mu_lognu + Xnu * gamma + nu_phylo[tip] + sd_lognu * z[species]
//
// gamma[k] is the association of histological character k with growth-curve
// SHAPE, *with phylogeny controlled* (nu_phylo) and ν-uncertainty propagated.
// gamma[k] whose posterior excludes 0 = a clear character->shape association.
// Designed for the feature model: drop in vascular/growthmark/tissue dummies as
// columns of Xnu. Everything else as in model 13.
data {
  int<lower=1> N; int<lower=1> S; int<lower=1> I; int<lower=1> D; int<lower=1> T;
  int<lower=0> K;                          // number of categorical predictors on nu
  vector<lower=0>[N] size;
  vector[N] age0;
  array[N] int<lower=1,upper=I> ind;
  array[N] int<lower=1,upper=D> dimr;
  array[I] int<lower=1,upper=S> sp_of_ind;
  array[S] int<lower=1,upper=T> tip_of_sp;
  matrix[T,T] Lphylo;
  matrix[S,K] Xnu;                         // species x predictors (centered dummies)
  array[I] int<lower=0,upper=1> efs;
  array[I] int<lower=0,upper=1> resorb;
  vector[I] logmax;
}
parameters {
  real mu_logA; real mu_logk; real mu_lognu; real mu_t0;
  real<lower=0> sd_logA_sp; real<lower=0> sd_logk; real<lower=0> sd_lognu; real<lower=0> sd_t0;
  real<lower=0> sd_logA_ind;
  real<lower=0> sigma_phy; vector[T] z_phy;
  vector[K] gamma;                         // character -> shape coefficients
  vector[S] zA; vector[S] zk; vector[S] znu; vector[S] zt0;
  vector[I] zAi;
  real<lower=0> tau_off;
  vector<lower=0>[I] miss_raw;
  vector<lower=0>[D] sigma;
}
transformed parameters {
  vector[T] nu_phylo = sigma_phy * (Lphylo * z_phy);
  vector[S] logA_sp  = mu_logA  + sd_logA_sp * zA;
  vector[S] logk_sp  = mu_logk  + sd_logk    * zk;
  vector[S] lognu_sp;
  vector[S] t0_sp    = mu_t0    + sd_t0      * zt0;
  vector[I] logA_i   = logA_sp[sp_of_ind] + sd_logA_ind * zAi;
  vector<lower=0>[I] miss;
  {
    vector[S] xb = (K > 0) ? (Xnu * gamma) : rep_vector(0.0, S);
    for (s in 1:S) lognu_sp[s] = mu_lognu + xb[s] + nu_phylo[tip_of_sp[s]] + sd_lognu * znu[s];
  }
  for (i in 1:I) miss[i] = (resorb[i] == 1 ? tau_off : 0.25) * miss_raw[i];
}
model {
  mu_logA ~ normal(5, 2); mu_logk ~ normal(-1, 1); mu_lognu ~ normal(0, 0.5); mu_t0 ~ normal(2, 3);
  sd_logA_sp ~ student_t(3, 0, 1); sd_logk ~ student_t(3, 0, 0.5);
  sd_lognu ~ student_t(3, 0, 0.5); sd_t0 ~ student_t(3, 0, 2); sd_logA_ind ~ student_t(3, 0, 0.5);
  sigma_phy ~ student_t(3, 0, 0.5); z_phy ~ std_normal();
  gamma ~ normal(0, 1);                    // weakly-informative on the association
  zA ~ std_normal(); zk ~ std_normal(); znu ~ std_normal(); zt0 ~ std_normal(); zAi ~ std_normal();
  tau_off ~ normal(0, 3); miss_raw ~ std_normal(); sigma ~ student_t(3, 0, 0.3);
  for (i in 1:I) {
    if (efs[i] == 1) logmax[i] ~ normal(logA_i[i], 0.10);
    else target += normal_lpdf(fmax(0.0, logmax[i] - logA_i[i]) | 0, 0.07);
  }
  {
    vector[N] mu;
    for (n in 1:N) {
      int s = sp_of_ind[ind[n]];
      real ta = age0[n] + miss[ind[n]];
      real nu = exp(lognu_sp[s]);
      mu[n] = logA_i[ind[n]] - (1.0 / nu) * log1p(nu * exp(-exp(logk_sp[s]) * (ta - t0_sp[s])));
    }
    size ~ lognormal(mu, sigma[dimr]);
  }
}
generated quantities {
  vector[K] gamma_out = gamma;             // character -> shape associations
  real lambda_nu = square(sigma_phy) / (square(sigma_phy) + square(sd_lognu));
}
