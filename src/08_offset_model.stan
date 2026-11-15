// Hierarchical generalized-logistic (Richards) growth with an explicit
// per-specimen missing-LAG miss and categorical histology information.
//
//   true_age = recorded_ring_index + miss[specimen],  miss >= 0
//   miss prior scale depends on resorption category (absent -> ~0)
//   EFS present -> outermost ring is at the asymptote (anchors A)
//   size ~ lognormal(logRichards, sigma[dimension])
//
// Shape (lognu), rate (logk), location (t0) are per-species (partial pooling);
// scale (logA) is per-specimen nested in species. This is the principled
// replacement for letting an unconstrained t0 silently absorb resorption.
data {
  int<lower=1> N;                         // observations
  int<lower=1> S;                         // species
  int<lower=1> I;                         // specimens (individuals)
  int<lower=1> D;                         // size dimensions (circumference / % length)
  vector<lower=0>[N] size;
  vector[N] age0;                         // recorded ring index (true age for perinate-anchored taxa)
  array[N] int<lower=1,upper=I> ind;
  array[N] int<lower=1,upper=D> dimr;
  array[I] int<lower=1,upper=S> sp_of_ind;
  array[I] int<lower=0,upper=1> efs;      // EFS present?
  array[I] int<lower=0,upper=1> resorb;   // 1 = inner rings resorbed (miss > 0); 0 = absent
  vector[I] logmax;                       // log of max observed size per specimen
}
parameters {
  real mu_logA; real mu_logk; real mu_lognu; real mu_t0;
  real<lower=0> sd_logA_sp; real<lower=0> sd_logk; real<lower=0> sd_lognu; real<lower=0> sd_t0;
  real<lower=0> sd_logA_ind;
  vector[S] zA; vector[S] zk; vector[S] znu; vector[S] zt0;
  vector[I] zAi;
  real<lower=0> tau_off;                  // missing-ring scale for resorbed specimens
  vector<lower=0>[I] miss_raw;          // half-normal raw
  vector<lower=0>[D] sigma;
}
transformed parameters {
  vector[S] logA_sp  = mu_logA  + sd_logA_sp * zA;
  vector[S] logk_sp  = mu_logk  + sd_logk    * zk;
  vector[S] lognu_sp = mu_lognu + sd_lognu   * znu;
  vector[S] t0_sp    = mu_t0    + sd_t0      * zt0;
  vector[I] logA_i   = logA_sp[sp_of_ind] + sd_logA_ind * zAi;
  vector<lower=0>[I] miss;
  for (i in 1:I) miss[i] = (resorb[i] == 1 ? tau_off : 0.25) * miss_raw[i];
}
model {
  mu_logA ~ normal(5, 2); mu_logk ~ normal(-1, 1); mu_lognu ~ normal(0, 0.5); mu_t0 ~ normal(2, 3);
  sd_logA_sp ~ student_t(3, 0, 1); sd_logk ~ student_t(3, 0, 0.5);
  sd_lognu ~ student_t(3, 0, 0.5); sd_t0 ~ student_t(3, 0, 2); sd_logA_ind ~ student_t(3, 0, 0.5);
  zA ~ std_normal(); zk ~ std_normal(); znu ~ std_normal(); zt0 ~ std_normal(); zAi ~ std_normal();
  tau_off ~ normal(0, 3);                 // half-normal (lower=0)
  miss_raw ~ std_normal();              // half-normal via lower=0
  sigma ~ student_t(3, 0, 0.3);

  // categorical histology -> asymptote constraint
  for (i in 1:I) {
    if (efs[i] == 1) logmax[i] ~ normal(logA_i[i], 0.10);          // growth ceased: A ~ max ring
    else target += normal_lpdf(fmax(0.0, logmax[i] - logA_i[i]) | 0, 0.07); // still growing: A >= max
  }

  // likelihood
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
  vector[I] miss_yr = miss;           // estimated missing inner rings per specimen
}
