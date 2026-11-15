// Branch-shift model for the shape parameter nu (plan Option B).
// Every branch of the tree gets its own increment to lognu, regularized by a
// REGULARIZED HORSESHOE so that most increments shrink to ~0 and the few real
// shifts stand out. tip-level shape deviation = sum of increments along the
// root->tip path (P * delta). Branches whose delta posterior excludes 0 are the
// "shifts" -- where growth-curve SHAPE changes on the tree.
//
//   lognu_sp = mu_lognu + (P * delta)[tip] + sd_lognu * z[species]
//
// Offset / EFS / per-specimen A / k / t0 / per-dimension sigma: as in model 13.
data {
  int<lower=1> N; int<lower=1> S; int<lower=1> I; int<lower=1> D;
  int<lower=1> T;                          // tips
  int<lower=1> B;                          // branches (tree edges)
  vector<lower=0>[N] size;
  vector[N] age0;
  array[N] int<lower=1,upper=I> ind;
  array[N] int<lower=1,upper=D> dimr;
  array[I] int<lower=1,upper=S> sp_of_ind;
  array[S] int<lower=1,upper=T> tip_of_sp;
  matrix[T,B] P;                           // root->tip path matrix (0/1)
  array[I] int<lower=0,upper=1> efs;
  array[I] int<lower=0,upper=1> resorb;
  vector[I] logmax;
  real<lower=0> slab_scale;                // slab SD for a real shift (lognu units)
  real<lower=0> tau0;                       // global shrinkage scale (sparsity)
}
parameters {
  real mu_logA; real mu_logk; real mu_lognu; real mu_t0;
  real<lower=0> sd_logA_sp; real<lower=0> sd_logk; real<lower=0> sd_lognu; real<lower=0> sd_t0;
  real<lower=0> sd_logA_ind;
  // regularized horseshoe on branch increments
  vector[B] z_delta;
  vector<lower=0>[B] lambda;
  real<lower=0> tau;
  real<lower=0> caux;
  vector[S] zA; vector[S] zk; vector[S] znu; vector[S] zt0;
  vector[I] zAi;
  real<lower=0> tau_off;
  vector<lower=0>[I] miss_raw;
  vector<lower=0>[D] sigma;
}
transformed parameters {
  real c2 = square(slab_scale) * caux;                                   // slab variance
  vector<lower=0>[B] lambda_t = sqrt(c2 * square(lambda) ./ (c2 + square(tau) * square(lambda)));
  vector[B] delta = z_delta .* lambda_t * tau;                           // sparse branch shifts
  vector[T] tip_dev = P * delta;                                         // accumulated along paths
  vector[S] logA_sp  = mu_logA  + sd_logA_sp * zA;
  vector[S] logk_sp  = mu_logk  + sd_logk    * zk;
  vector[S] lognu_sp;
  vector[S] t0_sp    = mu_t0    + sd_t0      * zt0;
  vector[I] logA_i   = logA_sp[sp_of_ind] + sd_logA_ind * zAi;
  vector<lower=0>[I] miss;
  for (s in 1:S) lognu_sp[s] = mu_lognu + tip_dev[tip_of_sp[s]] + sd_lognu * znu[s];
  for (i in 1:I) miss[i] = (resorb[i] == 1 ? tau_off : 0.25) * miss_raw[i];
}
model {
  mu_logA ~ normal(5, 2); mu_logk ~ normal(-1, 1); mu_lognu ~ normal(0, 0.5); mu_t0 ~ normal(2, 3);
  sd_logA_sp ~ student_t(3, 0, 1); sd_logk ~ student_t(3, 0, 0.5);
  sd_lognu ~ student_t(3, 0, 0.3); sd_t0 ~ student_t(3, 0, 2); sd_logA_ind ~ student_t(3, 0, 0.5);
  z_delta ~ std_normal();
  lambda ~ cauchy(0, 1);
  tau ~ cauchy(0, tau0);                    // global shrinkage -> sparsity
  caux ~ inv_gamma(2, 2);                    // slab (regularizes large shifts)
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
  vector[B] delta_out = delta;              // per-branch shape shift (lognu units)
  vector[I] miss_yr = miss;
}
