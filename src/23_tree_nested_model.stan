// Model 23: tree-nested multivariate growth model.
// The tree topology IS the nesting structure: every node value is
//   child = parent + change,   change_b ~ MVN(0, brlen_b^kappa * Sigma_phy)
// for the K=4 species-level growth traits (logA, logk, lognu, t0) jointly.
// Tip value = sum of changes on the root->tip path (P * delta). Species add a
// non-phylogenetic residual; individuals add a residual on logA.
//
//   dev[s, ]  = (P * delta)[tip_of_sp[s], ] + sd_res .* z_res[s, ]
//   logA_i    = mu_logA[dimension] + dev[sp, 1] + sd_logA_ind * z_i
//   logk_sp   = mu_logk  + dev[, 2];  lognu_sp = mu_lognu + dev[, 3];  t0_sp = mu_t0 + dev[, 4]
//
// brlen = 1 for every branch -> topology-only (speciational) model. With real
// branch lengths and est_kappa = 1, kappa in [0,1] interpolates between
// topology-only (0) and gradual Brownian change (1).
// sd_phy is on the TIP scale: mean phylogenetic variance across tips = sd_phy^2,
// so lambda = sd_phy^2 / (sd_phy^2 + sd_res^2) per trait.
// mu_logA is per size-dimension (mm circumference, % length, rescaled mass are
// not commensurable); the tree therefore sees logA only within a dimension.
// Offset / EFS / Richards likelihood: as in model 13.
data {
  int<lower=1> N; int<lower=1> S; int<lower=1> I; int<lower=1> D;
  int<lower=1> T; int<lower=1> B;
  vector<lower=0>[N] size;
  vector[N] age0;
  array[N] int<lower=1,upper=I> ind;
  array[N] int<lower=1,upper=D> dimr;
  array[I] int<lower=1,upper=D> dim_of_ind;
  array[I] int<lower=1,upper=S> sp_of_ind;
  array[S] int<lower=1,upper=T> tip_of_sp;
  matrix[T,B] P;                           // root->tip path matrix (0/1)
  vector<lower=0>[B] brlen;                // all 1 = topology only
  int<lower=0,upper=1> est_kappa;
  array[I] int<lower=0,upper=1> efs;
  array[I] int<lower=0,upper=1> resorb;
  vector[I] logmax;
}
transformed data {
  int K = 4;
  vector[B] logbr = log(brlen);
  vector[K] sd_scale = [1.0, 0.5, 0.5, 2.0]';   // prior scales: logA, logk, lognu, t0
}
parameters {
  vector[D] mu_logA; real mu_logk; real mu_lognu; real mu_t0;
  vector<lower=0>[K] sd_phy;               // tip-scale phylogenetic SD per trait
  vector<lower=0>[K] sd_res;               // species-specific (non-phylogenetic) SD
  cholesky_factor_corr[K] L_Om;            // evolutionary trait correlation
  matrix[B,K] z_delta;
  matrix[S,K] z_res;
  real<lower=0> sd_logA_ind; vector[I] zAi;
  real<lower=0,upper=1> kappa_raw;
  real<lower=0> tau_off;
  vector<lower=0>[I] miss_raw;
  vector<lower=0>[D] sigma;
}
transformed parameters {
  real kappa = est_kappa ? kappa_raw : 1.0;
  matrix[B,K] delta;                       // per-branch change, per trait
  matrix[T,K] phy;                         // tip value relative to root
  matrix[S,K] dev;                         // species deviation from root
  vector[S] logk_sp; vector[S] lognu_sp; vector[S] t0_sp;
  vector[I] logA_i;
  vector<lower=0>[I] miss;
  {
    vector[B] bvar = exp(kappa * logbr);
    real dbar = mean(P * bvar);            // mean root->tip variance units
    delta = diag_pre_multiply(sqrt(bvar / dbar), z_delta) * diag_pre_multiply(sd_phy, L_Om)';
  }
  phy = P * delta;
  dev = phy[tip_of_sp, ] + diag_post_multiply(z_res, sd_res);
  logk_sp  = mu_logk  + dev[, 2];
  lognu_sp = mu_lognu + dev[, 3];
  t0_sp    = mu_t0    + dev[, 4];
  logA_i   = mu_logA[dim_of_ind] + dev[sp_of_ind, 1] + sd_logA_ind * zAi;
  for (i in 1:I) miss[i] = (resorb[i] == 1 ? tau_off : 0.25) * miss_raw[i];
}
model {
  mu_logA ~ normal(4, 2); mu_logk ~ normal(-1, 1); mu_lognu ~ normal(0, 0.5); mu_t0 ~ normal(2, 3);
  sd_phy ~ student_t(3, 0, sd_scale);
  sd_res ~ student_t(3, 0, sd_scale);
  L_Om ~ lkj_corr_cholesky(2);
  to_vector(z_delta) ~ std_normal(); to_vector(z_res) ~ std_normal();
  sd_logA_ind ~ student_t(3, 0, 0.5); zAi ~ std_normal();
  kappa_raw ~ beta(1, 1);
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
  vector[K] lambda = square(sd_phy) ./ (square(sd_phy) + square(sd_res));
  corr_matrix[K] Omega = multiply_lower_tri_self_transpose(L_Om);
  vector[I] miss_yr = miss;
}
