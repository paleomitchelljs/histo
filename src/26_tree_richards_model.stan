// Model 26: tree-nested multivariate growth model, FULL Richards family.
// The tree topology is the nesting structure: for every branch
//   child = parent + change,   change_b ~ MVN(0, brlen_b^kappa * Sigma_phy)
// jointly for K=4 species-level traits:
//   1 logA   log asymptotic size (per size-dimension intercept)
//   2 logr   log mid-growth rate: proportion of adult size gained per year
//            between 25% and 75% of the asymptote, r = 0.5 / T(25->75)
//   3 nu     Richards shape on a LINEAR scale, free to be negative
//            (1 logistic | 0 Gompertz | -1/3 von Bertalanffy | -1 monomolecular)
//   4 t50    age at 50% of the asymptote
// Curve: M(t) = A * (1 + nu * exp(x))^(-1/nu),  x = -k (t - t50) + g(0.5, nu),
//   g(f, nu) = log((f^-nu - 1) / nu),  k = 2 r (g(0.25) - g(0.75)).
// (r, t50) describe the same curves as (k, t0) but stay comparable across
// shapes and are far less correlated with nu, which the sampler needs.
// For nu < 0 the curve reaches size 0 at a finite age; earlier ages are
// penalized smoothly (soft floor on the bracket), not rejected.
//
//   dev[s, ] = (P * delta)[tip_of_sp[s], ] + species residual (sd_res)
//   logA_i   ~ normal(mu_logA[dimension] + dev[sp, 1], sd_logA_ind)
//
// SAMPLING FORM. The branch changes are integrated out: the species values are
// sampled directly with the covariance the tree implies,
//   vec(theta) ~ MVN(mean, Sigma_phy (x) C + diag(sd_res^2) (x) I),  C = P D P'
// and the branch changes delta (hence every node value) are drawn exactly from
// their conditional distribution in generated quantities. Same model as the
// branch-by-branch form; the data pin species values far more tightly than the
// prior, which makes the non-centred branch form nearly unsampleable.
//
// Missing inner rings: true age = recorded ring + miss_i, miss_i >= 0 with a
// half-normal(tau_off) prior, as before. For resorbed specimens the data only
// identify (t50 - miss), so the sampled parameter is d_off = t50_sp - m with
// miss = softplus(c m) / c (identity away from 0, smooth floor at 0; c = 4 for
// resorbed specimens, 16 otherwise) and the Jacobian added: same model, without
// the t50-vs-offset ridge. Prior scale: tau_off if resorbed, 0.25 yr if not.
// brlen = 1 everywhere -> topology-only; est_kappa = 1 with real branch lengths
// lets kappa in [0,1] run from topology-only to gradual change.
// sd_phy is on the TIP scale; lambda = sd_phy^2 / (sd_phy^2 + sd_res^2).
functions {
  real gfun(real f, real nu) {               // x at which relative size = f
    real a = -nu * log(f);
    if (abs(a) < 1e-8) return log(-log(f));
    return log(expm1(a) / nu);
  }
  real rich_log(real nu, real x) {           // log(M/A)
    real ex = exp(fmin(x, 50));
    real u = nu * ex;
    real w;
    if (abs(u) < 1e-4) return -ex * (1 - 0.5 * u + u * u / 3.0);   // Gompertz limit
    w = 1 + u;
    return -log(0.5 * (w + sqrt(w * w + 1e-12))) / nu;             // soft floor at w <= 0
  }
  // covariance of vec(theta) (column-major: trait blocks of S species)
  matrix build_sigma(matrix Vphy, matrix C, vector sd_res) {
    int S = rows(C); int K = rows(Vphy);
    matrix[S * K, S * K] Sig;
    for (k1 in 1:K) for (k2 in 1:K)
      Sig[((k1 - 1) * S + 1):(k1 * S), ((k2 - 1) * S + 1):(k2 * S)] = Vphy[k1, k2] * C;
    for (k in 1:K) for (s in 1:S) Sig[(k - 1) * S + s, (k - 1) * S + s] += square(sd_res[k]) + 1e-9;
    return Sig;
  }
}
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
  vector[K] sd_scale = [1.0, 0.5, 1.0, 3.0]';   // prior scales: logA, logr, nu, t50
  matrix[S,B] Ps = P[tip_of_sp, ];              // root->species path matrix
  vector[I] csoft;
  for (i in 1:I) csoft[i] = resorb[i] == 1 ? 4 : 16;
}
parameters {
  vector[D] mu_logA; real mu_logr; real mu_nu; real mu_t50;
  vector<lower=0>[K] sd_phy;               // tip-scale phylogenetic SD per trait
  vector<lower=0>[K] sd_res;               // species-specific (non-phylogenetic) SD
  cholesky_factor_corr[K] L_Om;            // evolutionary trait correlation
  matrix[S,K] theta;                       // species values: [,1] logA deviation, [,2] logr, [,3] nu, [,4] t50
  real<lower=0> sd_logA_ind;
  vector[I] logA_i;                        // centred: each specimen's rings pin its asymptote
  real<lower=0,upper=1> kappa_raw;
  real<lower=0> tau_off;
  vector[I] d_off;                         // t50_sp - m: apparent t50 on the recorded-ring clock
  vector<lower=0>[D] sigma;
}
transformed parameters {
  real kappa = est_kappa ? kappa_raw : 1.0;
  vector[B] bscale;                        // relative variance per branch
  vector[S] logr_sp = theta[, 2]; vector[S] nu_sp = theta[, 3]; vector[S] t50_sp = theta[, 4];
  vector[S] k_sp; vector[S] g50_sp;
  vector[I] m_off = t50_sp[sp_of_ind] - d_off;
  vector<lower=0>[I] miss = log1p_exp(csoft .* m_off) ./ csoft;
  {
    vector[B] bvar = exp(kappa * logbr);
    bscale = bvar / mean(P * bvar);
  }
  for (s in 1:S) {
    g50_sp[s] = gfun(0.5, nu_sp[s]);
    k_sp[s] = 2 * exp(logr_sp[s]) * (gfun(0.25, nu_sp[s]) - gfun(0.75, nu_sp[s]));
  }
}
model {
  mu_logA ~ normal(4, 2); mu_logr ~ normal(-2.5, 1); mu_nu ~ normal(0, 1); mu_t50 ~ normal(5, 5);
  sd_phy ~ normal(0, sd_scale);
  sd_res ~ normal(0, sd_scale);
  L_Om ~ lkj_corr_cholesky(2);
  {
    matrix[K,K] Vphy = multiply_lower_tri_self_transpose(diag_pre_multiply(sd_phy, L_Om));
    matrix[S,S] C = Ps * diag_post_multiply(Ps, bscale)';
    vector[S * K] mvec = append_row(rep_vector(0, S), append_row(rep_vector(mu_logr, S),
                           append_row(rep_vector(mu_nu, S), rep_vector(mu_t50, S))));
    to_vector(theta) ~ multi_normal(mvec, build_sigma(Vphy, C, sd_res));
  }
  sd_logA_ind ~ normal(0, 0.5);
  logA_i ~ normal(mu_logA[dim_of_ind] + theta[sp_of_ind, 1], sd_logA_ind);
  kappa_raw ~ beta(1, 1);
  tau_off ~ normal(0, 3); sigma ~ student_t(3, 0, 0.3);
  for (i in 1:I) {
    target += normal_lpdf(miss[i] | 0, resorb[i] == 1 ? tau_off : 0.25) + log_inv_logit(csoft[i] * m_off[i]);
    if (efs[i] == 1) logmax[i] ~ normal(logA_i[i], 0.10);
    else target += normal_lpdf(fmax(0.0, logmax[i] - logA_i[i]) | 0, 0.07);
  }
  {
    vector[N] mu;
    for (n in 1:N) {
      int s = sp_of_ind[ind[n]];
      real ta = age0[n] + miss[ind[n]];
      mu[n] = logA_i[ind[n]] + rich_log(nu_sp[s], -k_sp[s] * (ta - t50_sp[s]) + g50_sp[s]);
    }
    size ~ lognormal(mu, sigma[dimr]);
  }
}
generated quantities {
  vector[K] lambda = square(sd_phy) ./ (square(sd_phy) + square(sd_res));
  corr_matrix[K] Omega = multiply_lower_tri_self_transpose(L_Om);
  vector[I] miss_yr = miss;
  matrix[S,K] dev;                         // species deviation from the root state
  matrix[B,K] delta;                       // per-branch change: exact conditional draw given theta
  {
    matrix[K,K] Lphy = diag_pre_multiply(sd_phy, L_Om);
    matrix[K,K] Vphy = multiply_lower_tri_self_transpose(Lphy);
    matrix[S,S] C = Ps * diag_post_multiply(Ps, bscale)';
    matrix[B,K] dstar; matrix[S,K] tstar; matrix[S,K] W;
    dev = theta;
    dev[, 2] -= rep_vector(mu_logr, S); dev[, 3] -= rep_vector(mu_nu, S); dev[, 4] -= rep_vector(mu_t50, S);
    for (b in 1:B) for (k in 1:K) dstar[b, k] = normal_rng(0, 1);
    dstar = diag_pre_multiply(sqrt(bscale), dstar) * Lphy';
    tstar = Ps * dstar;
    for (s in 1:S) for (k in 1:K) tstar[s, k] += normal_rng(0, sd_res[k]);
    W = to_matrix(build_sigma(Vphy, C, sd_res) \ to_vector(dev - tstar), S, K);
    delta = dstar + diag_pre_multiply(bscale, Ps' * W) * Vphy;
  }
}
