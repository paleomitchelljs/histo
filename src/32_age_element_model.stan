// Model 32: model 26 (tree-nested full Richards, marginal tree form) with
//   R1  one latent AGE per individual, elements nested in individuals
//   R2  a count-and-size module on the same curve
//   R3  three-state EFS and a one-sided asymptote pin
// Curve, species traits (logA, logr, nu, t50), tree covariance and the generated
// branch changes are unchanged from src/26_tree_richards_model.stan.
//
// INDIVIDUAL i has one age Tage[i]: the age at its outermost growth mark for
// count-based animals, or its age at death where the age is known/estimated.
//   age_type 0: Tage >= Lmax (largest count among its elements); every element's
//               lost marks = Tage - (marks it preserves) get a half-normal prior.
//   age_type 1: Tage ~ normal(age_est, age_sd)  (known age, or the source's estimate)
// Tage = Lmax + softplus(c m)/c. Two parametrizations of m (data flag tpar):
//   tpar 0  m = t50_sp - t_raw   resorbed ring series: the marks identify t50 - (lost
//           marks), so this removes the t50-vs-offset ridge (as in model 26)
//   tpar 1  m = t_loc + t_scale * t_raw   everything else (known or estimated ages,
//           count-and-size animals): the age is set by its prior or by the count, and
//           tying it to t50 would create a ridge instead of removing one.
//
// ELEMENT e (bone x individual) has its own asymptote
//   logA_e = mu_logA[dimension] + theta[sp,1] + beta_el[element type] + a_i + b_e
// a_i (individual) and b_e (element) split sd_logA_ind^2 by rho_ind; a_i is only
// sampled for individuals with more than one sized element.
//  ring element: sized marks at ages age0 + miss_e, miss_e = Tage - last_e, so that
//    marks preserved + marks lost = the individual's age for every bone of one animal.
//    If the outer cortex is eroded, outer marks may be missing too:
//    miss_e = f_e (Tage - last_e), the rest is outer loss (half-normal(tau_off)).
//  count-and-size element: optional count N_e (lost = Tage - N_e, same prior) and
//    optional size at death S_e, one point on the curve at age Tage + t_add with the
//    element effect b_e integrated out:
//    log S_e ~ normal(mu + rich_log(Tage + t_add), sqrt(sd_el^2 / nrep + sigma_cs^2 + size_logsd^2))
//    (nrep > 1: the mean of nrep animals, e.g. an age class).
//
// EFS, three states (0 absent, 1 present, 2 unknown). An EFS is a Bernoulli
// observation whose probability rises with the fraction of the asymptote attained:
//   P(EFS) = inv_logit((rel - efs_thr) / efs_slope)
// rel = (largest mark) / A_e for ring elements, the curve value at death for
// count-and-size elements. Present favours rel above efs_thr but does not pin the
// outer mark to A (growth continues after an EFS); absent favours rel below it;
// unknown contributes nothing. All ring elements keep the floor A_e >= largest mark.
functions {
  real gfun(real f, real nu) {
    real a = -nu * log(f);
    if (abs(a) < 1e-8) return log(-log(f));
    return log(expm1(a) / nu);
  }
  real rich_log(real nu, real x) {
    real ex = exp(fmin(x, 50));
    real u = nu * ex;
    real w;
    if (abs(u) < 1e-4) return -ex * (1 - 0.5 * u + u * u / 3.0);
    w = 1 + u;
    return -log(0.5 * (w + sqrt(w * w + 1e-12))) / nu;
  }
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
  int<lower=1> S; int<lower=1> D; int<lower=1> T; int<lower=1> B; int<lower=0> L;
  array[S] int<lower=1,upper=T> tip_of_sp;
  matrix[T,B] P; vector<lower=0>[B] brlen; int<lower=0,upper=1> est_kappa;
  // individuals
  int<lower=1> I;
  array[I] int<lower=1,upper=S> sp_of_ind;
  array[I] int<lower=0,upper=1> age_type;
  vector<lower=0>[I] Lmax; vector[I] age_est; vector<lower=0>[I] age_sd; vector<lower=0>[I] csoft;
  array[I] int<lower=0,upper=1> tpar; vector[I] t_loc; vector<lower=0>[I] t_scale;
  int<lower=0> IM; array[I] int<lower=0,upper=IM> im_of_ind;      // index among multi-element individuals (0 = none)
  // ring elements
  int<lower=0> ER;
  array[ER] int<lower=1,upper=I> r_ind; array[ER] int<lower=1,upper=D> r_dim; array[ER] int<lower=0,upper=L> r_elt;
  int<lower=1> DR; array[ER] int<lower=1,upper=DR> r_sdim;        // residual SD index: dimensions that have ring series
  array[ER] int<lower=0,upper=2> r_efs; array[ER] int<lower=0,upper=1> r_resorb;
  int<lower=0> NERO; array[ER] int<lower=0,upper=NERO> r_ero;     // index among eroded ring elements (0 = not eroded)
  vector[ER] r_last; vector[ER] r_logmax;
  int<lower=0> N; vector<lower=0>[N] size; vector[N] age0; array[N] int<lower=1,upper=ER> el;
  // count-and-size elements
  int<lower=0> EC;
  array[EC] int<lower=1,upper=I> c_ind; array[EC] int<lower=1,upper=D> c_dim; array[EC] int<lower=0,upper=L> c_elt;
  array[EC] int<lower=0,upper=2> c_efs; array[EC] int<lower=0,upper=1> c_resorb;
  array[EC] int<lower=0,upper=1> c_has_count; vector[EC] c_count;
  array[EC] int<lower=0,upper=1> c_has_size; vector[EC] c_logsize; vector<lower=0>[EC] c_logsd;
  vector<lower=1>[EC] c_nrep; vector<lower=0>[EC] c_tadd;
  array[S] int<lower=1,upper=D> main_dim;  // each species' best-sampled size dimension
  real<lower=0,upper=1> efs_thr; real<lower=0> efs_slope;
  real<lower=0> tau_fresh;                 // lost-mark scale where no resorption is reported
}
transformed data {
  int K = 4;
  vector[B] logbr = log(brlen);
  vector[K] sd_scale = [1.0, 0.5, 1.0, 3.0]';
  matrix[S,B] Ps = P[tip_of_sp, ];
}
parameters {
  vector[D] mu_logA; real mu_logr; real mu_nu; real mu_t50;
  vector<lower=0>[K] sd_phy; vector<lower=0>[K] sd_res;
  cholesky_factor_corr[K] L_Om;
  // theta_raw[,1] = mu_logA[main_dim] + logA deviation: the species' absolute log asymptote
  // in its main dimension, which is what the data pin. Sampling it (instead of the
  // deviation) makes the weakly identified dimension intercepts mu_logA their own axis.
  // Same model; the map is a shear with unit Jacobian.
  matrix[S,K] theta_raw;
  vector[L] beta_el;                       // element-type effect on log asymptote (femur and body measures = 0)
  real<lower=0> sd_logA_ind;               // SD of element asymptotes within a species (individual + element)
  real<lower=0,upper=1> rho_ind;           // share of that variance common to the bones of one animal
  vector[IM] z_ind;
  vector[ER] logA_e;
  real<lower=0,upper=1> kappa_raw;
  real<lower=0> tau_off;
  vector[I] t_raw;
  vector<lower=0,upper=1>[NERO] f_ero;
  vector<lower=0>[DR] sigma;
  real<lower=0> sigma_cs;
}
transformed parameters {
  real kappa = est_kappa ? kappa_raw : 1.0;
  matrix[S,K] theta = theta_raw;           // [,1] logA deviation, [,2] logr, [,3] nu, [,4] t50
  vector[B] bscale;
  vector[S] logr_sp = theta_raw[, 2]; vector[S] nu_sp = theta_raw[, 3]; vector[S] t50_sp = theta_raw[, 4];
  vector[S] k_sp; vector[S] g50_sp;
  vector[I] m_off;
  vector<lower=0>[I] Tage;
  vector[I] a_ind;                         // individual effect on log asymptote
  vector<lower=0>[I] sd_el;                // element SD given the individual
  vector[ER] miss;                         // lost INNER marks of ring elements
  theta[, 1] -= mu_logA[main_dim];
  {
    vector[B] bvar = exp(kappa * logbr);
    bscale = bvar / mean(P * bvar);
  }
  for (s in 1:S) {
    g50_sp[s] = gfun(0.5, nu_sp[s]);
    k_sp[s] = 2 * exp(logr_sp[s]) * (gfun(0.25, nu_sp[s]) - gfun(0.75, nu_sp[s]));
  }
  for (i in 1:I) {
    m_off[i] = tpar[i] == 0 ? t50_sp[sp_of_ind[i]] - t_raw[i] : t_loc[i] + t_scale[i] * t_raw[i];
    Tage[i] = Lmax[i] + log1p_exp(csoft[i] * m_off[i]) / csoft[i];
    if (im_of_ind[i] > 0) { a_ind[i] = sqrt(rho_ind) * sd_logA_ind * z_ind[im_of_ind[i]]; sd_el[i] = sqrt(1 - rho_ind) * sd_logA_ind; }
    else { a_ind[i] = 0; sd_el[i] = sd_logA_ind; }
  }
  for (e in 1:ER) {
    real R = Tage[r_ind[e]] - r_last[e];
    miss[e] = r_ero[e] > 0 ? f_ero[r_ero[e]] * R : R;
  }
}
model {
  mu_logA ~ normal(4, 2); mu_logr ~ normal(-2.5, 1); mu_nu ~ normal(0, 1); mu_t50 ~ normal(5, 5);
  sd_phy ~ normal(0, sd_scale); sd_res ~ normal(0, sd_scale);
  L_Om ~ lkj_corr_cholesky(2);
  {
    matrix[K,K] Vphy = multiply_lower_tri_self_transpose(diag_pre_multiply(sd_phy, L_Om));
    matrix[S,S] C = Ps * diag_post_multiply(Ps, bscale)';
    vector[S * K] mvec = append_row(rep_vector(0, S), append_row(rep_vector(mu_logr, S),
                           append_row(rep_vector(mu_nu, S), rep_vector(mu_t50, S))));
    to_vector(theta) ~ multi_normal(mvec, build_sigma(Vphy, C, sd_res));
  }
  beta_el ~ normal(0, 1);
  sd_logA_ind ~ normal(0, 0.5); rho_ind ~ beta(2, 2); z_ind ~ std_normal();
  kappa_raw ~ beta(1, 1);
  tau_off ~ normal(0, 3); sigma ~ student_t(3, 0, 0.3); sigma_cs ~ normal(0, 0.15);

  // individual ages
  for (i in 1:I) {
    target += log_inv_logit(csoft[i] * m_off[i]);                 // Jacobian d Tage / d t_raw (up to a constant)
    if (age_type[i] == 1) Tage[i] ~ normal(age_est[i], age_sd[i]);
  }
  // ring elements
  for (e in 1:ER) {
    int i = r_ind[e]; int s = sp_of_ind[i];
    real bel = r_elt[e] > 0 ? beta_el[r_elt[e]] : 0;
    real gap = r_logmax[e] - logA_e[e];
    real tau = r_resorb[e] == 1 ? tau_off : tau_fresh;
    logA_e[e] ~ normal(mu_logA[r_dim[e]] + theta[s, 1] + bel + a_ind[i], sd_el[i]);
    target += normal_lpdf(miss[e] | 0, tau);
    if (r_ero[e] > 0) {
      real R = Tage[i] - r_last[e];
      target += normal_lpdf(R - miss[e] | 0, tau_off) + log(R + 1e-8);   // outer loss + Jacobian of f -> miss
    }
    target += normal_lpdf(fmax(0.0, gap) | 0, 0.07);              // A_e >= largest mark
    if (r_efs[e] == 1)      target += log_inv_logit((exp(fmin(0.0, gap)) - efs_thr) / efs_slope);
    else if (r_efs[e] == 0) target += log1m_inv_logit((exp(fmin(0.0, gap)) - efs_thr) / efs_slope);
  }
  {
    vector[N] mu;
    for (n in 1:N) {
      int e = el[n]; int s = sp_of_ind[r_ind[e]];
      mu[n] = logA_e[e] + rich_log(nu_sp[s], -k_sp[s] * (age0[n] + miss[e] - t50_sp[s]) + g50_sp[s]);
    }
    size ~ lognormal(mu, sigma[r_sdim[el]]);
  }
  // count-and-size elements
  for (e in 1:EC) {
    int i = c_ind[e]; int s = sp_of_ind[i];
    real lrel = rich_log(nu_sp[s], -k_sp[s] * (Tage[i] + c_tadd[e] - t50_sp[s]) + g50_sp[s]);
    if (c_has_count[e] == 1) target += normal_lpdf(Tage[i] - c_count[e] | 0, c_resorb[e] == 1 ? tau_off : tau_fresh);
    if (c_has_size[e] == 1) {
      real bel = c_elt[e] > 0 ? beta_el[c_elt[e]] : 0;
      c_logsize[e] ~ normal(mu_logA[c_dim[e]] + theta[s, 1] + bel + a_ind[i] + lrel,
                            sqrt(square(sd_el[i]) / c_nrep[e] + square(sigma_cs) + square(c_logsd[e])));
    }
    if (c_efs[e] == 1)      target += log_inv_logit((exp(lrel) - efs_thr) / efs_slope);
    else if (c_efs[e] == 0) target += log1m_inv_logit((exp(lrel) - efs_thr) / efs_slope);
  }
}
generated quantities {
  vector[K] lambda = square(sd_phy) ./ (square(sd_phy) + square(sd_res));
  corr_matrix[K] Omega = multiply_lower_tri_self_transpose(L_Om);
  vector[EC] c_lost;                       // marks lost from count-and-size elements (inner or outer)
  vector[EC] c_rel;                        // fraction of the asymptote attained at death
  matrix[S,K] dev; matrix[B,K] delta;
  for (e in 1:EC) {
    int i = c_ind[e]; int s = sp_of_ind[i];
    c_lost[e] = c_has_count[e] == 1 ? Tage[i] - c_count[e] : 0;
    c_rel[e] = exp(rich_log(nu_sp[s], -k_sp[s] * (Tage[i] + c_tadd[e] - t50_sp[s]) + g50_sp[s]));
  }
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
