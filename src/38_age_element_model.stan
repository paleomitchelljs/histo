// Model 38: model 37 + four changes (JM 2026-11-16, in one fit by JM's choice; docs/allometry_prior_plan.md step 1b).
// (f) ELEMENT SHIFT: inner AND outer missing marks for every ring element. R_e = Tage - marks preserved is split into
//     inner = f_e R_e (shifts the element's marks along the curve) and outer = (1 - f_e) R_e (marks lost at the
//     surface: compaction in the OCL/EFS, outer cortex removed, not resolved). Density HN(inner | tau_in) HN(outer |
//     tau_outer) R_e (Jacobian). Model 37 had inner only (outer only for elements flagged eroded), which aligned every
//     element at its outer surface; Orthosuchus SAM-PK-K409 shows the counts differ mostly at the outside. tau_outer =
//     the resorption scale for eroded elements (as before), else tau_out ~ half-normal(0, 1). Count-and-size: an exact
//     count uses the convolution of the two half-normals (Tage - count = inner + outer); interval counts use a
//     half-normal with scale sqrt(tau^2 + tau_out^2) (+ s_open^2 for open counts), an approximation.
// (g) nu POOLED BY CLADE, OFF THE TREE: nu_sp = mu_nu + sd_nu_cl z_nu_cl[clade] + sd_nu_sp z_nu_sp[species]
//     (clades as in M7). The tree block now holds logA, logr, t50 only (K = 3; Omega 3 x 3; no nu-trait correlations).
//     nu was prior-dominated in most species and its tree SD drove the r37/r38 divergences.
// (h) Tighter scales for nu: sd_nu_cl ~ half-normal(0, 0.5), sd_nu_sp ~ half-normal(0, 0.3) (were sd_phy, sd_res ~
//     half-normal(0, 1)).
// (i) kappa_raw exists only when est_kappa = 1 (it was a free parameter with no effect on the likelihood).
//
// ---- model 37 header follows ----
// Model 37: model 36 + two fixes. Written 2026-11-15.
// (d) OPEN COUNTS ("at least lo"). Models 35/36 entered them as Pr(P >= lo) = F(Tage - lo), which tends to 1 as Tage
//     grows and replaced the lost-mark density, so individuals with only open counts had no proper age prior
//     (r36: Hofmeyria BP/1/4404 ~3600 yr, Theriognathus BP/1/719 ~74000 yr). Now Tage - lo = lost + extra preserved
//     marks beyond lo, approximated as half-normal with scale sqrt(tau^2 + s_open^2) (s_open data, default 2 marks).
//     Bounded ranges are unchanged.
// (e) CHANNEL EXPONENTS ONLY WHERE IDENTIFIED. In r36 the tibia circumference exponent (2.03) moved with the curve
//     shape of tibia-only species (Maiasaura t50, r = -0.49; Coelophysis t50) because for a species observed in one
//     channel rel^beta is nearly another Richards curve (docs/model_revision_plan.md M1, identifiability). The
//     free exponent beta_c now applies only in species with sized data in the reference channel of c's
//     co-occurrence component (femur circumference for component 1); in other species every channel has beta = 1,
//     i.e. the latent curve is that species' own channel. Per-element index r_bi / c_bi (0 = beta 1).
//
// ---- model 36 header follows ----
// Model 36: model 35 + (a) species-specific element offsets, (b) non-centred species shape nu,
// (c) periosteal sizes not counted as marks (data side, src/36_helpers.R). Written 2026-11-14.
// (a) ELEMENT OFFSETS BY SPECIES. Model 32/35 used one global offset per element type (beta_el), so species
//     with unusual limb proportions (Allosaurus forelimb; Jeholosaurus tibia vs fibula) pushed systematic
//     differences into the element-level SD. Now, for a species with two or more sized long-bone types,
//     every type other than the species' reference type (femur if present, else its best-sampled type) gets
//        beta_el[type] + gamma[species, type],   gamma ~ normal(0, sd_gam)
//     (exchangeable across species; a tree term is a later step). Relative channels are unaffected.
// (b) NON-CENTRED nu. The species trait block is still one multivariate normal with the tree covariance,
//     factorised as  p(logA, logr, t50) * p(nu | logA, logr, t50):  the first three are sampled as before
//     (centred, logA on the absolute scale); nu = conditional mean + chol(conditional covariance) * z_nu,
//     z_nu ~ std_normal(). Same joint distribution; the funnel between nu and its SDs is removed.
// (c) Data only: the size of the periosteal surface (Maiasaura) stays in a ring series at age marks + 0.5 but is
//     excluded from the number of marks preserved (r_last), so lost marks are no longer off by half a year.
//
// ---- model 35 header follows ----
// Model 35: model 32 + the changes of docs/model_revision_plan.md that current data support.
// Written 2026-11-14. Model 32 (src/32_age_element_model.stan) is the baseline; everything not
// listed here is unchanged from it (curve, species traits on the tree, marginal sampling form,
// individual age identity, count-and-size module, three-state EFS, generated branch changes).
//
// M1  SIZE CHANNELS WITH ALLOMETRIC EXPONENTS. A channel is element type x size dimension
//     (body measures: the dimension alone). Each element's size follows the species' latent
//     relative curve raised to its channel exponent:
//        log S_e(t) = logA_e + beta[ch] * log rel_i(t)
//     beta = 1 is isometry with the latent curve (= model 32 everywhere). beta is FIXED at 1 for
//     the reference channel (femur circumference) and for one channel in every group of channels
//     that never co-occurs with the reference in any species (there beta is not identified);
//     the others are free, log beta ~ normal(0, 0.3), one value per channel across all species
//     (no tree term yet). The asymptote intercepts are as in model 32 (mu_logA[dimension] +
//     beta_el[element type]).
//     RELATIVE channels (dimensions scaled to the largest specimen of a sample: pct_max_size,
//     pct_bone_length) do not inform the species asymptote: their level is a free intercept per
//     species x channel, o_rel ~ normal(log 100, 1). They inform shape, timing and individual
//     variation only. (Plan, data handling: "the largest sampled animal is an adult" is not assumed.)
// M2  INDIVIDUAL TIMING. t50_i = t50_sp + u_i, u_i ~ normal(0, sd_u). Explicit (non-centred) for
//     individuals with a ring series or more than one sized element. For individuals with a single
//     sized observation u_i is integrated out by the delta method: the size variance gains
//     (beta * d lrel / d t50)^2 sd_u^2 / nrep. For resorbed ring series (tpar 0) the age
//     reparametrization now uses t50_i (unit Jacobian shear, as before).
//     NOT included: an individual effect on logr (plan: t50 first).
// M5  INTERVAL COUNTS (count-and-size elements). A printed range [lo, hi] or "at least lo" is an
//     interval-censored observation of the preserved marks P = Tage - lost, lost ~ half-normal(tau):
//        Pr(lo <= P <= hi) = F(Tage - lo) - F(Tage - hi),   F = half-normal CDF (F(x<0) = 0);
//     exact counts keep the half-normal density. The age floor uses lo. NOT included: a counting-
//     error term (needs replicate-count data, plan M5).
// M6  SEX. logA_e += sexc_i * g_sex[species] for species with both sexes, sexc = -0.5 female,
//     +0.5 male, 0 unknown or pooled. g_sex ~ normal(gbar_sex, sd_sex). Sex-specific timing: not yet.
// M7  LOST-MARK SCALE BY ELEMENT TYPE AND CLADE. For resorbed elements
//        log tau_e = log_tau0 + sd_tel z_tel[element type] + sd_tcl z_tcl[clade]
//     (model 32: one tau_off). The outer-loss term of eroded ring elements uses the same tau_e.
//     Non-resorbed elements keep tau_fresh (data).
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
  // d rich_log / dx (used only in a variance term, so the floor is approximated)
  real rich_dlog_dx(real nu, real x) {
    real ex = exp(fmin(x, 50));
    return -ex / fmax(1 + nu * ex, 1e-3);
  }
  // log of the half-normal CDF difference F(a) - F(b), b < a, F(x) = 0 for x <= 0
  real hn_interval_lp(real a, real b, real tau) {
    real bb = fmax(b, 0.0);
    if (a <= bb + 1e-9) return log(1e-12);
    return log2() + log_diff_exp(normal_lcdf(a | 0, tau), normal_lcdf(bb | 0, tau));
  }
  // log density of d = X + Y, X ~ half-normal(a), Y ~ half-normal(b), d > 0:
  // 4 phi_s(d) [Phi(d b / (a s)) - Phi(-d a / (b s))], s^2 = a^2 + b^2
  real hn2_lpdf(real d, real a, real b) {
    real s2 = sqrt(square(a) + square(b));
    if (d <= 1e-9) return log(1e-12);
    return log(4) + normal_lpdf(d | 0, s2) + log_diff_exp(normal_lcdf(d * b / (a * s2) | 0, 1), normal_lcdf(-d * a / (b * s2) | 0, 1));
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
  int<lower=0> IM; array[I] int<lower=0,upper=IM> im_of_ind;
  int<lower=0> IU; array[I] int<lower=0,upper=IU> u_of_ind;        // M2: explicit timing effect (0 = integrated)
  vector[I] sexc; int<lower=0> NSX; array[S] int<lower=0,upper=NSX> sx_of_sp;   // M6
  // channels (M1)
  int<lower=1> C; int<lower=0> NB; array[C] int<lower=0,upper=NB> beta_idx;     // 0 = beta fixed at 1
  int<lower=0> NO;                                                              // relative-channel intercepts
  // lost-mark scale groups (M7)
  int<lower=1> NTE; int<lower=1> NTC;
  int<lower=0> NG;                                                 // (a) species x element-type offsets
  // ring elements
  int<lower=0> ER;
  array[ER] int<lower=1,upper=I> r_ind; array[ER] int<lower=1,upper=D> r_dim; array[ER] int<lower=0,upper=L> r_elt;
  array[ER] int<lower=1,upper=C> r_ch; array[ER] int<lower=0,upper=NO> r_relch;
  array[ER] int<lower=1,upper=NTE> r_te; array[ER] int<lower=1,upper=NTC> r_tc; array[ER] int<lower=0,upper=NG> r_g;
  int<lower=1> DR; array[ER] int<lower=1,upper=DR> r_sdim;
  array[ER] int<lower=0,upper=2> r_efs; array[ER] int<lower=0,upper=1> r_resorb;
  int<lower=0> NERO; array[ER] int<lower=0,upper=NERO> r_ero;
  vector[ER] r_last; vector[ER] r_logmax;
  int<lower=0> N; vector<lower=0>[N] size; vector[N] age0; array[N] int<lower=1,upper=ER> el;
  // count-and-size elements
  int<lower=0> EC;
  array[EC] int<lower=1,upper=I> c_ind; array[EC] int<lower=1,upper=D> c_dim; array[EC] int<lower=0,upper=L> c_elt;
  array[EC] int<lower=1,upper=C> c_ch; array[EC] int<lower=0,upper=NO> c_relch;
  array[EC] int<lower=1,upper=NTE> c_te; array[EC] int<lower=1,upper=NTC> c_tc; array[EC] int<lower=0,upper=NG> c_g;
  array[EC] int<lower=0,upper=2> c_efs; array[EC] int<lower=0,upper=1> c_resorb;
  array[EC] int<lower=0,upper=2> c_cnt_type;                       // 0 none, 1 exact, 2 interval (hi may be "at least")
  vector[EC] c_count; vector[EC] c_lo; vector[EC] c_hi; array[EC] int<lower=0,upper=1> c_open;   // c_open 1: hi = infinity
  array[EC] int<lower=0,upper=1> c_has_size; vector[EC] c_logsize; vector<lower=0>[EC] c_logsd;
  vector<lower=1>[EC] c_nrep; vector<lower=0>[EC] c_tadd;
  array[S] int<lower=1,upper=D> main_dim;
  real<lower=0,upper=1> efs_thr; real<lower=0> efs_slope;
  real<lower=0> tau_fresh;
  real<lower=0> s_open;                                            // (d)
  array[ER] int<lower=0,upper=NB> r_bi; array[EC] int<lower=0,upper=NB> c_bi;   // (e) 0 = beta 1
  array[S] int<lower=1,upper=NTC> cl_of_sp;                       // (g) clade of each species
}
transformed data {
  int K = 3;                                                       // (g) tree traits: logA, logr, t50
  vector[B] logbr = log(brlen);
  vector[K] sd_scale = [1.0, 0.5, 3.0]';
  matrix[S,B] Ps = P[tip_of_sp, ];
}
parameters {
  vector[D] mu_logA; real mu_logr; real mu_nu; real mu_t50;
  vector<lower=0>[K] sd_phy; vector<lower=0>[K] sd_res;
  cholesky_factor_corr[K] L_Om;
  matrix[S,3] theta_o;                                             // logA (absolute, main dimension), logr, t50
  real<lower=0> sd_nu_cl; real<lower=0> sd_nu_sp; vector[NTC] z_nu_cl; vector[S] z_nu_sp;   // (g)
  real<lower=0> sd_gam; vector[NG] z_gam;                          // (a)
  vector[L] beta_el;
  real<lower=0> sd_logA_ind; real<lower=0,upper=1> rho_ind;
  vector[IM] z_ind;
  vector[ER] logA_e;
  array[est_kappa] real<lower=0,upper=1> kappa_raw;               // (i)
  real log_tau0; real<lower=0> sd_tel; real<lower=0> sd_tcl; vector[NTE] z_tel; vector[NTC] z_tcl;   // M7
  vector[I] t_raw;
  vector<lower=0,upper=1>[ER] f_el;                                // (f) inner share of R_e
  real<lower=0> tau_out;                                           // (f) outer-loss scale, non-eroded elements
  vector<lower=0>[DR] sigma;
  real<lower=0> sigma_cs;
  vector[NB] log_beta;                                             // M1
  vector[NO] o_rel;                                                // M1 relative channels
  real<lower=0> sd_u; vector[IU] z_u;                              // M2
  real gbar_sex; real<lower=0> sd_sex; vector[NSX] z_sex;          // M6
}
transformed parameters {
  real kappa = est_kappa ? kappa_raw[1] : 1.0;
  vector[B] bscale;
  matrix[S,4] theta;                                               // logA, logr, nu, t50 (nu not on the tree)
  vector[S] logr_sp = theta_o[, 2]; vector[S] t50_sp = theta_o[, 3];
  vector[S] nu_sp = mu_nu + sd_nu_cl * z_nu_cl[cl_of_sp] + sd_nu_sp * z_nu_sp;   // (g)
  real lp_o;                                                       // log p(logA, logr, t50 | tree)
  vector[NG] gam = sd_gam * z_gam;
  vector[S] k_sp; vector[S] g50_sp;
  vector[C] beta_ch;
  vector[NTE] tau_el = exp(sd_tel * z_tel); vector[NTC] tau_cl = exp(sd_tcl * z_tcl);
  vector[NSX] g_sex = gbar_sex + sd_sex * z_sex;
  vector[I] u_ind; vector[I] t50_i;
  vector[I] m_off;
  vector<lower=0>[I] Tage;
  vector[I] a_ind; vector<lower=0>[I] sd_el;
  vector[ER] miss;
  {
    vector[B] bvar0 = exp(kappa * logbr);
    vector[B] bsc = bvar0 / mean(P * bvar0);
    matrix[K,K] Vphy = multiply_lower_tri_self_transpose(diag_pre_multiply(sd_phy, L_Om));
    matrix[S,S] Cm = Ps * diag_post_multiply(Ps, bsc)';
    matrix[S * K, S * K] Sig = build_sigma(Vphy, Cm, sd_res);
    vector[3 * S] m_o = append_row(rep_vector(0, S), append_row(rep_vector(mu_logr, S), rep_vector(mu_t50, S)));
    vector[3 * S] th_o = append_row(theta_o[, 1] - mu_logA[main_dim], append_row(theta_o[, 2], theta_o[, 3]));
    lp_o = multi_normal_cholesky_lpdf(th_o | m_o, cholesky_decompose(Sig));
  }
  theta[, 1] = theta_o[, 1] - mu_logA[main_dim]; theta[, 2] = theta_o[, 2]; theta[, 3] = nu_sp; theta[, 4] = theta_o[, 3];
  for (c in 1:C) beta_ch[c] = beta_idx[c] > 0 ? exp(log_beta[beta_idx[c]]) : 1.0;
  {
    vector[B] bvar = exp(kappa * logbr);
    bscale = bvar / mean(P * bvar);
  }
  for (s in 1:S) {
    g50_sp[s] = gfun(0.5, nu_sp[s]);
    k_sp[s] = 2 * exp(logr_sp[s]) * (gfun(0.25, nu_sp[s]) - gfun(0.75, nu_sp[s]));
  }
  for (i in 1:I) {
    u_ind[i] = u_of_ind[i] > 0 ? sd_u * z_u[u_of_ind[i]] : 0;
    t50_i[i] = t50_sp[sp_of_ind[i]] + u_ind[i];
    m_off[i] = tpar[i] == 0 ? t50_i[i] - t_raw[i] : t_loc[i] + t_scale[i] * t_raw[i];
    Tage[i] = Lmax[i] + log1p_exp(csoft[i] * m_off[i]) / csoft[i];
    if (im_of_ind[i] > 0) { a_ind[i] = sqrt(rho_ind) * sd_logA_ind * z_ind[im_of_ind[i]]; sd_el[i] = sqrt(1 - rho_ind) * sd_logA_ind; }
    else { a_ind[i] = 0; sd_el[i] = sd_logA_ind; }
    if (sx_of_sp[sp_of_ind[i]] > 0) a_ind[i] += sexc[i] * g_sex[sx_of_sp[sp_of_ind[i]]];
  }
  for (e in 1:ER) miss[e] = f_el[e] * (Tage[r_ind[e]] - r_last[e]);   // (f) inner share
}
model {
  mu_logA ~ normal(4, 2); mu_logr ~ normal(-2.5, 1); mu_nu ~ normal(0, 1); mu_t50 ~ normal(5, 5);
  sd_phy ~ normal(0, sd_scale); sd_res ~ normal(0, sd_scale);
  L_Om ~ lkj_corr_cholesky(2);
  target += lp_o;                                                  // (b) p(logA, logr, t50); unit-Jacobian shear for logA
  sd_nu_cl ~ normal(0, 0.5); sd_nu_sp ~ normal(0, 0.3); z_nu_cl ~ std_normal(); z_nu_sp ~ std_normal();   // (g, h)
  tau_out ~ normal(0, 1);                                          // (f)
  sd_gam ~ normal(0, 0.5); z_gam ~ std_normal();                   // (a)
  beta_el ~ normal(0, 1);
  sd_logA_ind ~ normal(0, 0.5); rho_ind ~ beta(2, 2); z_ind ~ std_normal();
  kappa_raw ~ beta(1, 1);
  log_tau0 ~ normal(log(2), 0.75); sd_tel ~ normal(0, 0.5); sd_tcl ~ normal(0, 0.5); z_tel ~ std_normal(); z_tcl ~ std_normal();
  sigma ~ student_t(3, 0, 0.3); sigma_cs ~ normal(0, 0.15);
  log_beta ~ normal(0, 0.3);
  o_rel ~ normal(log(100), 1);
  sd_u ~ normal(0, 1); z_u ~ std_normal();
  gbar_sex ~ normal(0, 0.2); sd_sex ~ normal(0, 0.2); z_sex ~ std_normal();

  for (i in 1:I) {
    target += log_inv_logit(csoft[i] * m_off[i]);
    if (age_type[i] == 1) Tage[i] ~ normal(age_est[i], age_sd[i]);
  }
  // ring elements
  for (e in 1:ER) {
    int i = r_ind[e]; int s = sp_of_ind[i];
    real bel = (r_elt[e] > 0 ? beta_el[r_elt[e]] : 0) + (r_g[e] > 0 ? gam[r_g[e]] : 0);
    real base = r_relch[e] > 0 ? o_rel[r_relch[e]] : mu_logA[r_dim[e]] + theta[s, 1] + bel;
    real gap = r_logmax[e] - logA_e[e];
    real tau = r_resorb[e] == 1 ? exp(log_tau0) * tau_el[r_te[e]] * tau_cl[r_tc[e]] : tau_fresh;
    logA_e[e] ~ normal(base + a_ind[i], sd_el[i]);
    {                                                              // (f) inner + outer, Jacobian R
      real R = Tage[i] - r_last[e];
      real tau_o = r_ero[e] > 0 ? exp(log_tau0) * tau_el[r_te[e]] * tau_cl[r_tc[e]] : tau_out;
      target += normal_lpdf(miss[e] | 0, tau) + normal_lpdf(R - miss[e] | 0, tau_o) + log(R + 1e-8);
    }
    target += normal_lpdf(fmax(0.0, gap) | 0, 0.07);
    if (r_efs[e] == 1)      target += log_inv_logit((exp(fmin(0.0, gap)) - efs_thr) / efs_slope);
    else if (r_efs[e] == 0) target += log1m_inv_logit((exp(fmin(0.0, gap)) - efs_thr) / efs_slope);
  }
  {
    vector[N] mu;
    for (n in 1:N) {
      int e = el[n]; int i = r_ind[e]; int s = sp_of_ind[i];
      mu[n] = logA_e[e] + (r_bi[e] > 0 ? exp(log_beta[r_bi[e]]) : 1.0) * rich_log(nu_sp[s], -k_sp[s] * (age0[n] + miss[e] - t50_i[i]) + g50_sp[s]);
    }
    size ~ lognormal(mu, sigma[r_sdim[el]]);
  }
  // count-and-size elements
  for (e in 1:EC) {
    int i = c_ind[e]; int s = sp_of_ind[i];
    real x = -k_sp[s] * (Tage[i] + c_tadd[e] - t50_i[i]) + g50_sp[s];
    real bc = c_bi[e] > 0 ? exp(log_beta[c_bi[e]]) : 1.0;
    real lrel = bc * rich_log(nu_sp[s], x);
    real tau = c_resorb[e] == 1 ? exp(log_tau0) * tau_el[c_te[e]] * tau_cl[c_tc[e]] : tau_fresh;
    if (c_cnt_type[e] == 1) target += hn2_lpdf(Tage[i] - c_count[e] | tau, tau_out);              // (f) inner + outer
    else if (c_cnt_type[e] == 2) {
      real tio = sqrt(square(tau) + square(tau_out));
      if (c_open[e] == 1) target += normal_lpdf(Tage[i] - c_lo[e] | 0, sqrt(square(tio) + square(s_open)));   // (d) lost + extra marks
      else target += hn_interval_lp(Tage[i] - c_lo[e], Tage[i] - c_hi[e], tio);
    }
    if (c_has_size[e] == 1) {
      real bel = (c_elt[e] > 0 ? beta_el[c_elt[e]] : 0) + (c_g[e] > 0 ? gam[c_g[e]] : 0);
      real base = c_relch[e] > 0 ? o_rel[c_relch[e]] : mu_logA[c_dim[e]] + theta[s, 1] + bel;
      real vu = u_of_ind[i] == 0 ? square(bc * k_sp[s] * rich_dlog_dx(nu_sp[s], x) * sd_u) / c_nrep[e] : 0;
      c_logsize[e] ~ normal(base + a_ind[i] + lrel,
                            sqrt(square(sd_el[i]) / c_nrep[e] + square(sigma_cs) + square(c_logsd[e]) + vu));
    }
    if (c_efs[e] == 1)      target += log_inv_logit((exp(lrel) - efs_thr) / efs_slope);
    else if (c_efs[e] == 0) target += log1m_inv_logit((exp(lrel) - efs_thr) / efs_slope);
  }
}
generated quantities {
  vector[K] lambda = square(sd_phy) ./ (square(sd_phy) + square(sd_res));
  corr_matrix[K] Omega = multiply_lower_tri_self_transpose(L_Om);
  real tau_off = exp(log_tau0);
  vector[EC] c_lost; vector[EC] c_rel;
  matrix[S,K] dev; matrix[B,K] delta;
  for (e in 1:EC) {
    int i = c_ind[e]; int s = sp_of_ind[i];
    c_lost[e] = c_cnt_type[e] > 0 ? Tage[i] - (c_cnt_type[e] == 1 ? c_count[e] : c_lo[e]) : 0;
    c_rel[e] = exp((c_bi[e] > 0 ? exp(log_beta[c_bi[e]]) : 1.0) * rich_log(nu_sp[s], -k_sp[s] * (Tage[i] + c_tadd[e] - t50_i[i]) + g50_sp[s]));
  }
  {
    matrix[K,K] Lphy = diag_pre_multiply(sd_phy, L_Om);
    matrix[K,K] Vphy = multiply_lower_tri_self_transpose(Lphy);
    matrix[S,S] Cm = Ps * diag_post_multiply(Ps, bscale)';
    matrix[B,K] dstar; matrix[S,K] tstar; matrix[S,K] W;
    dev[, 1] = theta[, 1]; dev[, 2] = theta[, 2] - mu_logr; dev[, 3] = theta[, 4] - mu_t50;   // (g) logA, logr, t50
    for (b in 1:B) for (k in 1:K) dstar[b, k] = normal_rng(0, 1);
    dstar = diag_pre_multiply(sqrt(bscale), dstar) * Lphy';
    tstar = Ps * dstar;
    for (s in 1:S) for (k in 1:K) tstar[s, k] += normal_rng(0, sd_res[k]);
    W = to_matrix(build_sigma(Vphy, Cm, sd_res) \ to_vector(dev - tstar), S, K);
    delta = dstar + diag_pre_multiply(bscale, Ps' * W) * Vphy;
  }
}
