// LAG-WIDTH growth-rate model (the native "how bones grow" target).
// Response: log radial LAG width = log(annual radial apposition, mm/yr), the
// least-derived growth signal (no mass imputation). Modeled as declining with
// age (exponential decline in width => von-Bertalanffy-like radius), hierarchical
// by specimen within species, with PHYLOGENETIC structure on the growth-rate
// LEVEL -- the mass-free analogue of "maximum growth rate" in the metabolism
// literature. EFS rings (growth ceased) are excluded upstream.
//
//   log w_ij = alpha + a_phylo[tip] + a_sp[species] + a_id[specimen]
//                    + (beta + b_sp[species]) * age_ij + e_ij
//
// lambda_rate = sigma_phy^2 / (sigma_phy^2 + sd_sp^2): phylogenetic signal in the
// growth-rate level. level@ref = fitted log radial rate at a reference age.
data {
  int<lower=1> N; int<lower=1> S; int<lower=1> I; int<lower=1> T;
  vector[N] logw;
  vector[N] age;
  array[N] int<lower=1,upper=S> sp;
  array[N] int<lower=1,upper=I> id;
  array[S] int<lower=1,upper=T> tip_of_sp;
  matrix[T,T] Lphylo;                      // lower Cholesky of phylo correlation
  real ref_age;                            // reference age for the reported rate level
}
parameters {
  real alpha; real beta;
  real<lower=0> sigma;
  real<lower=0> sigma_phy; vector[T] z_phy;
  real<lower=0> sd_sp;     vector[S] z_sp;
  real<lower=0> sd_id;     vector[I] z_id;
  real<lower=0> sd_beta;   vector[S] z_beta;
}
transformed parameters {
  vector[T] a_phylo = sigma_phy * (Lphylo * z_phy);
  vector[S] a_sp = sd_sp * z_sp;
  vector[S] b_sp = beta + sd_beta * z_beta;
  vector[I] a_id = sd_id * z_id;
}
model {
  alpha ~ normal(0, 3); beta ~ normal(0, 1); sigma ~ student_t(3, 0, 1);
  sigma_phy ~ student_t(3, 0, 1); z_phy ~ std_normal();
  sd_sp ~ student_t(3, 0, 1);     z_sp ~ std_normal();
  sd_id ~ student_t(3, 0, 0.5);   z_id ~ std_normal();
  sd_beta ~ student_t(3, 0, 0.5); z_beta ~ std_normal();
  {
    vector[N] mu;
    for (n in 1:N)
      mu[n] = alpha + a_phylo[tip_of_sp[sp[n]]] + a_sp[sp[n]] + a_id[id[n]] + b_sp[sp[n]] * age[n];
    logw ~ normal(mu, sigma);
  }
}
generated quantities {
  real lambda_rate = square(sigma_phy) / (square(sigma_phy) + square(sd_sp));
  vector[S] lograte_ref;                   // log radial apposition (mm/yr) at ref_age, per species
  for (s in 1:S)
    lograte_ref[s] = alpha + a_phylo[tip_of_sp[s]] + a_sp[s] + b_sp[s] * ref_age;
}
