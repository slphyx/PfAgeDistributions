functions {
  // White 1992: Gaussian reference brood plus shifted adjacent
  // generations, each weighted by the per-cycle multiplication factor.
  // Midpoint quadrature: common factor delta-age cancels in probabilities.
  vector age_bins(real mu, real sigma, real M,
                  real onset, real decay,
                  vector ages, array[] int bin) {
    int G = num_elements(ages);
    vector[7] masses = rep_vector(0.0, 7);
    for (g in 1:G) {
      real log_f = negative_infinity();
      for (r in -4:4) {
        log_f = log_sum_exp(log_f,
                r * log(M) + normal_lpdf(ages[g] | mu - 48.0*r, sigma));
      }
      masses[bin[g]] += exp(log_f - decay * fmax(0.0, ages[g] - onset));
    }
    return masses;
  }

  real circulation_fraction(real mu, real sigma, real M,
                            real onset, real decay,
                            vector ages, array[] int bin) {
    vector[7] visible = age_bins(mu, sigma, M, onset, decay, ages, bin);
    vector[7] total = age_bins(mu, sigma, M, 48.0, 0.0, ages, bin);
    return sum(visible) / sum(total);
  }

  // Orthonormal Helmert contrasts: 7 proportions -> 6 independent
  // coordinates; avoids double-counting the sum-to-one constraint.
  vector helmert(vector proportions) {
    vector[6] z;
    for (j in 1:6) {
      real subtotal = 0.0;
      for (h in 1:j) subtotal += proportions[h];
      z[j] = (subtotal - j * proportions[j+1]) / sqrt(1.0*j*(j+1));
    }
    return z;
  }
}

data {
  int<lower=1> N;
  int<lower=1> G;
  real<lower=1> M;
  vector[G] ages;
  array[G] int<lower=1, upper=7> bin;
  matrix[N, 7] observed_proportions;
  vector<lower=0>[N] observed_total_density;
}

parameters {
  vector<lower=0, upper=48>[N] mu;
  vector<lower=1, upper=18>[N] sigma;
  real<lower=0.001, upper=0.25> stage_error_sd;
  real<lower=0.05, upper=2.0> k_s; // shared by all 12 patients
}

model {
  // Working observation model on six orthonormal composition contrasts.
  // These are percentages, NOT made-up multinomial counts.
  mu ~ normal(26, 10);
  sigma ~ lognormal(log(7), 0.6);
  stage_error_sd ~ normal(0, 0.05);
  k_s ~ lognormal(log(0.6), 0.6);
  for (i in 1:N) {
    real onset = 26.0;
    real decay = k_s;
    vector[7] masses = age_bins(mu[i], sigma[i], M,
                                onset, decay, ages, bin);
    vector[7] p = masses / sum(masses);
    vector[7] y_i = to_vector(observed_proportions[i]');
    target += normal_lpdf(helmert(y_i) | helmert(p), stage_error_sd);
  }
}

generated quantities {
  vector[N] q_circulating;
  vector[N] total_density_N0;
  vector[N] sequestered_density;
  matrix[N,7] predicted_stage_proportions;
  vector[N] log_lik_patient;
  for (i in 1:N) {
    real onset = 26.0;
    real decay = k_s;
    vector[7] masses = age_bins(mu[i], sigma[i], M,
                                onset, decay, ages, bin);
    vector[7] p = masses / sum(masses);
    vector[7] y_i = to_vector(observed_proportions[i]');
    q_circulating[i] = circulation_fraction(mu[i], sigma[i], M,
                                onset, decay, ages, bin);
    total_density_N0[i] = observed_total_density[i] / q_circulating[i];
    sequestered_density[i] = total_density_N0[i] - observed_total_density[i];
    predicted_stage_proportions[i] = p';
    log_lik_patient[i] = normal_lpdf(helmert(y_i) | helmert(p),
                                     stage_error_sd);
  }
}
