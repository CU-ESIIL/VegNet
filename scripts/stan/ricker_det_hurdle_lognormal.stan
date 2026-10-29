// Deterministic-skeleton multispecies Ricker competition model
// with a hurdle (zero / log-normal mixture) observation model.
// This is ricker_ss_hurdle_lognormal.stan with process error fixed to
// zero (sigma_proc = 0): all variation is observation error.
//
// Process:  x[p,t] = x[p,t-1] + r - A * exp(x[p,t-1]) + B * X_env,  x = log N
// Observe:  y = 0                  with prob pi = inv_logit(theta0 - theta1 * x)
//           y ~ LN(x, sigma_obs)   with prob 1 - pi
//
// Notes:
// - Scale densities (e.g., divide by a typical max) so exp(x) and the A
//   priors are O(1); A is then on the scaled-density scale.
// - A is unbounded (negative entries = facilitation). With deterministic
//   dynamics, negative or near-zero A can make x drift or blow up.
// - Missing surveys: set obs = 0 and put any placeholder (e.g., 0) in y.
// - Latent N is the median of the observation; for mean-unbiased
//   observations use lognormal(x - 0.5 * square(sigma_obs), sigma_obs).
data {
  int<lower=1> P;                               // plots
  int<lower=2> T;                               // time steps
  int<lower=1> S;                               // species
  int<lower=1> K;                               // environmental covariate
  array[P, T, S] real<lower=0> y;               // observed densities (scaled)
  array[P, T, S] int<lower=0, upper=1> obs;     // 1 = observed, 0 = missing
  array[P, T, K] real X_env;                    // environmental variable
  // prior hyperparameters for initial log densities
  real mu0_mean;
  real<lower=0> mu0_sd;
}

parameters {
  vector[S] r;                                  // intrinsic growth (log scale)
  matrix<lower=0>[S, S] A;                               // A[i,j]: effect of j on i
  matrix[S, K] B;                               // B[i, j] effect of cov j on species i
  vector<lower=0>[S] sigma_obs;                 // log-normal observation sd
  vector[S] theta0;                             // logit Pr(zero) intercept
  vector<lower=0>[S] theta1;                    // decline in Pr(zero) with log density

  array[P] vector[S] x1;                        // initial log densities
}

transformed parameters {
  array[P, T] vector[S] x;                      // latent log densities
  for (p in 1:P) {
    x[p, 1] = x1[p];
    for (t in 2:T) {
      x[p, t] = x[p, t - 1] + r - A * exp(x[p, t - 1])
        + B * to_vector(X_env[p, t - 1]);
    }
  }
}

model {
  // priors
  r ~ normal(0.5, 1);
  to_vector(A) ~ normal(0, 0.5);
  to_vector(B) ~ normal(0, 0.5);
  sigma_obs ~ normal(0, 0.5);
  theta0 ~ normal(-2, 1.5);
  theta1 ~ normal(0, 1);
  for (p in 1:P) {
    x1[p] ~ normal(mu0_mean, mu0_sd);
  }

  // hurdle observation likelihood
  for (p in 1:P) {
    for (t in 1:T) {
      for (s in 1:S) {
        if (obs[p, t, s] == 1) {
          real eta = theta0[s] - theta1[s] * x[p, t, s];
          if (y[p, t, s] == 0) {
            target += log_inv_logit(eta);
          } else {
            target += log1m_inv_logit(eta)
                      + lognormal_lpdf(y[p, t, s] | x[p, t, s], sigma_obs[s]);
          }
        }
      }
    }
  }
}

generated quantities {
  array[P, T, S] real y_rep;                    // posterior predictive draws
  array[P, T, S] real log_lik;                  // pointwise log-lik (for loo)
  for (p in 1:P) {
    for (t in 1:T) {
      for (s in 1:S) {
        real eta = theta0[s] - theta1[s] * x[p, t, s];
        y_rep[p, t, s] = bernoulli_logit_rng(eta)
                         ? 0 : lognormal_rng(x[p, t, s], sigma_obs[s]);
        if (obs[p, t, s] == 1) {
          log_lik[p, t, s] = y[p, t, s] == 0
            ? log_inv_logit(eta)
            : log1m_inv_logit(eta)
              + lognormal_lpdf(y[p, t, s] | x[p, t, s], sigma_obs[s]);
        } else {
          log_lik[p, t, s] = 0;
        }
      }
    }
  }
}
