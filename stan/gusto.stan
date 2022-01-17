functions {
}
data {
  int<lower=1> N;  // total number of observations
  int Y[N];  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
  int n_risks; // number of baseline risk points to compute delta events
}
transformed data {
}
parameters {
  vector[K] b;  // population-level effects
}
transformed parameters {
}
model {
  // likelihood including constants
  if (!prior_only) {
    target += bernoulli_logit_glm_lpmf(Y | X, 0, b);
  }
  // priors including constants
  target += normal_lpdf(b[1] | -2.5, 0.75);
  target += normal_lpdf(b[2] | 0, 0.5);
  target += normal_lpdf(b[3] | 0, 0.5);
  target += normal_lpdf(b[4] | 0, 0.5);
  target += normal_lpdf(b[5] | 0, 0.5);
  target += normal_lpdf(b[6] | 0, 0.5);
  target += normal_lpdf(b[7] | 0, 0.5);
  target += normal_lpdf(b[8] | 0, 0.5);
  target += normal_lpdf(b[9] | 0, 0.5);
  target += normal_lpdf(b[10] | 0, 0.5);
  target += normal_lpdf(b[11] | 0, 0.5);
  target += normal_lpdf(b[12] | 0, 0.5);
}
generated quantities {
  vector[n_risks] baseline_risks;
  vector[n_risks] baseline_events;
  vector[n_risks] treated_risks;
  vector[n_risks] treated_events;
  for (i in 1:n_risks) {
    real r = i/100.0;
    real treated_odds = (r/(1-r)) * exp(b[2]);
    
    baseline_risks[i] = r;
    baseline_events[i] = binomial_rng(1000, r);
    
    treated_risks[i] = treated_odds / (1 + treated_odds);
    treated_events[i] = binomial_rng(1000, treated_risks[i]);
  }
}
