data {
  int<lower=1> N;
  int<lower=1> K;
  array[N] int<lower=1, upper=K> dose;
  array[N] real<lower=0> time;
  array[N] int<lower=0, upper=1> event;

  int<lower=0> N0;
  array[N0] int<lower=1, upper=K> dose0;
  array[N0] real<lower=0> time0;
  array[N0] int<lower=0, upper=1> event0;

  real<lower=0, upper=1> phi;
  real<lower=0, upper=1> cp;
  array[K] int<lower=0> xk;
  array[K] int<lower=0> mk;
  array[K] int<lower=0> x0k;
  array[K] int<lower=0> m0k;

  array[K] int<lower=0, upper=1> dose_eliminated;

  real<lower=0, upper=1> a0;
  int<lower=0, upper=1> use_historical;
}

parameters {
  vector[K-1] beta;
}

transformed parameters {
  vector[K] safety_prob;

  for (k in 1:K) {
    if (dose_eliminated[k] == 1) {
      safety_prob[k] = 1.0;
    } else {
      real historical_weight = use_historical == 1 ? a0 : 0.0;
      real alpha_post = phi + xk[k] + historical_weight * x0k[k];
      real beta_post = 1 - phi + mk[k] - xk[k] + historical_weight * (m0k[k] - x0k[k]);
      safety_prob[k] = 1 - beta_cdf(phi | alpha_post, beta_post);
    }
  }
}

model {
  beta ~ normal(0, sqrt(2));

  for (i in 1:N) {
    if (event[i] == 1) {
      real numerator;
      real denominator = 0.0;

      if (dose[i] == 1) numerator = 0.0;
      else numerator = beta[dose[i] - 1];

      for (j in 1:N) {
        if (time[j] >= time[i]) {
          real risk_contrib;
          if (dose[j] == 1) risk_contrib = 0.0;
          else risk_contrib = beta[dose[j] - 1];
          denominator += exp(risk_contrib);
        }
      }

      if (denominator > 0) target += numerator - log(denominator);
    }
  }

  if (use_historical == 1 && N0 > 0) {
    for (i in 1:N0) {
      if (event0[i] == 1) {
        real numerator;
        real denominator = 0.0;

        if (dose0[i] == 1) numerator = 0.0;
        else numerator = beta[dose0[i] - 1];

        for (j in 1:N0) {
          if (time0[j] >= time0[i]) {
            real risk_contrib;
            if (dose0[j] == 1) risk_contrib = 0.0;
            else risk_contrib = beta[dose0[j] - 1];
            denominator += exp(risk_contrib);
          }
        }

        if (denominator > 0) target += a0 * (numerator - log(denominator));
      }
    }
  }
}

generated quantities {
  vector[K] log_hr;
  vector[K-1] adj_log_hr;
  vector[K-1] hr_adj;
  vector[K] all_betas;
  int<lower=0, upper=K> best_dose;
  array[K] int<lower=0, upper=1> unsafe_dose;
  int<lower=0> num_unsafe_doses;

  log_hr[1] = 0.0;
  for (k in 2:K) log_hr[k] = beta[k - 1];

  for (k in 1:(K - 1)) {
    adj_log_hr[k] = log_hr[k + 1] - log_hr[k];
    hr_adj[k] = exp(adj_log_hr[k]);
  }

  for (k in 1:K) unsafe_dose[k] = (safety_prob[k] > cp) ? 1 : 0;
  num_unsafe_doses = sum(unsafe_dose);

  all_betas[1] = 0.0;
  for (k in 2:K) all_betas[k] = beta[k - 1];

  best_dose = 0;
  for (k in 1:K) {
    if (dose_eliminated[k] == 0 && unsafe_dose[k] == 0) {
      if (best_dose == 0 || all_betas[k] < all_betas[best_dose]) best_dose = k;
    }
  }
}
