
data {
  int<lower=1> N;              // Number of patients
  int<lower=1> K;              // Number of dose levels
  array[N] int<lower=1, upper=K> dose;  // dose level for each patient
  array[N] real<lower=0> time;       // Survival/censoring times
  array[N] int<lower=0, upper=1>  event; // Event indicator (1=event, 0=censored)
  
  // Safety monitoring parameters
  real<lower=0, upper=1> phi;        // Maximum acceptable toxicity rate
  real<lower=0, upper=1> cp;         // Safety threshold probability
  array[K] int<lower=0> xk;          // Number of toxicities at each dose
  array[K] int<lower=0> mk;          // Total patients treated at each dose
  
  // Dose elimination indicator
  array[K] int<lower=0, upper=1> dose_eliminated; // 1 if dose is eliminated due to safety
}

parameters {
  vector[K-1] beta;            // Coefficients for doses 2..K (dose 1 is reference)
}

transformed parameters {
  // Safety probabilities for each dose level
  vector[K] safety_prob;
  
  // Calculate posterior probability that toxicity rate exceeds phi
  for (k in 1:K) {
    if (dose_eliminated[k] == 1) {
      // If dose is already eliminated, set safety probability to 1
      safety_prob[k] = 1.0;
    } else {
      // Calculate Pr(pk > phi) using Beta distribution
      real alpha_post = phi + xk[k];
      real beta_post = 1 - phi + mk[k] - xk[k];
      
      // Use complementary CDF: 1 - Beta_cdf(phi | alpha_post, beta_post)
      safety_prob[k] = 1 - beta_cdf(phi | alpha_post, beta_post);
    }
  }
}

model {
  // Priors for beta coefficients
  beta ~ normal(0, 2^0.5);         // Weakly informative priors
  
  // Partial likelihood for Cox model - implemented directly
  for (i in 1:N) {
    if (event[i] == 1 && dose_eliminated[dose[i]] == 0) {  // Only contribute for observed events
      real numerator;
      real denominator = 0.0;
      
      // Linear predictor for the event case (numerator)
      if (dose[i] == 1) {
        numerator = 0.0;  // reference category
      } else {
        numerator = beta[dose[i] - 1];
      }
      
      // Sum over risk set (denominator)
      for (j in 1:N) {
        if (time[j] >= time[i] && dose_eliminated[dose[j]] == 0) {  // j is in risk set at time[i]
          real risk_contrib;
          if (dose[j] == 1) {
            risk_contrib = 0.0;  // reference category
          } else {
            risk_contrib = beta[dose[j] - 1];
          }
          denominator += exp(risk_contrib);
        }
      }
      
      // Add partial likelihood contribution
      if (denominator > 0) {
        target += numerator - log(denominator);
      }
    }
  }
}

generated quantities {
  vector[K] log_hr;           // Log hazard ratios (relative to dose 1)
  vector[K-1] adj_log_hr;     // Adjacent dose log hazard ratios
  vector[K-1] hr_adj;         // Adjacent dose hazard ratios

  vector[K] all_betas;        // All beta coefficients including reference
  int<lower=0, upper=K> best_dose; // Dose with lowest beta (lowest hazard)
  
  // Safety indicators
  array[K] int<lower=0, upper=1> unsafe_dose; // 1 if dose is unsafe (Pr(pk > phi) > cp)
  int<lower=0> num_unsafe_doses;              // Number of unsafe doses
  
  // Calculate hazard ratios relative to reference (dose 1)
  log_hr[1] = 0.0;            // Reference dose
  for (k in 2:K) {
    log_hr[k] = beta[k-1];
  }
  
  // Calculate adjacent dose hazard ratios
  for (k in 1:(K-1)) {
    adj_log_hr[k] = log_hr[k+1] - log_hr[k];
    hr_adj[k] = exp(adj_log_hr[k]);
  }
  
  // Calculate which doses are unsafe
  for (k in 1:K) {
    unsafe_dose[k] = (safety_prob[k] > cp) ? 1 : 0;
  }
  num_unsafe_doses = sum(unsafe_dose);
  
  // Calculate which dose has the lowest beta (lowest hazard)
  // Create vector of all betas (dose 1 has beta = 0)
  all_betas[1] = 0.0;  // Reference dose
  for (k in 2:K) {
    all_betas[k] = beta[k-1];
  }
  
  // Find the best dose among non-eliminated and safe doses
  best_dose = 0;  // 0 indicates no safe dose found
  for (k in 1:K) {
    if (dose_eliminated[k] == 0 && unsafe_dose[k] == 0) {
      if (best_dose == 0 || all_betas[k] < all_betas[best_dose]) {
        best_dose = k;
      }
    }
  }
}