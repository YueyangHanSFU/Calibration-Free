# Real-trial-context simulation example for the CFH design.
# Context: CA209-010 / CheckMate 010, nivolumab in metastatic renal cell carcinoma.
# The patient-level data generated below are simulated; only the clinical context
# and calibration targets are taken from the published/registered trial.

setwd("C:/Users/y277h/Desktop/CFH Project/My Code and Results")

library(cmdstanr)
library(survival)
library(dplyr)
library(ggplot2)

## ---------------------------
## Real-trial context settings
## ---------------------------

trial_name <- "CA209-010 / CheckMate 010"
dose_label <- c("0.3 mg/kg", "2 mg/kg", "10 mg/kg")
K <- length(dose_label)

# Registered/published median PFS values, in months, for nivolumab 0.3, 2, and
# 10 mg/kg. These values create a plateau-like dose-efficacy pattern.
median_pfs <- c(2.7, 4.0, 4.2)

# Cox-model parameters relative to the lowest dose. Under an exponential working
# calibration, hazard_k / hazard_1 = median_pfs[1] / median_pfs[k].
beta.true <- log(median_pfs[1] / median_pfs)
base_hazard <- log(2) / median_pfs[1]

# Illustrative toxicity probabilities. Here toxicity is defined as a
# drug-related adverse event leading to discontinuation, using aggregate counts
# reported on ClinicalTrials.gov: 4/59, 5/54, and 4/54.
true.toxicity.rate <- c(4 / 59, 5 / 54, 4 / 54)

# Ordered dose scores used only for the minimum noninferior dose penalty.
dose_score <- c(1.0, 1.1, 1.2)

## ---------------------------
## Design settings
## ---------------------------

nsim <- 1
ncohort <- 56       # Set to 56 to mimic the original trial size of 168 patients.
cohortsize <- 3
followup_time <- 24 # Months; the registered PFS analysis window was about 2 years.

eta <- 1.2
lambda <- 0.7
phi <- 0.33
c_T <- 0.6
phi_a <- phi
phi_b <- 1 - phi

gamma_L <- 1
gamma_R <- 1
prob_cutoff <- 0.5

stan_file <- "Stan Code Partial Likelihood.stan"
mod <- cmdstan_model(stan_file)

extract_beta_matrix <- function(params, K) {
  beta_names <- paste0("beta[", seq_len(K - 1), "]")
  if (!all(beta_names %in% names(params))) {
    stop("Cannot find beta[1], ..., beta[K-1] in the Stan output.")
  }
  beta_mat <- matrix(0, nrow = nrow(params), ncol = K)
  beta_mat[, 2:K] <- as.matrix(params[, beta_names, drop = FALSE])
  beta_mat
}

posterior_best_dose <- function(beta_mat, safe_doses) {
  safe_beta <- beta_mat[, safe_doses, drop = FALSE]
  safe_doses[max.col(-safe_beta, ties.method = "random")]
}

modal_dose <- function(draws, K) {
  tab <- tabulate(as.integer(draws), nbins = K)
  which.max(tab)
}

choose_by_hazard_rule <- function(beta_mat, cidx, max_feasible_dose) {
  if (max_feasible_dose == 1) {
    return(1)
  }

  if (cidx > 1 && cidx < max_feasible_dose) {
    hr_left <- exp(beta_mat[, cidx] - beta_mat[, cidx - 1])
    hr_right <- exp(beta_mat[, cidx + 1] - beta_mat[, cidx])

    prob_left_gt <- mean(hr_left > gamma_L)
    prob_left_lt <- mean(hr_left < gamma_L)
    prob_right_gt <- mean(hr_right > gamma_R)
    prob_right_lt <- mean(hr_right < gamma_R)

    if (prob_left_gt > prob_cutoff && prob_right_gt > prob_cutoff) {
      return(cidx - 1)
    }
    if (prob_left_lt > prob_cutoff && prob_right_lt > prob_cutoff) {
      return(cidx + 1)
    }
    if (prob_left_gt > prob_cutoff && prob_right_lt > prob_cutoff) {
      deescalate_prob <- mean(hr_left) / (mean(hr_left) + mean(1 / hr_right))
      return(ifelse(runif(1) < deescalate_prob, cidx - 1, cidx + 1))
    }

    weight_left <- mean(1 / exp(beta_mat[, cidx - 1]))
    weight_stay <- mean(1 / exp(beta_mat[, cidx]))
    weight_right <- mean(1 / exp(beta_mat[, cidx + 1]))
    probs <- c(weight_left, weight_stay, weight_right) /
      (weight_left + weight_stay + weight_right)
    return(sample(c(cidx - 1, cidx, cidx + 1), size = 1, prob = probs))
  }

  if (cidx == 1) {
    hr_right <- exp(beta_mat[, 2] - beta_mat[, 1])
    if (mean(hr_right < gamma_R) > prob_cutoff) {
      return(2)
    }
    weight_stay <- mean(1 / exp(beta_mat[, 1]))
    weight_right <- mean(1 / exp(beta_mat[, 2]))
    return(sample(c(1, 2), size = 1, prob = c(weight_stay, weight_right)))
  }

  hr_left <- exp(beta_mat[, cidx] - beta_mat[, cidx - 1])
  if (mean(hr_left > gamma_L) > prob_cutoff) {
    return(cidx - 1)
  }
  weight_left <- mean(1 / exp(beta_mat[, cidx - 1]))
  weight_stay <- mean(1 / exp(beta_mat[, cidx]))
  sample(c(cidx - 1, cidx), size = 1, prob = c(weight_left, weight_stay))
}

number_of_patients <- rep(0, K)
dose_selection <- rep(0, K)
dose_selection_MND <- rep(0, K)
early_stop <- 0

# Rows are simulated trials and columns are cohorts. Entry (s, c) records the
# dose level assigned to cohort c in simulated trial s. NA means the trial had
# already stopped before that cohort was enrolled.
cohort_dose_history <- matrix(NA_integer_, nrow = nsim, ncol = ncohort)
colnames(cohort_dose_history) <- paste0("cohort_", seq_len(ncohort))

set.seed(20260811)

for (simu in seq_len(nsim)) {
  cidx <- 1
  survival_outcomes <- vector("list", K)
  censor_status <- vector("list", K)
  toxicity_indicator <- vector("list", K)
  dose_eliminated <- rep(0, K)
  params <- NULL
  beta_mat <- NULL

  for (cohort in seq_len(ncohort)) {
    cohort_dose_history[simu, cohort] <- cidx

    beta_c <- beta.true[cidx]
    true_times <- rexp(cohortsize, rate = base_hazard * exp(beta_c))
    censoring_times <- rep(followup_time, cohortsize)
    observed_time <- pmin(true_times, censoring_times)
    event <- as.numeric(true_times <= censoring_times)
    toxicity <- rbinom(cohortsize, 1, true.toxicity.rate[cidx])

    survival_outcomes[[cidx]] <- c(survival_outcomes[[cidx]], observed_time)
    censor_status[[cidx]] <- c(censor_status[[cidx]], event)
    toxicity_indicator[[cidx]] <- c(toxicity_indicator[[cidx]], toxicity)

    dose_level <- unlist(lapply(seq_len(K), function(k) {
      rep(k, length(survival_outcomes[[k]]))
    }))
    t_array <- unlist(survival_outcomes)
    c_array <- unlist(censor_status)
    xk <- sapply(toxicity_indicator, sum)
    mk <- sapply(toxicity_indicator, length)

    tox_post_prob <- 1 - pbeta(phi, phi_a + xk, phi_b + mk - xk)
    first_overdose <- which(tox_post_prob > c_T)
    if (length(first_overdose) > 0) {
      dose_eliminated[min(first_overdose):K] <- 1
    }

    max_feasible_dose <- ifelse(any(dose_eliminated == 1),
                                min(which(dose_eliminated == 1)) - 1,
                                K)
    if (max_feasible_dose == 0) {
      early_stop <- early_stop + 1
      break
    }
    cidx <- min(cidx, max_feasible_dose)

    stan_data <- list(
      N = length(t_array),
      K = K,
      dose = dose_level,
      time = t_array,
      event = c_array,
      phi = phi,
      cp = c_T,
      xk = xk,
      mk = mk,
      dose_eliminated = dose_eliminated
    )

    fit_mcmc <- mod$sample(
      data = stan_data,
      seed = 10000 + simu * 100 + cohort,
      chains = 4,
      parallel_chains = 4,
      refresh = 0,
      iter_sampling = 5000,
      iter_warmup = 2000,
      thin = 1
    )

    params <- fit_mcmc$draws(format = "df")
    beta_mat <- extract_beta_matrix(params, K)

    safe_doses <- seq_len(max_feasible_dose)
    best_dose_draws <- posterior_best_dose(beta_mat, safe_doses)
    nearby_doses <- intersect(c(cidx - 1, cidx, cidx + 1), safe_doses)
    psi <- mean(best_dose_draws %in% nearby_doses)

    # This matches the manuscript description: hazard comparison with probability
    # psi, and Thompson sampling with probability 1 - psi.
    if (runif(1) < psi) {
      next_dose <- choose_by_hazard_rule(beta_mat, cidx, max_feasible_dose)
    } else {
      next_dose <- sample(best_dose_draws, size = 1)
    }

    cidx <- min(max(next_dose, 1), max_feasible_dose)
  }

  if (is.null(beta_mat)) {
    next
  }

  final_safe <- which(dose_eliminated == 0)
  if (length(final_safe) == 0) {
    next
  }

  final_best_draws <- posterior_best_dose(beta_mat, final_safe)
  best_dose <- modal_dose(final_best_draws, K)

  mean_hazard_ratio <- colMeans(exp(beta_mat))
  mean_hazard_ratio[dose_eliminated == 1] <- Inf
  noninferior <- mean_hazard_ratio <= eta * min(mean_hazard_ratio[final_safe])
  mnd_metric <- mean_hazard_ratio * dose_score^lambda
  MND_best_dose <- which.min(ifelse(noninferior, mnd_metric, Inf))

  number_of_patients <- number_of_patients + sapply(survival_outcomes, length)
  dose_selection[best_dose] <- dose_selection[best_dose] + 1
  dose_selection_MND[MND_best_dose] <- dose_selection_MND[MND_best_dose] + 1
}

real_trial_summary <- data.frame(
  dose = dose_label,
  median_pfs_months = median_pfs,
  beta_true = round(beta.true, 3),
  toxicity_rate = round(true.toxicity.rate, 3),
  selection_percent = round(100 * dose_selection / nsim, 1),
  mnd_selection_percent = round(100 * dose_selection_MND / nsim, 1),
  allocation_percent = round(100 * number_of_patients /
                               sum(number_of_patients), 1)
)

print(real_trial_summary)
cat("Early stopping percentage:", round(100 * early_stop / nsim, 1), "\n")

cohort_dose_record <- data.frame(
  simulation = rep(seq_len(nsim), each = ncohort),
  cohort = rep(seq_len(ncohort), times = nsim),
  dose_level = as.vector(t(cohort_dose_history))
)
cohort_dose_record <- cohort_dose_record[!is.na(cohort_dose_record$dose_level), ]
cohort_dose_record$dose <- factor(
  cohort_dose_record$dose_level,
  levels = seq_len(K),
  labels = dose_label
)
cohort_dose_record$dose_index <- cohort_dose_record$dose_level - 1

cohort_total <- cohort_dose_record %>%
  group_by(cohort) %>%
  summarise(total_trials_not_stopped = n(), .groups = "drop")

cohort_dose_trend <- cohort_dose_record %>%
  group_by(cohort, dose, dose_level) %>%
  summarise(n = n(), .groups = "drop") %>%
  left_join(cohort_total, by = "cohort") %>%
  mutate(percent = 100 * n / total_trials_not_stopped)

trend_grid <- expand.grid(
  cohort = seq_len(ncohort),
  dose = factor(dose_label, levels = dose_label)
)
trend_grid$dose_level <- as.integer(trend_grid$dose)

cohort_dose_trend <- trend_grid %>%
  left_join(cohort_dose_trend, by = c("cohort", "dose", "dose_level")) %>%
  left_join(cohort_total, by = "cohort", suffix = c("", "_grid")) %>%
  mutate(
    total_trials_not_stopped = ifelse(
      is.na(total_trials_not_stopped),
      total_trials_not_stopped_grid,
      total_trials_not_stopped
    ),
    n = ifelse(is.na(n), 0, n),
    percent = ifelse(is.na(percent), 0, percent)
  ) %>%
  select(cohort, dose_level, dose, n, total_trials_not_stopped, percent)

write.csv(
  cohort_dose_record,
  file = "real_trial_cohort_dose_history.csv",
  row.names = FALSE
)
write.csv(
  cohort_dose_trend,
  file = "real_trial_cohort_dose_trend.csv",
  row.names = FALSE
)

# Plot one simulated trial path. The y-axis is coded as 0 for the lowest dose,
# 1 for the middle dose, and 2 for the highest dose.
plot_simulation <- 1
cohort_dose_path <- cohort_dose_record %>%
  filter(simulation == plot_simulation) %>%
  arrange(cohort)

cohort_dose_path_plot <- ggplot(
  cohort_dose_path,
  aes(x = cohort, y = dose_index)
) +
  geom_point(size = 2, color = "#1f78b4") +
  scale_x_continuous(breaks = seq(1, ncohort, by = 5)) +
  scale_y_continuous(
    breaks = 0:(K - 1),
    labels = paste0(dose_label),
    limits = c(0, K - 1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    x = "Cohort",
    y = "Dose level",
    title = "Dose movement across cohorts"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15)
  )

print(cohort_dose_path_plot)
