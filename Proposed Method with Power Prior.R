setwd("C:/Users/y277h/Desktop/CFH Project/My Code and Results")

library(cmdstanr)
library(StanHeaders)
library(rlist)
library(survival)
library(dplyr)
library(ggplot2)

ncohort <- 20
cohortsize <- 4
number_of_simulations <- 500

beta.true <- c(0, -1, -2.7, -3, -2, -1, 0)
d <- c(1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6)
eta <- 1.2
lambda <- 0.7

true.toxicity.rate <- c(0.01, 0.03, 0.05, 0.07, 0.10, 0.12, 0.15)

base_hazard <- 1
followup_time <- 50

gamma_L <- 1
gamma_R <- 1
prob_cutoff <- 0.5
phi <- 0.33
cp <- 0.8
K <- length(beta.true)

a0_fixed <- 0.5
use_historical <- 1

beta.historical.true <- beta.true
toxicity.rate.historical <- true.toxicity.rate
n_historical_per_dose <- 10

file <- file.path("Stan Code Power Prior.stan")
mod <- cmdstan_model(file)

posterior_best_dose <- function(draws, max_feasible_dose = K) {
  best_dose <- draws[["best_dose"]]
  best_dose[best_dose > 0 & best_dose <= max_feasible_dose]
}

select_by_posterior_mode <- function(draws, max_feasible_dose = K) {
  best_dose <- posterior_best_dose(draws, max_feasible_dose)
  if (length(best_dose) == 0) return(0)
  best_dose_table <- table(best_dose)
  as.integer(names(best_dose_table)[which.max(best_dose_table)])
}

hr_adj_draws <- function(draws, k) {
  draws[[paste0("hr_adj[", k, "]")]]
}

posterior_mean_adj_hr <- function(draws, k) {
  exp(mean(draws[[paste0("adj_log_hr[", k, "]")]]))
}

posterior_mean_hazard_ratio <- function(draws, k) {
  exp(mean(draws[[paste0("log_hr[", k, "]")]]))
}

posterior_mean_inverse_hazard <- function(draws, k) {
  1 / posterior_mean_hazard_ratio(draws, k)
}

posterior_mean_log_hr <- function(draws) {
  sapply(1:K, function(k) mean(draws[[paste0("log_hr[", k, "]")]]))
}

generate_historical_data <- function(n_patients_per_dose,
                                     beta.historical,
                                     toxicity.rate.historical,
                                     base_hazard_hist = base_hazard) {
  if (length(beta.historical) != K) {
    stop("beta.historical must have length K.")
  }
  if (length(toxicity.rate.historical) != K) {
    stop("toxicity.rate.historical must have length K.")
  }
  
  historical_survival <- vector("list", K)
  historical_censor <- vector("list", K)
  historical_toxicity <- vector("list", K)
  
  for (k in 1:K) {
    true_times <- rexp(n_patients_per_dose, base_hazard_hist * exp(beta.historical[k]))
    censoring_times <- rep(followup_time, n_patients_per_dose)
    historical_survival[[k]] <- pmin(true_times, censoring_times)
    historical_censor[[k]] <- as.numeric(true_times <= censoring_times)
    historical_toxicity[[k]] <- rbinom(n_patients_per_dose, 1, toxicity.rate.historical[k])
  }
  
  list(
    N0 = sum(sapply(historical_survival, length)),
    dose0 = unlist(lapply(1:K, function(k) rep(k, length(historical_survival[[k]]))), use.names = FALSE),
    time0 = unlist(historical_survival, use.names = FALSE),
    event0 = unlist(historical_censor, use.names = FALSE),
    x0k = sapply(historical_toxicity, sum),
    m0k = sapply(historical_toxicity, length),
    beta_true = beta.historical
  )
}

update_discount_parameter <- function(current_log_hr, historical_log_hr, iteration) {
  discrepancy <- mean(abs(current_log_hr - historical_log_hr), na.rm = TRUE)
  a0_new <- a0_fixed * exp(-discrepancy) * (1 - min(iteration / ncohort, 0.7))
  min(max(a0_new, 0), a0_fixed)
}

number_of_patients <- rep(0, K)
dose_selection <- rep(0, K)
dose_selection_MND <- rep(0, K)
a0_trace <- c()

historical_data <- generate_historical_data(
  n_patients_per_dose = n_historical_per_dose,
  beta.historical = beta.historical.true,
  toxicity.rate.historical = toxicity.rate.historical
)

set.seed(123)

for (simu in 1:number_of_simulations) {
  cidx <- 1
  a0_current <- a0_fixed
  
  survival_outcomes <- vector("list", K)
  censor_status <- vector("list", K)
  toxicity_indicator <- vector("list", K)
  dose_eliminated <- rep(0, K)
  
  for (i in 1:ncohort) {
    beta_c <- beta.true[cidx]
    true_times <- rexp(cohortsize, base_hazard * exp(beta_c))
    censoring_times <- rep(followup_time, cohortsize)
    observed_times <- pmin(true_times, censoring_times)
    toxicity <- rbinom(cohortsize, 1, true.toxicity.rate[cidx])
    
    survival_outcomes[[cidx]] <- c(survival_outcomes[[cidx]], observed_times)
    censor_status[[cidx]] <- c(censor_status[[cidx]], as.numeric(true_times <= censoring_times))
    toxicity_indicator[[cidx]] <- c(toxicity_indicator[[cidx]], toxicity)
    
    dose_level <- unlist(lapply(1:K, function(k) rep(k, length(survival_outcomes[[k]]))), use.names = FALSE)
    t_array <- unlist(survival_outcomes, use.names = FALSE)
    c_array <- unlist(censor_status, use.names = FALSE)
    xk <- sapply(toxicity_indicator, sum)
    mk <- sapply(toxicity_indicator, length)
    
    stan_data <- list(
      N = length(t_array),
      K = K,
      dose = dose_level,
      time = t_array,
      event = c_array,
      N0 = historical_data$N0,
      dose0 = historical_data$dose0,
      time0 = historical_data$time0,
      event0 = historical_data$event0,
      phi = phi,
      cp = cp,
      xk = xk,
      mk = mk,
      x0k = historical_data$x0k,
      m0k = historical_data$m0k,
      dose_eliminated = dose_eliminated,
      a0 = a0_current,
      use_historical = use_historical
    )
    
    fit_mcmc <- mod$sample(
      data = stan_data,
      seed = 400000 + 100 * simu + i,
      chains = 4,
      parallel_chains = 4,
      refresh = 0,
      iter_sampling = 5000,
      iter_warmup = 2000,
      thin = 1
    )
    
    params <- fit_mcmc$draws(format = "df")
    
    unsafe_prob <- sapply(1:K, function(k) mean(params[[paste0("unsafe_dose[", k, "]")]] == 1))
    unsafe_dose_positions <- which(unsafe_prob > 0.5)
    if (length(unsafe_dose_positions) > 0) {
      first_unsafe_dose <- min(unsafe_dose_positions)
      dose_eliminated[first_unsafe_dose:K] <- 1
    }
    
    eliminated_dose_positions <- which(dose_eliminated == 1)
    if (length(eliminated_dose_positions) > 0) {
      max_feasible_dose_positions <- min(eliminated_dose_positions) - 1
    } else {
      max_feasible_dose_positions <- K
    }
    
    if (cidx > max_feasible_dose_positions) {
      cidx <- max_feasible_dose_positions
    }
    if (max_feasible_dose_positions == 0) {
      break
    }
    
    if (use_historical == 1 && i > 2) {
      current_log_hr <- posterior_mean_log_hr(params)
      a0_current <- update_discount_parameter(current_log_hr, historical_data$beta_true, i)
    }
    a0_trace <- c(a0_trace, a0_current)
    
    best_dose_draws <- posterior_best_dose(params, max_feasible_dose_positions)
    if (length(best_dose_draws) == 0) {
      break
    }
    neighbor_doses <- intersect(c(cidx - 1, cidx, cidx + 1), 1:max_feasible_dose_positions)
    psi <- mean(best_dose_draws %in% neighbor_doses)
    
    if (runif(1, 0, 1) < psi) {
      if (cidx > 1 && cidx < max_feasible_dose_positions) {
        hr_left_samples <- hr_adj_draws(params, cidx - 1)
        hr_right_samples <- hr_adj_draws(params, cidx)
        
        prob_left_gt <- mean(hr_left_samples > gamma_L)
        prob_left_lt <- mean(hr_left_samples < gamma_L)
        prob_right_gt <- mean(hr_right_samples > gamma_R)
        prob_right_lt <- mean(hr_right_samples < gamma_R)
        
        if (prob_left_gt > prob_cutoff && prob_right_gt > prob_cutoff) {
          next_dose <- cidx - 1
        } else if (prob_left_lt > prob_cutoff && prob_right_lt > prob_cutoff) {
          next_dose <- cidx + 1
        } else if (prob_left_gt > prob_cutoff && prob_right_lt > prob_cutoff) {
          left_ratio <- posterior_mean_adj_hr(params, cidx - 1)
          right_ratio <- 1 / posterior_mean_adj_hr(params, cidx)
          deescalate_prob <- left_ratio / (left_ratio + right_ratio)
          
          if (runif(1) < deescalate_prob) {
            next_dose <- cidx - 1
          } else {
            next_dose <- cidx + 1
          }
        } else {
          left <- posterior_mean_inverse_hazard(params, cidx - 1)
          stay <- posterior_mean_inverse_hazard(params, cidx)
          right <- posterior_mean_inverse_hazard(params, cidx + 1)
          
          deescalate_prob <- left / (left + stay + right)
          escalate_prob <- right / (left + stay + right)
          
          random_number <- runif(1)
          if (random_number < deescalate_prob) {
            next_dose <- cidx - 1
          } else if (random_number > 1 - escalate_prob) {
            next_dose <- cidx + 1
          } else {
            next_dose <- cidx
          }
        }
        
      } else if (cidx == 1) {
        hr_right_samples <- hr_adj_draws(params, cidx)
        prob_right_lt <- mean(hr_right_samples < gamma_R)
        
        if (prob_right_lt > prob_cutoff) {
          next_dose <- 2
        } else {
          stay <- posterior_mean_inverse_hazard(params, cidx)
          right <- posterior_mean_inverse_hazard(params, cidx + 1)
          escalate_prob <- right / (stay + right)
          
          random_number <- runif(1)
          if (random_number > 1 - escalate_prob) {
            next_dose <- cidx + 1
          } else {
            next_dose <- cidx
          }
        }
        
      } else {
        hr_left_samples <- hr_adj_draws(params, cidx - 1)
        prob_left_gt <- mean(hr_left_samples > gamma_L)
        
        if (prob_left_gt > prob_cutoff) {
          next_dose <- cidx - 1
        } else {
          left <- posterior_mean_inverse_hazard(params, cidx - 1)
          stay <- posterior_mean_inverse_hazard(params, cidx)
          deescalate_prob <- left / (left + stay)
          
          random_number <- runif(1)
          if (random_number < deescalate_prob) {
            next_dose <- cidx - 1
          } else {
            next_dose <- cidx
          }
        }
      }
    } else {
      next_dose <- sample(best_dose_draws, 1)
    }
    
    if (max_feasible_dose_positions == 1) {
      next_dose <- 1
    }
    cidx <- next_dose
    
    if (next_dose == 0) {
      break
    }
  }
  
  selected_dose <- select_by_posterior_mode(params, max_feasible_dose_positions)
  MND_best_dose <- 0
  if (max_feasible_dose_positions > 0) {
    admissible_doses <- 1:max_feasible_dose_positions
    MND <- rep(Inf, K)
    mean_hazard_ratio <- rep(NA_real_, K)
    mean_hazard_ratio[admissible_doses] <- sapply(admissible_doses, function(k) posterior_mean_hazard_ratio(params, k))
    MND[admissible_doses] <- mean_hazard_ratio[admissible_doses] * d[admissible_doses]^lambda
    
    min_hazard_ratio <- min(mean_hazard_ratio[admissible_doses], na.rm = TRUE)
    noninferior_doses <- admissible_doses[
      mean_hazard_ratio[admissible_doses] <= eta * min_hazard_ratio
    ]
    
    if (length(noninferior_doses) > 0) {
      MND_best_dose <- noninferior_doses[which.min(MND[noninferior_doses])]
    }
  }
  
  for (i in 1:K) {
    number_of_patients[i] <- number_of_patients[i] + length(survival_outcomes[[i]])
  }
  if (selected_dose > 0) {
    dose_selection[selected_dose] <- dose_selection[selected_dose] + 1
  }
  if (MND_best_dose > 0) {
    dose_selection_MND[MND_best_dose] <- dose_selection_MND[MND_best_dose] + 1
  }
  
  if (simu %% 50 == 0) {
    cat("Simulation", simu, "completed. Current mean a0:", round(mean(tail(a0_trace, 10)), 3), "\n")
  }
}

simulation_summary <- data.frame(
  Dose_Level = 1:K,
  Patients = number_of_patients,
  Allocation_Percent = round(100 * number_of_patients / sum(number_of_patients), 1),
  Selection_Count = dose_selection,
  Selection_Percent = round(100 * dose_selection / number_of_simulations, 1),
  MND_Selection_Count = dose_selection_MND,
  MND_Selection_Percent = round(100 * dose_selection_MND / number_of_simulations, 1)
)

discount_summary <- data.frame(
  Mean_a0 = round(mean(a0_trace), 3),
  Median_a0 = round(median(a0_trace), 3),
  Min_a0 = round(min(a0_trace), 3),
  Max_a0 = round(max(a0_trace), 3)
)

print(simulation_summary)
print(discount_summary)
