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
eta = 1.2
lambda = 0.7

true.toxicity.rate = c(0.01,0.03,0.05,0.07,0.10,0.12,0.15)

base_hazard <- 1

followup_time = 50

gamma_L = 1
gamma_R = 1
prob_cutoff = 0.5
K = length(beta.true)

number_of_patients = rep(0, K)
dose_selection = rep(0, K)
dose_selection_MND = rep(0, K)

file <- file.path("Stan Code Partial Likelihood.stan")
mod <- cmdstan_model(file)

posterior_best_dose <- function(draws) {
  best_dose <- draws[["best_dose"]]
  best_dose[best_dose > 0]
}

select_by_posterior_mode <- function(draws) {
  best_dose <- posterior_best_dose(draws)
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

set.seed(123)

for (simu in 1:number_of_simulations){
  cidx = 1
  
  survival_outcomes <- vector("list", K)
  censor_status <- vector("list", K)
  toxicity_indicator <- vector("list", K)
  dose_eliminated = rep(0,K)
  
  for (i in 1:ncohort){
    beta_c <- beta.true[cidx]
    true_times <- rexp(cohortsize, base_hazard*exp(beta_c))
    censoring_times <- rep(followup_time, cohortsize)
    t = pmin(true_times, censoring_times)
    toxicity = rbinom(cohortsize,1,true.toxicity.rate[cidx])
    survival_outcomes[[cidx]] = c(survival_outcomes[[cidx]],t)
    censor_status[[cidx]] = c(censor_status[[cidx]], as.numeric(true_times <= censoring_times))
    toxicity_indicator[[cidx]] = c(toxicity_indicator[[cidx]], toxicity)
    
    dose_level = c()
    for (k in 1:length(beta.true)){dose_level = c(dose_level, rep(k,length(survival_outcomes[[k]])))}
    t_array = c()
    for (k in 1:length(beta.true)){t_array = c(t_array, survival_outcomes[[k]])}
    c_array = c()
    for (k in 1:length(beta.true)){c_array = c(c_array, censor_status[[k]])}
    xk = c()
    for (k in 1:length(beta.true)){xk = c(xk, sum(toxicity_indicator[[k]]))}
    mk = c()
    for (k in 1:length(beta.true)){mk = c(mk, length(toxicity_indicator[[k]]))}
    
    stan_data <- list(
      N = length(t_array),
      K = K,
      dose = dose_level,
      time = t_array,
      event = c_array,
      phi = 0.33,
      cp = 0.8,
      xk = xk,
      mk = mk,
      dose_eliminated = dose_eliminated
    )
    
    fit_mcmc <- mod$sample(
      data = stan_data,
      seed = 200000 + 100 * simu + i,
      chains = 4,
      parallel_chains = 4,
      refresh = 0, 
      iter_sampling = 5000,
      iter_warmup = 2000, 
      thin = 1
    )
    
    # Extract posterior samples
    params <- fit_mcmc$draws(format = "df")
    
    unsafe_prob <- sapply(1:K, function(k) mean(params[[paste0("unsafe_dose[", k, "]")]] == 1))
    unsafe_dose_positions <- which(unsafe_prob > 0.5)
    if (length(unsafe_dose_positions) > 0) {
      first_unsafe_dose <- min(unsafe_dose_positions)
      dose_eliminated[first_unsafe_dose:K] = 1
    }
    eliminated_dose_positions <- which(dose_eliminated == 1)
    if (length(eliminated_dose_positions) > 0) {
      max_feasible_dose_positions = min(eliminated_dose_positions)-1
    } else {
      max_feasible_dose_positions = K
    }
    
    if (cidx>max_feasible_dose_positions){cidx=max_feasible_dose_positions}
    if (max_feasible_dose_positions==0){break}
    
    best_dose_draws <- posterior_best_dose(params)
    if (length(best_dose_draws) == 0) {break}
    neighbor_doses <- intersect(c(cidx - 1, cidx, cidx + 1), 1:max_feasible_dose_positions)
    psi <- mean(best_dose_draws %in% neighbor_doses)
    
    if (runif(1,0,1)<psi){
      if (cidx > 1 && cidx < max_feasible_dose_positions) {
        # For doses not at boundaries
        hr_left_samples <- hr_adj_draws(params, cidx - 1)
        hr_right_samples <- hr_adj_draws(params, cidx)
        
        # Calculate probabilities
        prob_left_gt <- mean(hr_left_samples > gamma_L)
        prob_left_lt <- mean(hr_left_samples < gamma_L)
        prob_right_gt <- mean(hr_right_samples > gamma_R)
        prob_right_lt <- mean(hr_right_samples < gamma_R)
        
        # Apply decision rules
        if (prob_left_gt > prob_cutoff && prob_right_gt > prob_cutoff) {
          next_dose <- cidx - 1
          decision_type <- "De-escalation"
        } else if (prob_left_lt > prob_cutoff && prob_right_lt > prob_cutoff) {
          next_dose <- cidx + 1
          decision_type <- "Escalation"
        } else if (prob_left_gt > prob_cutoff && prob_right_lt > prob_cutoff) {
          # Random decision case one
          left_ratio <- posterior_mean_adj_hr(params, cidx - 1)
          right_ratio <- 1 / posterior_mean_adj_hr(params, cidx)  # exp(beta_k)/exp(beta_k+1)
          
          deescalate_prob <- left_ratio / (left_ratio + right_ratio)
          
          if (runif(1) < deescalate_prob) {
            next_dose <- cidx - 1
            decision_type <- "Random de-escalation"
          } else {
            next_dose <- cidx + 1
            decision_type <- "Random escalation"
          }
        } else {
          # Random decision case two
          left <- posterior_mean_inverse_hazard(params, cidx - 1)
          stay <- posterior_mean_inverse_hazard(params, cidx)
          right <- posterior_mean_inverse_hazard(params, cidx + 1)
          
          deescalate_prob <- left / (left + stay + right)
          escalate_prob <- right / (left + stay + right)
          
          random_number = runif(1)
          if (random_number < deescalate_prob) {
            next_dose <- cidx - 1
            decision_type <- "Random de-escalation"
          } else if (random_number > 1-escalate_prob){
            next_dose <- cidx + 1
            decision_type <- "Random escalation"
          } else {
            next_dose <- cidx
            decision_type <- "Random stay"
          }
        }
        
      } else if (cidx == 1) {
        # At lowest dose
        hr_right_samples <- hr_adj_draws(params, cidx)
        prob_right_lt <- mean(hr_right_samples < gamma_R)
        
        if (prob_right_lt > prob_cutoff) {
          next_dose <- 2
          decision_type <- "Escalation"
        } else {
          # Random decision case two
          stay <- posterior_mean_inverse_hazard(params, cidx)
          right <- posterior_mean_inverse_hazard(params, cidx + 1)
          
          escalate_prob <- right / (stay + right)
          
          random_number = runif(1)
          if (random_number > 1-escalate_prob){
            next_dose <- cidx + 1
            decision_type <- "Random escalation"
          } else {
            next_dose <- cidx
            decision_type <- "Random stay"
          }
        }
        
      } else { # cidx == max_feasible_dose_positions
        # At highest dose
        hr_left_samples <- hr_adj_draws(params, cidx - 1)
        prob_left_gt <- mean(hr_left_samples > gamma_L)
        
        if (prob_left_gt > prob_cutoff) {
          next_dose <- cidx - 1
          decision_type <- "De-escalation"
        } else {
          # Random decision case two
          left <- posterior_mean_inverse_hazard(params, cidx - 1)
          stay <- posterior_mean_inverse_hazard(params, cidx)
          
          deescalate_prob <- left / (left + stay)
          
          random_number = runif(1)
          if (random_number < deescalate_prob) {
            next_dose <- cidx - 1
            decision_type <- "Random de-escalation"
          } else {
            next_dose <- cidx
            decision_type <- "Random stay"
          }
        }
      }
    } else{
      next_dose = sample(best_dose_draws, 1)
    }
    
    if (max_feasible_dose_positions==1){next_dose=1}
    cidx = next_dose
    
    if (next_dose==0){break}
  }
  
  selected_dose <- select_by_posterior_mode(params)
  
  # minimum noninferiority dose
  MND = sapply(1:K, function(k) posterior_mean_hazard_ratio(params, k) * d[k]^lambda)
  
  mean_hazard_ratio = sapply(1:K, function(k) posterior_mean_hazard_ratio(params, k))
  
  noninferior_indicator = rep(0, K)
  for (i in 1:K){if (mean_hazard_ratio[i] <= eta * min(mean_hazard_ratio)){noninferior_indicator[i] = 1}}
  MND_best_dose = 0
  if (any(noninferior_indicator == 1)) {
    noninferior_doses <- which(noninferior_indicator == 1)
    MND_best_dose = noninferior_doses[which.min(MND[noninferior_doses])]
  }
  
  for (i in 1:K){number_of_patients[i] = number_of_patients[i] + length(survival_outcomes[[i]])}
  if (selected_dose > 0) {dose_selection[selected_dose] = dose_selection[selected_dose]+1}
  if (MND_best_dose > 0) {dose_selection_MND[MND_best_dose] = dose_selection_MND[MND_best_dose]+1}
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

print(simulation_summary)


