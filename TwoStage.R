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

phi <- 0.33
cp <- 0.8
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

posterior_mean_hazard_ratio <- function(draws, k) {
  mean(exp(draws[[paste0("log_hr[", k, "]")]]))
}

treat_cohort <- function(dose_index, survival_outcomes, censor_status, toxicity_indicator) {
  beta_c <- beta.true[dose_index]
  true_times <- rexp(cohortsize, base_hazard * exp(beta_c))
  censoring_times <- rep(followup_time, cohortsize)
  observed_times <- pmin(true_times, censoring_times)
  toxicity <- rbinom(cohortsize, 1, true.toxicity.rate[dose_index])
  
  survival_outcomes[[dose_index]] <- c(survival_outcomes[[dose_index]], observed_times)
  censor_status[[dose_index]] <- c(censor_status[[dose_index]], as.numeric(true_times <= censoring_times))
  toxicity_indicator[[dose_index]] <- c(toxicity_indicator[[dose_index]], toxicity)
  
  list(
    survival_outcomes = survival_outcomes,
    censor_status = censor_status,
    toxicity_indicator = toxicity_indicator,
    toxicity_count = sum(toxicity)
  )
}

set.seed(123)

for (simu in 1:number_of_simulations){
  cidx = 1
  
  survival_outcomes <- vector("list", K)
  censor_status <- vector("list", K)
  toxicity_indicator <- vector("list", K)
  
  ncohort_stage1 = 0
  
  for (i in 1:ncohort){
    ncohort_stage1 = ncohort_stage1 + 1
    
    cohort_result <- treat_cohort(cidx, survival_outcomes, censor_status, toxicity_indicator)
    survival_outcomes <- cohort_result$survival_outcomes
    censor_status <- cohort_result$censor_status
    toxicity_indicator <- cohort_result$toxicity_indicator
    toxicity_count <- cohort_result$toxicity_count
    
    if (toxicity_count == 0 && cidx < K) {
      cidx <- cidx + 1
    } else {
      if (toxicity_count >= 2) {
        cidx <- cidx - 1
      }
      break
    }
  }
  
  max_stage1_dose <- max(cidx, 0)
  ncohort_stage2 = ncohort - ncohort_stage1
  
  if (ncohort_stage2 > 0 && max_stage1_dose > 0) {
    stage2_cohorts_vector <- rep(1:max_stage1_dose, length.out = ncohort_stage2)
    stage2_cohorts_vector <- sample(stage2_cohorts_vector)
    
    for (stage2_dose in stage2_cohorts_vector){
      cohort_result <- treat_cohort(stage2_dose, survival_outcomes, censor_status, toxicity_indicator)
      survival_outcomes <- cohort_result$survival_outcomes
      censor_status <- cohort_result$censor_status
      toxicity_indicator <- cohort_result$toxicity_indicator
    }
  }
  
  dose_eliminated <- as.integer(1:K > max_stage1_dose)
  
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
    phi = phi,
    cp = cp,
    xk = xk,
    mk = mk,
    dose_eliminated = dose_eliminated
  )
  
  fit_mcmc <- mod$sample(
    data = stan_data,
    seed = 300000 + simu,
    chains = 4,
    parallel_chains = 4,
    refresh = 0, 
    iter_sampling = 5000,
    iter_warmup = 2000, 
    thin = 1
  )
  
  params <- fit_mcmc$draws(format = "df")
  
  selected_dose <- select_by_posterior_mode(params)
  
  final_unsafe_prob <- sapply(1:K, function(k) mean(params[[paste0("unsafe_dose[", k, "]")]] == 1))
  admissible_doses <- which(dose_eliminated == 0 & final_unsafe_prob <= 0.5)
  
  MND_best_dose = 0
  if (length(admissible_doses) > 0) {
    mean_hazard_ratio <- rep(NA_real_, K)
    mean_hazard_ratio[admissible_doses] <- sapply(admissible_doses, function(k) posterior_mean_hazard_ratio(params, k))
    
    MND <- rep(Inf, K)
    MND[admissible_doses] <- mean_hazard_ratio[admissible_doses] * d[admissible_doses]^lambda
    
    min_hazard_ratio <- min(mean_hazard_ratio[admissible_doses], na.rm = TRUE)
    noninferior_doses <- admissible_doses[
      mean_hazard_ratio[admissible_doses] <= eta * min_hazard_ratio
    ]
    
    if (length(noninferior_doses) > 0) {
      MND_best_dose = noninferior_doses[which.min(MND[noninferior_doses])]
    }
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
