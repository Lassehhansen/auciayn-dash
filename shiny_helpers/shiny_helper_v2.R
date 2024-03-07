library(purrr)
library(dplyr)
library(yardstick)

set.seed(1997)


convert_ev_confidence_to_alpha_beta <- function(ev, confidence) {
  alpha <- ev * confidence
  beta <- (1 - ev) * confidence
  return(list(alpha = alpha, beta = beta))
}

# Function to simulate a subgroup
simulate_subgroup <- function(total_samples, alpha_pos, beta_pos, alpha_neg, beta_neg, prevalence) {
  
  n_pos <- round(total_samples * prevalence)
  n_neg <- total_samples - n_pos
  
  # Simulate positives
  positives <- rbeta(n_pos, alpha_pos, beta_pos)
  
  # Simulate negatives
  negatives <- rbeta(n_neg, alpha_neg, beta_neg)
  
  probability <- c(positives, negatives)
  labels <- c(rep(1, n_pos), rep(0, n_neg))
  
  return(list(probability = probability, labels = labels))
}

calculate_metrics <- function(dataset) {
  # Ensure labels are a factor and the levels are correctly ordered to reflect the event of interest
  dataset$labels <- factor(dataset$labels, levels = c(1, 0))
  
  # Calculate AUROC
  auroc_value <- yardstick::roc_auc(dataset, labels, probability)$.estimate
  
  # Calculate AUPRC
  auprc_value <- yardstick::pr_auc(dataset, labels, probability)$.estimate
  
  brier_value <- yardstick::brier_class(dataset, labels, probability)$.estimate
  
  return(list(auroc = auroc_value, auprc = auprc_value, brier = brier_value))
}

sample_population <- function(total_samples, attribute_ratio, 
                              alpha_pos_a, beta_pos_a, alpha_neg_a, beta_neg_a,
                              alpha_pos_b, beta_pos_b, alpha_neg_b, beta_neg_b,
                              p_pos_space) {
  
  # Total samples for simulation
  total_samples <- total_samples # Adjust as needed
  
  # Calculate samples for each subgroup based on the specified ratio
  total_samples_a <- round(total_samples * attribute_ratio)
  total_samples_b <- total_samples - total_samples_a
  
  # Simulate subgroup A
  sim_a <- simulate_subgroup(total_samples_a, alpha_pos_a, beta_pos_a, alpha_neg_a, beta_neg_a, p_pos_space)
  
  # Simulate subgroup B
  sim_b <- simulate_subgroup(total_samples_b, alpha_pos_b, beta_pos_b, alpha_neg_b, beta_neg_b, p_pos_space)
  
  # Merge the simulated data
  data_combined <- c(sim_a$probability, sim_b$probability)
  labels_combined <- c(sim_a$labels, sim_b$labels)
  subgroup_combined <- c(rep("a", length(sim_a$probability)), rep("b", length(sim_b$probability)))
  
  dataset_total <- data.frame(
    probability = data_combined,
    labels = as.factor(labels_combined),
    subgroup = subgroup_combined
  )
  
  return(dataset_total)
}


# Adjusted get_metrics function to include parameters
get_metrics <- function(ev_pos_a, ev_pos_b, confidence_pos_a, confidence_pos_b, 
                        confidence_neg_a, confidence_neg_b, ev_neg_a, ev_neg_b, 
                        p_pos_space, attribute_ratio, total_samples) {
  # Convert expected values and confidence to alpha and beta
  pdf_pos_a_params <- convert_ev_confidence_to_alpha_beta(ev_pos_a, confidence_pos_a)
  pdf_neg_a_params <- convert_ev_confidence_to_alpha_beta(ev_neg_a, confidence_neg_a)
  
  pdf_pos_b_params <- convert_ev_confidence_to_alpha_beta(ev_pos_b, confidence_pos_b)
  pdf_neg_b_params <- convert_ev_confidence_to_alpha_beta(ev_neg_b, confidence_neg_b)
  
  # Sample population
  dataset_total <- sample_population(total_samples, attribute_ratio, 
                                     ev_pos_a, ev_neg_a, ev_pos_b, ev_neg_b, 
                                     confidence_pos_a, confidence_neg_a, confidence_pos_b, confidence_neg_b, 
                                     p_pos_space)
  
  # Calculate metrics
  dataset_a <- filter(dataset_total, subgroup == "a")
  dataset_b <- filter(dataset_total, subgroup == "b")
  
  metrics_a <- calculate_metrics(dataset_a)
  metrics_b <- calculate_metrics(dataset_b)
  metrics_total <- calculate_metrics(dataset_total)
  
  # Create a dataframe for results including parameters
  data <- data.frame(
    auroc_a = metrics_a$auroc,
    auprc_a = metrics_a$auprc,
    auroc_b = metrics_b$auroc,
    auprc_b = metrics_b$auprc,
    auroc_total = metrics_total$auroc,
    auprc_total = metrics_total$auprc,
    brier_a = metrics_a$brier,
    brier_b = metrics_b$brier,
    brier_total = metrics_total$brier,
    
    alpha_pos_a = pdf_pos_a_params$alpha,
    beta_pos_a = pdf_pos_a_params$beta,
    alpha_neg_a = pdf_neg_a_params$alpha,
    beta_neg_a = pdf_neg_a_params$beta,
    alpha_pos_b = pdf_pos_b_params$alpha,
    beta_pos_b = pdf_pos_b_params$beta,
    alpha_neg_b = pdf_neg_b_params$alpha,
    beta_neg_b = pdf_neg_b_params$beta,
    p_pos_space = p_pos_space,
    attribute_ratio = attribute_ratio,
    total_samples, total_samples,
    ev_pos_a = ev_pos_a, 
    ev_pos_b = ev_pos_b, 
    confidence_pos_a = confidence_pos_a, 
    confidence_pos_b = confidence_pos_b,
    ev_neg_a = ev_neg_a, 
    ev_neg_b = ev_neg_b, 
    confidence_neg_a = confidence_neg_a, 
    confidence_neg_b = confidence_neg_b
  )
  
  return(data)
}