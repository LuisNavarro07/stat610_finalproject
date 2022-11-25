################################################################################
################################################################################
### STAT-610 
### Final Project 
### Matched Difference-in-Difference 
### Author: Luis Navarro
### This Version: November 24, 2022
### Script: Functions for the Implementation 
################################################################################
################################################################################

### Function 1: Randomization Inference 
################################################################################
### Randomization Inference: Estimate the Regression Model (with the weights fixed), by create a fake treatment variable
### that randomly assigns the treatment. 
did_randomization <- function(outcome,treat_var,treat_period,time_var,unit_var,weights,data,samples){
### Baseline Estimation
baseline_model <- did_model(outcome,treat_var,treat_period,time_var,unit_var,weights,data)
### Randomization Inference
### Empty Data frame to store the results 
empirical_distribution <- matrix(nrow = samples, ncol = 1)
i  = 1
while (i <= samples){
  ### Create Random Sample 
  data_sample <- random_treatment(unit_var,data)
  ### Estimate the model with random treatment variable 
  model_estimate <- did_model(outcome = outcome, 
                              treat_var = "random_treat",
                              treat_period = treat_period, 
                              time_var = time_var, 
                              unit_var = unit_var, 
                              weights = weights, 
                              data = data_sample)
  empirical_distribution[i,1] <- coef(model_estimate)
  i = i + 1
}
### Average Treatment Effect 
ate <- coef(baseline_model)
### Standard Error: standard deviation from the empirical distribution 
std_error <- sd(empirical_distribution)
### P-value: proportion of times the placebo treatment effect is larger than the estimated treatment effect 
pval <- sum(empirical_distribution >= ate)/samples
model_results <- data.frame(coefficient = ate, std_error = std_error, pval = pval)
return(model_results)
}

### Function 2: Estimate the Difference-in-Difference Model 
################################################################################
### Difference in difference regression model 
did_model <- function(outcome,treat_var,treat_period,time_var,unit_var,weights,data){
  ### Build the Variables to run the Difference-in-Difference Model
  ### Rename Variables For Simplicity
  data_model <- data
  ### Rename variables for consistency 
  names(data_model)[which(colnames(data_model) == outcome)] <- "outcome"
  names(data_model)[which(colnames(data_model) == treat_var)] <- "treat_var"
  names(data_model)[which(colnames(data_model) == time_var)] <- "time_var"
  names(data_model)[which(colnames(data_model) == unit_var)] <- "unit_var"
  ### If missing weights 
  if(is.null(weights)) {
    uniform_weights <- rep(1, length(data_model$unit_var))
    data_model <- mutate(data_model, weights = uniform_weights)
  } else {
    names(data_model)[which(colnames(data_model) == weights)] <- "weights"
  }
  ### Create Variables to Run Difference-in-Difference Model 
  data_model <- mutate(data_model, time_post = case_when(data_model$time_var <= treat_period ~ 0,
                                                         data_model$time_var > treat_period ~ 1),
                       did_var = treat_var*time_post)

  ### Create Post Variable for the Difference - in - Difference Model
  data_model <- mutate(data_model, time_post = case_when(data_model$time_var <= treat_period ~ 0,
                                                         data_model$time_var > treat_period ~ 1),
                       did_var = treat_var*time_post)
  ### Estimate the Regression Model with IP weights 
  model_did <- feols(outcome ~ did_var | unit_var + time_var, data = data_model, weights = ~ weights) 
  return(model_did)
}
################################################################################

### Function 3: Create a Random Assignment in Treatment 
################################################################################
random_treatment <- function(unit_var,data){
  data_treat <- data 
  names(data_treat)[which(colnames(data_treat) == unit_var)] <- "unit_var"
  ### Get A list of unique treatment units
  treat_units <- unique(data_treat$unit_var)
  ### Create a fake random treatment 
  random_treat <- round(runif(length(treat_units)))
  ### Map of Randomized Treatment
  random_assignment <- data.frame(unit_var = treat_units, random_treat = random_treat)
  ### Merge with the data set 
  data_treat <- merge(data_treat, random_assignment, by = "unit_var")
  return(data_treat)
}
################################################################################

### Function 4: Estimate the Propensity Score Model   
################################################################################
psm_weights <- function(outcome,predictors,data) {
### Get the model with highest accuracy 
model <- model_selection(outcome,predictors,data)
model_pred <- c("outcome",model[[1]])
### Get the Data set with the chosen variables 
data_build_psm <- data_build(outcome,predictors,data)
data_build_psm <- subset(data_build_psm, select = model_pred)
### Estimate the Inverse Probability Weights 
### Propensity Score Model: Numerator
psm_model_num <- glm(outcome ~ 1, 
                       data = data_build_psm, 
                       family = "binomial")
### Propensity Score Model: Denominator
psm_model_den <- glm(outcome ~ ., 
                       data = data_build_psm, 
                       family = "binomial")
### Compute Propensity Score Weights
psm_weights <- data %>% mutate(prop_num = psm_model_num$fitted.values, 
                                     prop_den = psm_model_den$fitted.values,
                                     prop_num_out = ifelse(treat_det == 1, prop_num, 1 - prop_num),
                                     prop_den_out = ifelse(treat_det == 1, prop_den, 1 - prop_den),
                                     ipw = prop_num_out/prop_den_out) 
psm_weights <- data.frame(fips = psm_weights$fips, ipw = psm_weights$ipw)  
### Return Data Frame with Fips and IPW weights 
psm_results <- list(model, psm_weights)
return(psm_results)

}

################################################################################
### Function 5: Model Selection for Propensity Score Matching  
################################################################################
model_selection <- function(outcome,predictors,data){
### Build data set with organized variables
train_data <- data_build(outcome,predictors,data)
### Run a loop that will estimate 
models <- list()
kappa <- matrix(data = NA, nrow = ncol(train_data)-1, ncol = 2)
#### Estimate all Models in a for loop 
for(i in 2:ncol(train_data)) {
  ### Create a subset with only the outcome and the predictors for each lap 
  predictors_i <- colnames(train_data)[2:i] 
  data_subset <- subset(train_data, select = c("outcome",predictors_i))
  ### Estimate the Logit Models Using Cross Validation
  models[[i - 1]] <-  train(outcome ~ . , data = data_subset, 
                           trControl = cv_options, 
                           method = "glm", 
                           family = "binomial")
  ### Store the Results 
  kappa[i-1,1] <- i -1 
  kappa[i-1,2] <- models[[i - 1]]$results$Kappa
}
### Choose the one with highest accuracy
max_kappa <- max(kappa[,2])
best_model <- models[[kappa[which(kappa[,2] == max_kappa),1]]]
### Store Best Model Predictors 
best_model_predictors <- colnames(best_model$finalModel$data)
best_model_predictors <- best_model_predictors[-which(best_model_predictors == ".outcome")]
results <- list(best_model_predictors,kappa)
return(results)
}
################################################################################

### Function 6: Build Data set To Estimate Iteratively the Propensity Score Model 
################################################################################
data_build <- function(outcome,predictors,data){
  ### Partition the data set in order to do the cross validation 
  full_data <- data %>% na.omit()
  ### Rename outcome
  names(full_data)[which(colnames(full_data) == outcome)] <- "outcome"
  full_data$outcome <- as.factor(full_data$outcome)
  ### Do the partition between training and testing data
  partition_dataset <- createDataPartition(y = full_data$outcome, p = 0.5 , list = FALSE, times = 1)
  training <- full_data[partition_dataset,]
  ### Outcome has the name of the outcome variable 
  y <- subset(full_data, select = outcome)
  ### Predictors is a vector with the names of all the variables that will be used as predictors 
  x <- subset(full_data, select = predictors)
  #### Squared Predictors 
  x2 <- x^2
  colnames(x2) <- paste0(colnames(x),'.sq')
  ### Interaction Terms 
  interactions <- model.matrix( ~ .^2, data = x)
  interactions <- as.data.frame(interactions[,-which(colnames(interactions) == "(Intercept)")])
  ### Return Dataframe with ordered variables 
  regdata <- data.frame(y, x, x2, interactions)
  return(regdata)
}
################################################################################



