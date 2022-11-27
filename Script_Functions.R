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
psm_weights <- function(outcome,predictors,id_var,data) {
### Get the model with highest accuracy 
model <- model_selection(outcome,predictors,id_var,data)
### Get the Data set with the chosen variables 
data_build_psm <- subset(data_build(outcome, predictors, id_var, data), 
                        select = c(id_var, "outcome", model[[1]]))
### Propensity Score Model
reg_form <- as.formula(paste("outcome",
                             paste(model[[1]], collapse = " + "), 
                             sep = " ~ "))
psm_model <- glm(reg_form, 
                 data = data_build_psm, 
                 family = binomial(link = "logit"))
# Generate propensity scores and IPWs
psm_weights <- data_build_psm %>% mutate(propensity = psm_model$fitted.values) %>% 
                        mutate(ipw = ifelse(outcome == 1, 1/propensity, 1/(1-propensity)))
### Return Data Frame with Fips and IPW weights 
psm_results <- list(model, psm_weights, psm_model)
return(psm_results)
}

################################################################################
### Function 5: Model Selection for Propensity Score Matching  
################################################################################
model_selection <- function(outcome,predictors,id_var,data){
### Build data set with organized variables. This data set should only include the chosen predictors
### the predictors squared, interactions across predictors, and the outcome 
train_data <- data_build(outcome,predictors,id_var,data)
### Initial Set of Predictors
predictors_i <- colnames(train_data)[which(colnames(train_data) != "outcome")] 
pred_no <- length(predictors_i) - 1
predictor_star <- NULL
kappa_master  <- matrix(data = NA, nrow = pred_no, ncol = 3)
models_master <- list()
### Counters 
delta_kappa = 0.00001
kappa_old = 0 
j = 1
while(delta_kappa > 0) {
models <- list()
kappa <- matrix(data = NA, nrow = pred_no, ncol = 3)
### Find the variable with highest stand alone prediction accuracy
for(i in 1:pred_no) {
### Model Definition 
data_subset <- subset(train_data, select = c("outcome",
                                             predictor_star,
                                             predictors_i[i]))
### Estimate the Logit Models Using Cross Validation
model_cv <-  train(outcome ~ . , data = data_subset, 
                          trControl = cv_options, 
                          method = "glm", 
                          family = binomial(link = "logit"))
### Store the Results 
kappa[i,1] <- i
### Store the Accuracy of the model. 
kappa[i,2] <- model_cv$results$Kappa
### Compute the accuracy improvement relative to the last iteration
kappa[i,3] <- model_cv$results$Kappa - kappa_old
models[[i]] <- model_cv
}
### Update the variables for the other lap of the loop
kappa_master[j,1] <- j
kappa_master[j,2] <- max(kappa[,2])
kappa_master[j,3] <- max(kappa[,3])
models_master[[j]] <- models[[which(kappa[,2] == max(kappa[,2]))]]
### Choose the predictor with highest accuracy
predictor_star <- c(predictor_star, 
                    predictors_i[kappa[which(kappa[,2] == max(kappa[,2])),1]])
### New Predictor Vector 
predictors_i <- predictors_i[-kappa[which(kappa[,2] == max(kappa[,2])),1]]
pred_no <- length(predictors_i) - 1
### Accuracy Improvement
kappa_old <- max(kappa[,2])
delta_kappa <- max(kappa[,3])
### Add one to the counter 
j = j + 1
}
### Remove the last predictor (i.e. the one that derived in decreasing the accuracy)
predictor_star <- predictor_star[-length(predictor_star)]
results <- list(predictor_star,kappa_master,models_master)
return(results)
}


################################################################################

### Function 6: Build Data set To Estimate Iteratively the Propensity Score Model 
################################################################################
data_build <- function(outcome,predictors,id_var,data){
  ### Partition the data set in order to do the cross validation 
  full_data <- data %>% na.omit()
  ### Rename outcome
  names(full_data)[which(colnames(full_data) == outcome)] <- "outcome"
  full_data$outcome <- as.factor(full_data$outcome)
  ### Do the partition between training and testing data
  partition_dataset <- createDataPartition(y = full_data$outcome, p = 0.5 , list = FALSE, times = 1)
  training <- full_data[partition_dataset,]
  ### id_var
  id_var <- subset(full_data, select = id_var)
  ### Outcome has the name of the outcome variable 
  y <- subset(full_data, select = outcome)
  ### Predictors is a vector with the names of all the variables that will be used as predictors 
  x <- subset(full_data, select = predictors)
  #### Squared Predictors 
  x2 <- x^2
  colnames(x2) <- paste0(colnames(x),'.sq')
  ### Interaction Terms and Intercept  
  interactions <- model.matrix( ~ .^2, data = x)
  inter_var <- interactions[,-which(colnames(interactions) == "(Intercept)")] %>% as.data.frame()
  ### Weird thing happening with the variable names when there is only one predictor 
  if(length(predictors) == 1){
    names(inter_var) <- predictors
  }
  ### Return Dataframe with ordered variables 
  regdata <- data.frame(id_var, y, inter_var, x2)
  return(regdata)
}
################################################################################



