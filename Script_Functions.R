################################################################################
################################################################################
### STAT-610 
### Final Project 
### Matched Difference-in-Difference 
### Author: Luis Navarro
### This Version: December 15, 2022
### Script: Functions for the Implementation 
################################################################################
################################################################################

### Function 1: Randomization Inference 
################################################################################
### Randomization Inference: Estimate the Regression Model (with the weights fixed), by create a fake treatment variable
### that randomly assigns the treatment. 
#' @export
did_randomization <- function(outcome,treat_var,treat_period,time_var,unit_var,weights,data,samples){
if(missing(samples)) {
  warning('number of samples not defined. Default is 100')
  samples <- 100 
}
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
ate <- coef(baseline_model) %>% as.double()
### Standard Error: standard deviation from the empirical distribution 
std_error <- sd(empirical_distribution) %>% as.double()
### P-value: proportion of times the placebo treatment effect is larger than the estimated treatment effect 
pval <- sum(empirical_distribution >= ate)/samples %>% as.double()
model_results <- data.frame(coefficient = ate, std_error = std_error, pval = pval)
### Include the Average Treatment Effect as a variable in the empirical distribution function
empirical_distribution <- empirical_distribution %>% as.data.frame() %>% mutate(ATE = ate)
ouput <- list(model_results, empirical_distribution)
return(ouput)
}

################################################################################
### Function 1.1 Show the Empirical Density of the Placebo Distribution 
#' @export
ate_inference_graph <- function(empirical_distribution) {
### Do a graph showing the hypothesis graphically 
### Density Plot of the Empirical Distribution 
density_plot <- ggplot(empirical_distribution, aes(x = V1)) + 
  geom_density(alpha=0.25) + 
  geom_vline(aes(xintercept=mean(ATE)), color = "red", linetype = "dashed") +
  theme_bw() + geom_density(color="darkblue") +
  labs(title = "Placebo Distribution of the Average Treatment Effect", y="Density", x = "Parameter Estimate") + 
  scale_color_manual(name = "Estimate")
print(density_plot)
return(density_plot)
}

### Function 2: Estimate the Difference-in-Difference Model 
################################################################################
### Difference in difference regression model 
#' @export
did_model <- function(outcome,treat_var,treat_period,time_var,unit_var,weights,data){
  ### Stop if inputs are not defined 
  if(missing(outcome)) {stop('outcome variable not defined')}
  if(missing(treat_var)) {stop('treatment variable not defined')}
  if(missing(treat_period)) {stop('intervention period not specified')}
  if(missing(time_var)) {stop('time_var not defined')}
  if(missing(unit_var)) {stop('unit_var not defined')}
  if(missing(weights)) {warning('ipw not specified. Equal weights assumed.')
    weights = NULL}
  if(missing(data)) {stop('data not found')}
  #### Stop if time var is not double 
  if(is.double(treat_period) == FALSE) {stop('treat_period is not double')}
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
#' @export
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
#' @export
psm_weights <- function(outcome,predictors,unit_var,data) {
### Get the model with highest accuracy 
model <- model_selection(outcome, predictors, unit_var, data)
### Optimal set of predictor variables 
predictor_star <- model[[1]]
### Build the data set
reg_data <- data_build(outcome, predictors, unit_var, data)
### Get the Data set with the chosen variables 
data_build_psm <- subset(reg_data, 
                         select = c("outcome", unit_var, predictor_star))
### Propensity Score Model Formual
reg_form <- as.formula(paste("outcome",
                             paste(predictor_star, collapse = " + "), 
                             sep = " ~ "))
### Estimate the Propensity Score Model
psm_model <- glm(reg_form, 
                 data = data_build_psm, 
                 family = binomial(link = "logit"))
# Generate propensity scores and IPWs
psm_weights <- data_build_psm %>% mutate(propensity = psm_model$fitted.values) %>% 
                        mutate(ipw = ifelse(outcome == 1, 1/propensity, 1/(1-propensity)))
### Compute the Prediction Error of the Model 
data_results <- data_build_psm %>% mutate(predicted = psm_model$fitted.values) %>% 
  mutate(outcome = as.double(outcome)) %>% 
  mutate(sq_error = (outcome - predicted)^2) 

### Accuracy Metric: Root Mean Squared Prediction Error 
rmspe <- data_results$sq_error %>% mean() %>% sqrt()
### Return Data Frame with Fips and IPW weights 
psm_results <- list(model, psm_weights, rmspe)
return(psm_results)
}

################################################################################
### Function 5: Model Selection for Propensity Score Matching  
################################################################################
#' @export
model_selection <- function(outcome,predictors,unit_var,data){
### Build data set with organized variables. This data set should only include the chosen predictors
### the predictors squared, interactions across predictors, and the outcome 
reg_data <- data_build(outcome,predictors,unit_var,data)
### Initial Set of Predictors
predictors_i <- colnames(reg_data)[-which(colnames(reg_data) == c("outcome") | colnames(reg_data) == c(unit_var))] 
pred_no <- length(predictors_i) 
predictor_star <- NULL
rmspe_master  <- matrix(data = NA, nrow = pred_no, ncol = 2)
### Counters 
rmspe_old = 100 
j = 1
### The Loop will run until we ran out of possible predictor variables 
### This is an forward stepwise selection algorithm. It begins with a simple univariate regression
### and iteratively augments the model by choosing the variable that improves more the RMSPE. 
### If there is no improvement, then break the loop and return the model. 
while(pred_no > 0) {
### Empty matrix store RMSPEs 
rmspe <- matrix(data = NA, nrow = pred_no, ncol = 3)
### Find the variable with highest stand alone prediction accuracy
for(i in 1:pred_no) {
### Model Definition 
data_subset <- subset(reg_data, select = c("outcome",
                                             predictor_star,
                                             predictors_i[i]))
### Estimate the logit model 
psm_model <- glm(outcome ~ ., data = data_subset, family = binomial(link = "logit"))
### Store Predicted Values and compute the Squared Prediction Error 
data_results <- data_subset %>% mutate(predicted = psm_model$fitted.values) %>% 
                               mutate(outcome = as.double(outcome)) %>% 
                               mutate(sq_error = (outcome - predicted)^2) 
### Accuracy Metric: Root Mean Squared Prediction Error 
rmspe_new <- data_results$sq_error %>% mean() %>% sqrt()
### Store the Results of the RMSPE for each estimated model 
rmspe[i,1] <- i
### Store the Accuracy of the model. 
rmspe[i,2] <- rmspe_new
### Compute the Accuracy Improvement relative to the previous lap of the loop 
rmspe[i,3] <- rmspe_old - rmspe_new
}

### Store Outer Loop Results
rmspe_master[j,1] <- j
rmspe_master[j,2] <- min(rmspe[,2])
### Choose the Model with highest accuracy improvement
model_index <- rmspe[which(rmspe[,3] == max(rmspe[,3])),1]
### Choose the predictor with highest accuracy
predictor_star <- c(predictor_star, predictors_i[model_index])
### Update the Number of Remaining Predictors 
predictors_i <- predictors_i[-model_index]
pred_no <- length(predictors_i)
### Break if adding a variable did not improved prediction accuracy 
delta_rmspe <- rmspe_old -  min(rmspe[,2])
if(delta_rmspe <= 0){
  break 
}
### Define rmspe old for the new loop 
rmspe_old <- min(rmspe[,2])
### Add one to the counter 
j = j + 1
}
results <- list(predictor_star,rmspe_master)
return(results)
}


################################################################################

### Function 6: Build Data set To Estimate Iteratively the Propensity Score Model 
################################################################################
#' @export
data_build <- function(outcome,predictors,unit_var,data){
  ### Stop if inputs are not defined 
  if(missing(outcome)) {stop('outcome variable not defined')}
  if(missing(predictors)) {stop('predictor variables not defined')}
  if(missing(unit_var)) {stop('unit_var not defined')}
  if(missing(data)) {stop('data not found')}
  ### Exlcude Ommited Ovservations 
  full_data <- data %>% na.omit()
  ### Rename outcome
  names(full_data)[which(colnames(full_data) == outcome)] <- "outcome"
  full_data$outcome <- as.factor(full_data$outcome)
  unit_var <- subset(full_data, select = unit_var)
  ### Outcome has the name of the outcome variable 
  y <- subset(full_data, select = outcome)
  ### Predictors is a vector with the names of all the variables that will be used as predictors 
  x <- subset(full_data, select = predictors)
  #### Squared Predictors 
  x2 <- x^2
  colnames(x2) <- paste0(colnames(x),'.sq')
  ### Predictors in Inverse Hyperbolic Sin Function  
  xasinh <- asinh(x)
  colnames(xasinh) <- paste0(colnames(x),'asinh')
  ### Interaction Terms and Intercept  
  interactions <- model.matrix( ~ .^2, data = x)
  inter_var <- interactions[,-which(colnames(interactions) == "(Intercept)")] %>% as.data.frame()
  ### Weird thing happening with the variable names when there is only one predictor 
  if(length(predictors) == 1){
    names(inter_var) <- predictors
  }
  ### Return Dataframe with ordered variables 
  regdata <- data.frame(unit_var, y, inter_var, x2, xasinh)
  return(regdata)
}
################################################################################



