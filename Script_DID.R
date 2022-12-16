################################################################################
################################################################################
### STAT-610 
### Final Project 
### Matched Difference-in-Difference 
### Author: Luis Navarro
### This Version: November 24, 2022
### Script: Implementation 
################################################################################
################################################################################
### Clean and load libraries 
rm(list = ls()) 
pacman::p_load(rmarkdown, formatR, kableExtra, reshape, gtable, PtProcess, dplyr, broom, readr, tidyverse, ggplot2, gridExtra, lattice, leaps, fixest, cli, testthat)

################################################################################
### Set the environment
directory <- "C:/Users/luise/OneDrive - Indiana University/Statistical Computing/Final_Project"
file_psm <- paste(directory, "/Data/county_data_project_psm.csv", sep = "")
file_did <- paste(directory, "/Data/county_data_project_did.csv", sep = "")
file_functions <- paste(directory, "/Code/Script_Functions.R", sep = "")
set.seed(1)
### Load the functions 
source(file_functions)
################################################################################
#### Parameters of the Simulation 
treatment_year = 2015 
cutoff_treatment_assignment <- rnorm(1, mean = 0.4, sd = 0.15)
sample_no <- 1000
ate_beta = 7
### Vector of Predictor Variables to choose from 
psm_dgp <- c("college","popestimate","female","taxes_pc", "exp_pc" ,"unemployment_rate","age4564","currexp")
################################################################################
### Data Cleaning 
### Load Data for Propensity Score Matching 
data_psm <- read.csv(file_psm, header = TRUE, sep = ",")
nfips <- length(data_psm$fips)
data_psm <- mutate(data_psm, female = popest_fem/popestimate,
                   some_college = some_college/100,
                   less_hs = less_hs/100,
                   hsd = hsd/100, 
                   college = college/100,
                   age513 = age513_tot/popestimate,
                   age1824 = age1824_tot/popestimate,
                   age2544 = age2544_tot/popestimate,
                   age4564 = age4564_tot/popestimate,
                   taxes_pc = 1000*total_taxes/popestimate,
                   exp_pc = 1000*total_expenditure/popestimate,
                   currexp = total_current_expend/total_expenditure,
                   popestimate = popestimate/1000,
                   unemployment_rate = unemployment_rate/100,
                   ### Non Random Treatment Assignment -- Complex Non-Linear Data Generating Process
                   ### Data Generating Process: function of college, female, unemployment_rate, logpop, age513, age4564, taxes_pc, currexp
                   treat_formula = unemployment_rate*female + sqrt(age4564) + 
                     sqrt(taxes_pc/exp_pc) + currexp^2 + 
                     log(1 + college) + log(popestimate) + 
                     rnorm(nfips, mean = 0, sd = 1),
                   treat_aux = percent_rank(treat_formula),
                   ### Deterministic Treatment
                   treat_det = case_when(treat_aux <= cutoff_treatment_assignment ~ 0,
                                         treat_aux > cutoff_treatment_assignment ~ 1),)

### Store Non Random Treatment Assignment 
treat_det <- data.frame(fips = data_psm$fips, treat_det = data_psm$treat_det)
##########
### Load Data for Difference in Difference Analysis 
data_did <- read.csv(file_did, header = TRUE, sep = ",")
### Create Fixed Effects 
fips_unique <- unique(data_did$fips)
fipsfe <- rnorm(length(fips_unique))
yr_unique <- unique(data_did$year)
yrfe <- rnorm(length(yr_unique))
### Map of Fixed Effects
fips_fe <- data.frame(fips = fips_unique, fips_fe = fipsfe)
year_fe <- data.frame(year = yr_unique, year_fe = yrfe)
### Merge Fixed Effects 
data_did <- merge(data_did, fips_fe, by = "fips")
data_did <- merge(data_did, year_fe, by = "year")
### Non-Random Treatment Assignment 
data_did <- merge(data_did,  data.frame(fips = data_psm$fips, treat_det = data_psm$treat_det), by = "fips")

### Create Dataset to Run the Regressions 
data_did <- mutate(data_did,  logproptax = log(1+property_tax),
                   prop_tax = property_tax/total_revenue,
                   post = case_when(data_did$year <= treatment_year ~ 0,
                                    data_did$year > treatment_year ~ 1),
                   did = treat_det*post,
                   epsilon = rnorm(length(data_did$fips)),
                   ### Create Fake Outcome To do the Simulation
                   fake_outcome = logproptax + ate_beta*(did) + treat_det + 
                     post + fips_fe + year_fe + epsilon)

################################################################################
### Simulation Study: Adding Predictors  
### Estimate the Un-weighted Difference in Difference Model
unweighted_model <- did_model(outcome = "fake_outcome", 
                              treat_var = "treat_det", treat_period = 2015,
                              time_var = "year", unit_var = "fips",
                              weights = NULL, data = data_did)
################################################################################
#### Baseline Scenario 
### Simulation Study: how the model performs when I add one variable to the set of predictors 
simulation_results <- list()
### Do the Simulation 
for(k in 1:length(psm_dgp)) {
  ###########################
  ### Compute the IPW using the Cross Validation Algorithm
  ipw <- psm_weights(outcome = "treat_det", 
                     predictors = psm_dgp[1:k],
                     unit_var = "fips", 
                     data = data_psm)
  ### IPW Weights 
  ipw_weights <- subset(ipw[[2]], select = c("fips", "propensity", "ipw"))
  
  ### Merge IPW in the DID data set to estimate the model with and without weights 
  data_did_model <- merge(data_did, ipw_weights, by = "fips")
  
  ### Estimate the Matched Difference-in-Difference Model. Statistical Inference is done 
  ### using the placebo distribution of the treatment variable under the specified model. 
  did_estimation <- did_randomization(outcome = "fake_outcome", 
                                      treat_var = "treat_det", treat_period = 2015,
                                      time_var = "year", unit_var = "fips",
                                      weights = "ipw", data = data_did_model, samples = sample_no)
  
  simulation_results[[k]] <- list(ipw, did_estimation)
}
##########################################################3

## Show the Results of the Simulation
ate <- matrix(data = NA, ncol = length(psm_dgp), nrow = 3)
for(k in 1:length(psm_dgp)) {
  ate[1,k] <- simulation_results[[k]][[2]][[1]][1,1]
  ate[2,k] <- simulation_results[[k]][[2]][[1]][1,2]
  
  ### Which one is the best model? Retrieve the RMSPE from each one 
  ate[3,k] <- simulation_results[[k]][[1]][[3]] 
}

### Diagnostics of the Simulation 
### Results from the Simulation
results_baseline <- t(ate) %>% as.data.frame() %>% na.omit()
names(results_baseline) <- c("ATE", "SE", "RMSPE")


### Which one is the best model? 
final_comp <- matrix(data = NA, ncol = 2, nrow = 2, 
                     dimnames = list(c("ATE","SE"), c("Unweighted","IPW")))
final_comp[1,1] <- unweighted_model$coefficients[[1]]
final_comp[2,1] <- unweighted_model$se[[1]]
final_comp[1,2] <- ate[1,ate[3,] == min(ate[3,])] 
final_comp[2,2] <- ate[2,ate[3,] == min(ate[3,])] 


#### Convergence of the Stepwise Regression 
combinations <- data_build(outcome = "treat_det", 
                           predictors = psm_dgp, 
                           unit_var = "fips", 
                           data = data_psm) %>% length() - 2
### All RMSPE results 
### k - model_selection (1) - rmspe_master (2)
rmspe_all <- list()


### Do the graphs the ATE hypothesis Test
emp_dist <- simulation_results[[7]][[2]][[2]]
ate_graph <- ate_inference_graph(emp_dist)

### Do the graphs of the convergence 
rmspe_results2 <- simulation_results[[1]][[1]][[1]][[2]] %>% as.data.frame() %>% na.omit()
convergence_plot2 <- ggplot(rmspe_results2, aes(x = V1, y = V2)) + theme_bw() + geom_line(color="darkblue") + labs(title = "RMSPE Convergence - 2 Predictors", y="RMSPE", x = "Number of Predictors")

rmspe_results4 <- simulation_results[[4]][[1]][[1]][[2]] %>% as.data.frame() %>% na.omit()
convergence_plot4 <- ggplot(rmspe_results4, aes(x = V1, y = V2)) + theme_bw() + geom_line(color="darkblue") + labs(title = "RMSPE Convergence - 4 Predictors", y="RMSPE", x = "Number of Predictors")

rmspe_results6 <- simulation_results[[6]][[1]][[1]][[2]] %>% as.data.frame() %>% na.omit()
convergence_plot6 <- ggplot(rmspe_results6, aes(x = V1, y = V2)) + theme_bw() + geom_line(color="darkblue") + labs(title = "RMSPE Convergence - 6 Predictors", y="RMSPE", x = "Number of Predictors")

rmspe_results8 <- simulation_results[[8]][[1]][[1]][[2]] %>% as.data.frame() %>% na.omit()
convergence_plot8 <- ggplot(rmspe_results8, aes(x = V1, y = V2)) + theme_bw() + geom_line(color="darkblue") + labs(title = "RMSPE Convergence - 8 Predictors", y="RMSPE", x = "Number of Predictors")

#### Simulation Test of the Precision of the Algorithm across different values of the ATE. 

#################
#### Simulation based tests: how well the algorithm works for different values of the ATE 
#### Simulation Parameters 
### Vector of Predictor Variables to choose from 

### Load Data for Difference in Difference Analysis 
data_did <- read.csv(file_did, header = TRUE, sep = ",")
### Create Fixed Effects 
fips_unique <- unique(data_did$fips)
fipsfe <- rnorm(length(fips_unique))
yr_unique <- unique(data_did$year)
yrfe <- rnorm(length(yr_unique))
### Map of Fixed Effects
fips_fe <- data.frame(fips = fips_unique, fips_fe = fipsfe)
year_fe <- data.frame(year = yr_unique, year_fe = yrfe)
### Merge Fixed Effects 
data_did <- merge(data_did, fips_fe, by = "fips")
data_did <- merge(data_did, year_fe, by = "year")
### Non-Random Treatment Assignment 
data_did <- merge(data_did,  data.frame(fips = data_psm$fips, treat_det = data_psm$treat_det), by = "fips")
data_did <- mutate(data_did,  logproptax = log(1+property_tax),
                   prop_tax = property_tax/total_revenue,
                   post = case_when(data_did$year <= treatment_year ~ 0,
                                    data_did$year > treatment_year ~ 1),
                   did = treat_det*post)
##########
#### Different Values for the ATE Beta 
ate_random <- runif(100, min = 1, max = 10) %>% round(digits = 2)
## Do the Simulation
test_sim_res <- list()
for(k in 1:length(ate_random)) {
  ### Create Fake Outcome To do the Simulation 
  ate_beta <- ate_random[k]
  
  data_did_sim <- data_did %>% mutate(epsilon = rnorm(length(data_did$fips)),
                                      fake_outcome = logproptax + ate_beta*(did) + treat_det + 
                                        post + fips_fe + year_fe + epsilon)
  
  ### Unweighted Model 
  unweighted_model <- did_model(outcome = "fake_outcome", 
                                treat_var = "treat_det", treat_period = 2015,
                                time_var = "year", unit_var = "fips",
                                weights = NULL, data = data_did_sim)
  ### Unweighted Beta 
  unw_beta <- coef(unweighted_model)
  ### Compute the IPWs    
  ipw <- psm_weights(outcome = "treat_det", 
                     predictors = psm_dgp[1:4],
                     unit_var = "fips", 
                     data = data_psm)
  ### IPW Weights 
  ipw_weights <- subset(ipw[[2]], select = c("fips", "propensity", "ipw"))
  ### Merge IPW in the DID data set to estimate the model with and without weights 
  data_did_model <- merge(data_did_sim, ipw_weights, by = "fips")
  ### Estimate the Matched Difference-in-Difference Model. Statistical Inference is done 
  ### using the placebo distribution of the treatment variable under the specified model. 
  did_estimation <- did_randomization(outcome = "fake_outcome", 
                                      treat_var = "treat_det", treat_period = 2015,
                                      time_var = "year", unit_var = "fips",
                                      weights = "ipw", data = data_did_model, samples = 1)
  ### IPW BETA 
  weight_beta <- did_estimation[[1]][[1]]
  
  prediction_error1 <- (unw_beta - weight_beta)^2
  prediction_error2 <- (ate_beta - weight_beta)^2
  test_sim_res[[k]] <- list(ipw, did_estimation, prediction_error1, prediction_error2)
}
##########################################################3

## Show the Results of the Simulation
ate_error <- matrix(data = NA, nrow = length(ate_random), ncol = 2)
for(k in 1:length(ate_random)) {
  ate_error[k,1] <- test_sim_res[[k]][[3]][[1]]
  ate_error[k,2] <- test_sim_res[[k]][[4]][[1]]
}


#### Graphs Showing the Distribution of the Prediction Error 
emp_dist_error1 <- ate_error[,1] %>% as.data.frame()
names(emp_dist_error1) <- c("V1")
rmspe1_dist <- ggplot(emp_dist_error1, aes(x = V1)) + 
  geom_density(alpha=0.25) + 
  theme_bw() + geom_density(color="darkblue") +
  labs(title = "Prediction Error 1: Unweighted vs Weighted", y="Density", x = "Squared Error")


emp_dist_error2 <- ate_error[,2] %>% as.data.frame()
names(emp_dist_error2) <- c("V1")
rmspe2_dist <- ggplot(emp_dist_error2, aes(x = V1)) + 
  geom_density(alpha=0.25) + 
  theme_bw() + geom_density(color="darkblue") +
  labs(title = "Prediction Error 2: True ATE vs Weighted Estimate", y="Density", x = "Squared Error")




#### Comparing the Results 
rmspe1 <- ate_error[,1] %>% mean() %>% sqrt()
rmspe2 <- ate_error[,2] %>% mean() %>% sqrt()

pred_errors_results <- data.frame(RMSPE1 = rmspe1, RMSPE2 = rmspe2)