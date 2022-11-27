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
library(rmarkdown)
library(knitr)
library(formatR)
knitr::opts_chunk$set(echo = TRUE)
library(kableExtra)
library(dplyr)
library(jtools)
library(fixest)
library(broom)
library(readr)
library(tidyverse)
library(ggplot2)
library(lattice)
library(caret)
library(leaps)
################################################################################
### Set the environment
directory <- "C:/Users/luise/OneDrive - Indiana University/Statistical Computing/Final_Project"
file_psm <- paste(directory, "/Data/county_data_project_psm.csv", sep = "")
file_did <- paste(directory, "/Data/county_data_project_did.csv", sep = "")
file_functions <- paste(directory, "/Code/Script_Functions.R", sep = "")
set.seed(1234)
### Load the functions 
source(file_functions)

################################################################################
#### Parameters of the Simulation 
treatment_year = 2015 
ate_beta = 2 
cutoff_treatment_assignment <- rnorm(1, mean = 0.5, sd = 0.1)
### Variables to use in the simulation
covariates <- c("fips", "less_hs","hsd","some_college","college",
                "logpop","unemployment_rate","female",
                "age513","age1824","age2544","age4564",
                "health","justice","human_serv","community",
                "total_revenue", "taxes_pc","exp_pc","currexp")
### Cross Validation Parameters
cv_options <- trainControl(method = "cv", number = 10, savePredictions = TRUE)
################################################################################
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
                             logpop = log(popestimate/1000),
                             unemployment_rate = unemployment_rate/100,
                             ### Non Random Treatment Assignment -- Complex Non-Linear Data Generating Process
                             ### Data Generating Process: function of college, female, unemployment_rate, logpop, age513, age4564, taxes_pc, currexp
                             treat_formula = sqrt(sin(unemployment_rate + female + age4564)) + 
                                             sqrt(taxes_pc/exp_pc) + cos(currexp)^2 + 
                                             rnorm(nfips, mean = college, sd = college/2) + 
                                             sqrt(rnorm(nfips, mean = abs(logpop), sd = abs(logpop)/4)),
                             treat_aux = percent_rank(treat_formula),
                             ### Deterministic Treatment
                             treat_det = case_when(treat_aux <= cutoff_treatment_assignment/2 ~ 0,
                                                   treat_aux > cutoff_treatment_assignment/2 & treat_aux <= 1.5*cutoff_treatment_assignment ~ 1,
                                                   treat_aux > 1.5*cutoff_treatment_assignment ~ 0),)

### Balance between treatment and control 
data_psm$treat_det %>% table()
### Simplify Data set 
data_psm <- subset(data_psm, select = c("treat_det",covariates))
### Store Non Random Treatment Assignment 
treat_det <- data.frame(fips = data_psm$fips, treat_det = data_psm$treat_det)

###############################################################################
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
                              fake_outcome = logproptax + ate_beta*(did) + treat_det + post + fips_fe + year_fe + epsilon)

################################################################################
### Simulation Study: Adding Predictors  

### Estimate the Un-weighted Difference in Difference Model
unweighted_model <- did_model(outcome = "fake_outcome", 
                              treat_var = "treat_det", treat_period = 2015,
                              time_var = "year", unit_var = "fips",
                              weights = NULL, data = data_did)
################################################################################
### Simulation Study: how the model performs when I add one variable to the set of predictors 
psm_dgp <- c("college","logpop","female","taxes_pc", "unemployment_rate","age513","currexp")
simulation_results <- list()
### Do the Simulation 
for(k in 1:length(psm_dgp)) {
### Compute the IPW using the Cross Validation Algorithm
ipw <- psm_weights(outcome = "treat_det", 
                   predictors = psm_dgp[1:k],
                   id_var = "fips", 
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
              weights = "ipw", data = data_did_model, samples = 1000)

simulation_results[[k]] <- list(ipw, did_estimation)
}
##########################################################3


### Show the Results of the Simulation
ate <- matrix(data = NA, nrow = length(psm_dgp), ncol = 3)
for(k in 1:length(psm_dgp)) {
ate[k,1] <- simulation_results[[k]][[2]][1,1]
ate[k,2] <- simulation_results[[k]][[2]][1,2]
ate[k,3] <- k
}

unweighted <- unweighted_model$coefficients[[1]]

### Result from the Unweighted Regression 
unweighted
### Results from the Simulation
ate

ate_comp <- ate %>% as.data.frame() %>% mutate( ,unweighted = unweighted,  
                                                 real_ate = ate_beta)

sim_results_plot <- ggplot(ate_comp, aes(x=V3)) +
  geom_line(aes(y = V1), color = "black") +
  geom_line(aes(y = unweighted), color = "red") +
  geom_line(aes(y = real_ate), color = "blue") +
  theme_void() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size=10, angle=0), 
        axis.text.y = element_text(size=10, angle=0)) + 
  labs(title = "Simulation Results", y="Average Treatment Effect", x = "Number of Predictors")
print(sim_results_plot)



#### Compare the Results 
matrix_coef <- matrix(nrow = 3, ncol = length(psm_dgp) + 1, dimnames = list(c("ATE","SE","pval"), 
                                                                            c("Unweighted",seq(1:7))))

matrix_coef[1,1] <- unweighted_model$coeftable[1] %>% as.double()
matrix_coef[2,1] <- unweighted_model$coeftable[2] %>% as.double()
matrix_coef[3,1] <- unweighted_model$coeftable[4] %>% as.double()

for(k in 1:length(psm_dgp)) {
## Weighted Regression: beta, se, pval
  matrix_coef[1,k + 1] <- simulation_results[[k]][[2]][1]  %>% as.double()
  matrix_coef[2,k + 1] <- simulation_results[[k]][[2]][2]  %>% as.double()
  matrix_coef[3,k + 1] <- simulation_results[[k]][[2]][3]  %>% as.double()
} 

matrix_coef %>% print()

### As the number of predictors increase, the accuracy of the estimate increases as well. 
