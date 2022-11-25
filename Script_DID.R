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
ate = 2 
### Variables to use in the simulation
covariates <- c("fips", "less_hs","hsd","some_college","college",
                "logpop","unemployment_rate","female",
                "age513","age1824","age2544","age4564",
                "health","justice","human_serv","community",
                "total_revenue", "taxes_pc","exp_pc","currexp")
################################################################################
### Load Data for Propensity Score Matching 
data_psm <- read.csv(file_psm, header = TRUE, sep = ",")
data_psm <- data_psm %>% mutate(data_psm, female = popest_fem/popestimate,
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
                                          treat_aux = percent_rank((sqrt(female*college + logpop^2 + unemployment_rate) + 
                                                                     sqrt(cos(some_college)^2 + sin(age2544)^2)) + 
                                                                     abs(cos(less_hs)) + 
                                                                     log(1 + taxes_pc + currexp)),
                                          ### Deterministic Treatment
                                          treat_det = round(treat_aux,0))

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
                              fake_outcome = logproptax + ate*(did) + treat_det + post + fips_fe + year_fe + epsilon)

################################################################################
### Compute the Inverse Probabilty Weights 
cv_options <- trainControl(method = "cv", number = 10, savePredictions = TRUE)
ipw <- psm_weights(outcome = "treat_det", 
                   predictors = c("female","logpop","unemployment_rate","college","age2544","taxes_pc"),
                   data = data_psm)
### IPW Weights 
ipw_weights <- ipw[[2]]
### Merge IPW in the DID dataset to estimate the model with and without weights 
data_did <- merge(data_did, ipw_weights, by = "fips")

###############################################################################
### Simulation 
did_estimation <- did_randomization(outcome = "fake_outcome", 
              treat_var = "treat_det", treat_period = 2015,
              time_var = "year", unit_var = "fips",
              weights = "ipw", data = data_did, samples = 1000)


### Unweighted
unweighted_model <- did_model(outcome = "fake_outcome", 
                              treat_var = "treat_det", treat_period = 2015,
                              time_var = "year", unit_var = "fips",
                              weights = NULL, data = data_did)

#### Compare the Results 
compare_results <- matrix(nrow = 3, ncol =3, dimnames = list(c("ATE","SE","pval"),c("Real","Unweighted","IPW")))
### Real 
compare_results[1,1] <- ate
compare_results[2,1] <- 0
compare_results[3,1] <- 0

### Store unweighted Results 
compare_results[1,2] <- unweighted_model$coeftable$Estimate
compare_results[2,2] <- unweighted_model$coeftable$`Std. Error`
compare_results[3,2] <- unweighted_model$coeftable$`Pr(>|t|)`

### Weighted Regression 
compare_results[1,3] <- did_estimation$coefficient
compare_results[2,3] <- did_estimation$std_error
compare_results[3,3] <- did_estimation$pval

compare_results

