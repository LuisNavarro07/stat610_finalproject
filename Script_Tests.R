################################################################################
################################################################################
### STAT-610 
### Final Project 
### Matched Difference-in-Difference 
### Author: Luis Navarro
### This Version: December 15, 2022
### Script: Testing the Functions 
################################################################################
################################################################################
### Clean and load libraries 
rm(list = ls()) 
pacman::p_load(dplyr, broom, readr, tidyverse, ggplot2, lattice, leaps, fixest, cli, testthat, devtools, roxygen2)
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
ate_beta = 7 
cutoff_treatment_assignment <- rnorm(1, mean = 0.4, sd = 0.15)
sample_no <- 100
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
                   ### Fake Unitary Weights 
                   test_weights = 1, 
                   ### Create Fake Outcome To do the Simulation
                   fake_outcome = logproptax + ate_beta*(did) + treat_det + post + fips_fe + year_fe + epsilon)

################################################################################
################################################################################
### Test function: Did Randomization
### Warning if missing samples  
test_that("Show warning of missing samples", {
  expect_warning(did_randomization(outcome = "fake_outcome",treat_var = "treat_det", treat_period = 2015, 
                         time_var = "year", unit_var  = "fips", weights = NULL ,data = data_did), 'number of samples not defined. Default is 100')
})

### Informal Testing for number of samples 
random_without_samples <- did_randomization(outcome = "fake_outcome",treat_var = "treat_det", treat_period = 2015, 
                                            time_var = "year", unit_var  = "fips", weights = NULL ,data = data_did)[[1]]


random_with_samples <- did_randomization(outcome = "fake_outcome",treat_var = "treat_det", treat_period = 2015, samples = 100, 
                                            time_var = "year", unit_var  = "fips", weights = NULL ,data = data_did)[[1]]

test_that("Verify that the function leads to the same statistical inference", {
  expect_identical(random_without_samples, random_with_samples)
})

#### Check error messages appear even if they are written at the did model function 

### Error if the treatment period is not a number 
test_that("Show errors that stop the function", {
  ### Error if the treatment period is not a number 
  expect_error(did_randomization(outcome = "fake_outcome",treat_var = "treat_det", treat_period = "2015", 
                                 time_var = "year", unit_var  = "fips", weights = NULL ,data = data_did, samples = 100), 'treat_period is not double')
  ### Errors  
  expect_error(did_randomization(treat_var = "treat_det", treat_period = 2015, time_var = "year", 
                                 unit_var  = "fips", weights = NULL ,data = data_did, samples = 100), 'outcome variable not defined')
  
  expect_error(did_randomization(outcome = "fake_outcome", treat_period = 2015, time_var = "year", 
                                 unit_var  = "fips", weights = NULL ,data = data_did, samples = 100), 'treatment variable not defined')
  
  expect_error(did_randomization(outcome = "fake_outcome", treat_var = "treat_det", time_var = "year", 
                                 unit_var  = "fips", weights = NULL ,data = data_did, samples = 100), 'intervention period not specified')
  
  expect_error(did_randomization(outcome = "fake_outcome", treat_var = "treat_det", treat_period = 2015, 
                                 unit_var  = "fips", weights = NULL ,data = data_did, samples = 100), 'time_var not defined')
  
  expect_error(did_randomization(outcome = "fake_outcome", treat_var = "treat_det", treat_period = 2015, 
                                 time_var = "year" , weights = NULL ,data = data_did, samples = 100), 'unit_var not defined')
  
  expect_error(did_randomization(outcome = "fake_outcome", treat_var = "treat_det", treat_period = 2015, 
                         unit_var  = "fips", time_var = "year" , weights = NULL, samples = 100), 'data not found')
})


################################################################################
################################################################################
### Test Function: Did model  

### Same Coefficient Test
### This is a rather informal test 
did_without_weights <- did_model(outcome = "fake_outcome",treat_var = "treat_det", 
                                 treat_period = 2015, time_var = "year", 
                                 unit_var  = "fips", data = data_did) %>% coef()

did_with_weights <- did_model(outcome = "fake_outcome",treat_var = "treat_det", 
                              treat_period = 2015, time_var = "year", weights = NULL, 
                              unit_var  = "fips", data = data_did) %>% coef()


#### Shows the warning when weights are not specified
test_that("Shows warnings when weights are not specified", {
  ### When weights are not specified 
  expect_warning(did_model(outcome = "fake_outcome",treat_var = "treat_det", treat_period = 2015, 
                         time_var = "year", unit_var  = "fips", data = data_did), 'ipw not specified. Equal weights assumed')
})

### Tests the error messages appear 
test_that("Show errors that stop the function", {
### Error if the treatment period is not a number 
expect_error(did_model(outcome = "fake_outcome",treat_var = "treat_det", treat_period = "2015", 
                         time_var = "year", unit_var  = "fips", weights = NULL ,data = data_did), 'treat_period is not double')
### Errors if missing arguments 
expect_error(did_model(treat_var = "treat_det", treat_period = 2015, time_var = "year", 
                       unit_var  = "fips", weights = NULL ,data = data_did), 'outcome variable not defined')

expect_error(did_model(outcome = "fake_outcome", treat_period = 2015, time_var = "year", 
                       unit_var  = "fips", weights = NULL ,data = data_did), 'treatment variable not defined')

expect_error(did_model(outcome = "fake_outcome", treat_var = "treat_det", time_var = "year", 
                       unit_var  = "fips", weights = NULL ,data = data_did), 'intervention period not specified')

expect_error(did_model(outcome = "fake_outcome", treat_var = "treat_det", treat_period = 2015, 
                       unit_var  = "fips", weights = NULL ,data = data_did), 'time_var not defined')

expect_error(did_model(outcome = "fake_outcome", treat_var = "treat_det", treat_period = 2015, 
                       time_var = "year" , weights = NULL ,data = data_did), 'unit_var not defined')

expect_error(did_model(outcome = "fake_outcome", treat_var = "treat_det", treat_period = 2015, 
                       unit_var  = "fips", time_var = "year" , weights = NULL), 'data not found')
})

#### Estimates the same model when weights are set to null or not set at all 
test_that("coefficient estimates are the same", {
  expect_identical(did_with_weights, did_without_weights)
})


################################################################################
################################################################################
### Test Function: Data Build 
pred_test  <- c("college","popestimate","female","taxes_pc", "exp_pc" ,"unemployment_rate","age4564","currexp")

### Missing Input Tests 
test_that("Shows missing input eror", {
  expect_error(data_build(predictors = pred_test, unit_var = "fips", data = data_psm), 'outcome variable not defined')  
  expect_error(data_build(outcome = "treat_det",  unit_var = "fips", data = data_psm), 'predictor variables not defined')
  expect_error(data_build(outcome = "treat_det", predictors = pred_test, data = data_psm), 'unit_var not defined')
  expect_error(data_build(outcome = "treat_det", predictors = pred_test, unit_var = "fips"), 'data not found')
})


### Test that it is delimiting the model space correctly. The data frame should have 2 + 3n + 2^n
dim_fun <- function(vector){
  n <- length(vector)
  if (n == 1) {
  dimension <- 2 + 3   
  } else {
  dimension <- 2 + 3*n + (factorial(n)/(2*factorial(n - 2)))
  }
  return(dimension)
}

### Test the Dimensions of the model space are well defined 
test_that("Dimensions of the model space are not correct", {
  expect_equal(length(data_build(outcome = "treat_det", predictors = pred_test[1:1], unit_var = "fips", data = data_psm)), dim_fun(pred_test[1:1]))
  expect_equal(length(data_build(outcome = "treat_det", predictors = pred_test[1:2], unit_var = "fips", data = data_psm)), dim_fun(pred_test[1:2]))
  expect_equal(length(data_build(outcome = "treat_det", predictors = pred_test[1:3], unit_var = "fips", data = data_psm)), dim_fun(pred_test[1:3]))
  expect_equal(length(data_build(outcome = "treat_det", predictors = pred_test[1:4], unit_var = "fips", data = data_psm)), dim_fun(pred_test[1:4]))
  expect_equal(length(data_build(outcome = "treat_det", predictors = pred_test[1:5], unit_var = "fips", data = data_psm)), dim_fun(pred_test[1:5]))
  expect_equal(length(data_build(outcome = "treat_det", predictors = pred_test[1:6], unit_var = "fips", data = data_psm)), dim_fun(pred_test[1:6]))
  expect_equal(length(data_build(outcome = "treat_det", predictors = pred_test[1:7], unit_var = "fips", data = data_psm)), dim_fun(pred_test[1:7]))
  expect_equal(length(data_build(outcome = "treat_det", predictors = pred_test[1:8], unit_var = "fips", data = data_psm)), dim_fun(pred_test[1:8]))
})
