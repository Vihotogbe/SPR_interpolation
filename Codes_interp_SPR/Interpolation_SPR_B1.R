########################### SPR interpolation function ############################

## Loading the necessary packages
library(bigstatsr)
library(bigmemory)

# SPR interpolation method
spr_interp <- function(data_rcm_grid, n_spatterns = 10, data_obs_grid) {
  ## data_rcm_grid is the data from which the spatial patterns are calculated : n x p
  ## n is the number of days (rows) and p is the number of stations (columns)
  ## n_spatterns is the number of spatial patterns with p rows (stations)
  ## data_obs_grid is the data where interpolation is needed : m x p
  ## m is the number of days (rows) and p is the number of stations (columns)
  ## seed is the seed for the random number generator

  # Calculating the spatial patterns from RCM DATA
  data_rcm_FBM <- as_FBM(data_rcm_grid)
  pca_rcm <- big_randomSVD(data_rcm_FBM, k = n_spatterns, # ncores = 15,
                 fun.scaling = big_scale(center = TRUE, scale = FALSE))
  ## ncore = NULL: use all available cores in your system
  spatterns <- pca_rcm$v # spatial patterns : p x n_spatterns

  # Getting the indexes of the columns (stations) with NA values
  indexes_NA_col <- which(is.na(data_obs_grid[1, ]))
  # Removing the columns (stations) with NA values
  data_obs_valid <- data_obs_grid[, - indexes_NA_col]
  data_obs_valid <- as.matrix(data_obs_valid)
  nrow_obs <- nrow(data_obs_grid)
  ncol_obs <- ncol(data_obs_grid)
  ## spatterns : is restricted to the valid columns (stations) to match
  ## the data_test_valid
  spatterns_model <- spatterns[- indexes_NA_col, ]
  centres_pca <- pca_rcm$center
  # Results of this function : elements of the list to return by the function
  ## Add annytihng to be returned here
  data_interp_spr <- matrix(NA, nrow = nrow_obs, ncol = ncol_obs)

  # Function of PSPR interpolation
  for (i in 1:nrow_obs){
    # Fitting the linear model for each day
    data_Y <- data_obs_valid[i, ]
    data_Y_centered <- data_Y - centres_pca[- indexes_NA_col] #mean(data_Y) #- climatology_model[indice_mois, ]
    model <- lm(data_Y_centered ~ spatterns_model - 1)
    coefs_val <- as.matrix(coef(model))

    # Adding variability
    # set.seed(seed) # Remove the seed after to make it completely random
    # variability <- sample(residuals(model), nrow(spatterns), replace = T)
    # In case we use this, add the seed as an argument to the function

    # Generating the interpolated data : the gridded meteorological data
    ## Generated data for the complete RCM grid
    data_obs_interp <- spatterns %*% coefs_val + centres_pca # + variability
    data_obs_interp <- as.vector(data_obs_interp)

    # Tables to return
    data_interp_spr[i, ] <- data_obs_interp
  }
  return(data_interp_spr) # Return the interpolated : m x p
}
################################# END OF THE FUNCTION #################################


########################## SIMULATIONS ###########################################

# LOADING DATA
load(file = "Data/data.coord.sf.RData") # rotated coord
load(file = "Data/data.pr.RData") # precipitation
load(file = "Data/data.tasmin.RData") # tasmin
load(file = "Data/data.tasmax.RData") # tasmax
# load(file = "Data/climatology_pr.RData") # climatology precipitation
# load(file = "Data/climatology_tasmin.RData") # climatology tasmin
# load(file = "Data/climatology_tasmax.RData") # climatology tasmax


# REGIONS A and B DETERMINATION : Showing the 3 sizes of region
rlon <- data.coord.sf[, "lon"]
rlat <- data.coord.sf[, "lat"]
region_A1 <- which(rlat >= -2 & rlat <= 4 & rlon >= 13 & rlon <= 19, arr.ind = TRUE)
region_A2 <- which(rlat >= -1 & rlat <= 3 & rlon >= 14 & rlon <= 18, arr.ind = TRUE)
region_A3 <- which(rlat >= 0 & rlat <= 2 & rlon >= 15 & rlon <= 17, arr.ind = TRUE)
region_B1 <- which(rlat >= 7 & rlat <= 13 & rlon >= 14 & rlon <= 20, arr.ind = TRUE)
region_B2 <- which(rlat >= 8 & rlat <= 12 & rlon >= 15 & rlon <= 19, arr.ind = TRUE)
region_B3 <- which(rlat >= 9 & rlat <= 11 & rlon >= 16 & rlon <= 18, arr.ind = TRUE)


# Defining the period of RCM and observations data
rcm_period <- 1:10950 # 30 years of data (1980-2009)
obs_period_1 <- 7301:10950 # 10 years of data (2000-2009)
obs_period_2 <- 10951:14600 # 10 years of data (2010-2019)
obs_period_3 <- 14601:18250 # 10 years of data (2020-2029)


###############################################################################
########## HYPERPARAMETERS SELECTION FOR THE METHODS ##########################
########## Select the best hyperparameters for each method ####################
########## Do this for each variable, region, and period ######################
########## The % of missing data are 10% - 30% - 50% - 70% - 90% ##############
###############################################################################
# 1. Choose the region and the period
# 2. In each region, do for the 5 % of NA, 3 variables concerned
# 3. Save the best hyperparameters for each variable in OVERLEAF
# 3. Return to 1. and change the region and/or the period
###############################################################################


######################################################################
########## The following part is to be run for each region ###########
####### It does for all the 5 densities before change the region #####
######################################################################

## Choose the period and the region
obs_period <- obs_period_1 # choose the period
region <- region_B1 # choose the region

seed <- 2244677 # Define the seed

##############################################################
## Load packages for parralel computing
library(doParallel)
library(foreach)
cl <- makeCluster(15) # 15 cores
registerDoParallel(cl)
###############################################################

#######################################################################################
################################## 10% of missing data ################################

### RCM data
data_pr_grid <- data.pr[rcm_period, region]
data_tasmin_grid <- data.tasmin[rcm_period, region]
data_tasmax_grid <- data.tasmax[rcm_period, region]
# climatology.pr <- climatology_pr[, region]
# climatology.tasmin <- climatology_tasmin[, region]
# climatology.tasmax <- climatology_tasmax[, region]


### Observations data
data_pr_obs <- data.pr[obs_period, region]
### Different transformations for PRECIPITATION DATA
data_pr_obs[data_pr_obs <= 10^(-5)] <- 10^(-5) # lower bundary for positivity transformation
data_pr_obs <- log(exp(data_pr_obs) - 1) # ensuring positive values for precipitations
data_tasmin_obs <- data.tasmin[obs_period, region]
data_tasmax_obs <- data.tasmax[obs_period, region]


# # Extend the climatology over 10 years
# climatology_fun <- function(climatology, nyears) {
#   days_by_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
#   month_indexe <- rep(1:12, times = days_by_month)
#   climatology_year <- climatology[month_indexe, ]
#   climatology_nyears <- matrix(rep(climatology_year, times = nyears), nrow = 365 * nyears, byrow = TRUE)
#   return(climatology_nyears)
# }

# # Calculate anomalies for grid data
# anomalies_pr_grid <- data_pr_grid - climatology_fun(climatology.pr, nyears = 30)
# anomalies_tasmin_grid <- data_tasmin_grid - climatology_fun(climatology.tasmin, nyears = 30)
# anomalies_tasmax_grid <- data_tasmax_grid - climatology_fun(climatology.tasmax, nyears = 30)
  
# # Calculate anomalies for observations
# anomalies_pr_obs <- data_pr_obs - climatology_fun(climatology.pr, nyears = 10)
# anomalies_tasmin_obs <- data_tasmin_obs - climatology_fun(climatology.tasmin, nyears = 10)
# anomalies_tasmax_obs <- data_tasmax_obs - climatology_fun(climatology.tasmax, nyears = 10)

### Train and Test set
perc <- 0.1 # Change this percentage to represent each density
custom_round <- function(x) { # Custom rounding function
  ifelse(x - floor(x) == 0.5, ceiling(x), round(x))
}
n_test <- custom_round(perc * ncol(data_pr_obs)) # number of cells in test set
set.seed(seed) # set the seed
test_indexes <- sample(1:ncol(data_pr_obs), size = n_test, replace = FALSE)

### Extract data for test and train sets
#### data_..._obs contains the observed data with test indexes set to NA
data_pr_obs[, test_indexes] <- NA # set the test columns indexes to NA
data_tasmin_obs[, test_indexes] <- NA # set the test columns indexes to NA
data_tasmax_obs[, test_indexes] <- NA # set the test columns indexes to NA
# anomalies_pr_obs[, test_indexes] <- NA # set the test columns indexes to NA
# anomalies_tasmin_obs[, test_indexes] <- NA # set the test columns indexes to NA
# anomalies_tasmax_obs[, test_indexes] <- NA # set the test columns indexes to NA


### Now, apply SPR to data_..._val with different values of n_spatterns and
### return the number with the lower RMSE


## Precipitation
data_pr_test <- data_pr_obs[, - test_indexes]
n_pr <- min(ncol(data_pr_grid), ncol(data_pr_test))
prop_NA_pr <- 0.50
n_spatterns_pr <- custom_round(n_pr * prop_NA_pr)
spr_interp_pr_B1_10 <- foreach(i = 1:length(n_spatterns_pr), .combine = rbind,
                            .packages = c("bigmemory", "bigstatsr")) %dopar% {
  # Perform SPR for each n_spatterns
  number <- n_spatterns_pr
  interp_spr <- spr_interp(data_rcm_grid = data_pr_grid, n_spatterns = number,
                         data_obs_grid = data_pr_obs)
  # Back transformation for positivity
  interp_spr <- log(1 + exp(interp_spr)) # back transformation for positivity
  interp_spr[interp_spr < 0.5] <- 0 # lower bundary to remove noise
  return(interp_spr)
}

# Save the results in appropriate folder
save(spr_interp_pr_B1_10, file = "Projet_1/Interpolation_finale_3/Simulations_SPR/Region_B1/spr_interp_pr_B1_10.RData")


## Minimum temperature
data_tasmin_test <- data_tasmin_obs[, - test_indexes]
n_tasmin <- min(ncol(data_tasmin_grid), ncol(data_tasmin_test))
prop_NA_tasmin <- 0.30
n_spatterns_tasmin <- custom_round(n_tasmin * prop_NA_tasmin)
spr_interp_tasmin_B1_10 <- foreach(i = 1:length(n_spatterns_tasmin), .combine = rbind,
                            .packages = c("bigmemory", "bigstatsr")) %dopar% {
  # Perform SPR for each n_spatterns
  number <- n_spatterns_tasmin
  interp_spr <- spr_interp(data_rcm_grid = data_tasmin_grid, n_spatterns = number,
                         data_obs_grid = data_tasmin_obs)
  return(interp_spr)
}

# Save the results in appropriate folder
save(spr_interp_tasmin_B1_10, file = "Projet_1/Interpolation_finale_3/Simulations_SPR/Region_B1/spr_interp_tasmin_B1_10.RData")


## Maximum temperature
data_tasmax_test <- data_tasmax_obs[, - test_indexes]
n_tasmax <- min(ncol(data_tasmax_grid), ncol(data_tasmax_test))
prop_NA_tasmax <- 0.40
n_spatterns_tasmax <- custom_round(n_tasmax * prop_NA_tasmax)
spr_interp_tasmax_B1_10 <- foreach(i = 1:length(n_spatterns_tasmax), .combine = rbind,
                            .packages = c("bigmemory", "bigstatsr")) %dopar% {
  # Perform SPR for each n_spatterns
  number <- n_spatterns_tasmax
  interp_spr <- spr_interp(data_rcm_grid = data_tasmax_grid, n_spatterns = number,
                         data_obs_grid = data_tasmax_obs)
  return(interp_spr)
}

# Save the results in appropriate folder
save(spr_interp_tasmax_B1_10, file = "Projet_1/Interpolation_finale_3/Simulations_SPR/Region_B1/spr_interp_tasmax_B1_10.RData")

#######################################################################################
################################## 30% of missing data ################################


### RCM data
data_pr_grid <- data.pr[rcm_period, region]
data_tasmin_grid <- data.tasmin[rcm_period, region]
data_tasmax_grid <- data.tasmax[rcm_period, region]
# climatology.pr <- climatology_pr[, region]
# climatology.tasmin <- climatology_tasmin[, region]
# climatology.tasmax <- climatology_tasmax[, region]


### Observations data
data_pr_obs <- data.pr[obs_period, region]
### Different transformations for PRECIPITATION DATA
data_pr_obs[data_pr_obs <= 10^(-5)] <- 10^(-5) # lower bundary for positivity transformation
data_pr_obs <- log(exp(data_pr_obs) - 1) # ensuring positive values for precipitations
data_tasmin_obs <- data.tasmin[obs_period, region]
data_tasmax_obs <- data.tasmax[obs_period, region]


# # Extend the climatology over 10 years
# climatology_fun <- function(climatology, nyears) {
#   days_by_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
#   month_indexe <- rep(1:12, times = days_by_month)
#   climatology_year <- climatology[month_indexe, ]
#   climatology_nyears <- matrix(rep(climatology_year, times = nyears), nrow = 365 * nyears, byrow = TRUE)
#   return(climatology_nyears)
# }

# # Calculate anomalies for grid data
# anomalies_pr_grid <- data_pr_grid - climatology_fun(climatology.pr, nyears = 30)
# anomalies_tasmin_grid <- data_tasmin_grid - climatology_fun(climatology.tasmin, nyears = 30)
# anomalies_tasmax_grid <- data_tasmax_grid - climatology_fun(climatology.tasmax, nyears = 30)

# # Calculate anomalies for observations
# anomalies_pr_obs <- data_pr_obs - climatology_fun(climatology.pr, nyears = 10)
# anomalies_tasmin_obs <- data_tasmin_obs - climatology_fun(climatology.tasmin, nyears = 10)
# anomalies_tasmax_obs <- data_tasmax_obs - climatology_fun(climatology.tasmax, nyears = 10)

### Train and Test set
perc <- 0.3 # Change this percentage to represent each density
custom_round <- function(x) { # Custom rounding function
  ifelse(x - floor(x) == 0.5, ceiling(x), round(x))
}
n_test <- custom_round(perc * ncol(data_pr_obs)) # number of cells in test set
set.seed(seed) # set the seed
test_indexes <- sample(1:ncol(data_pr_obs), size = n_test, replace = FALSE)

### Extract data for test and train sets
#### data_..._obs contains the observed data with test indexes set to NA
data_pr_obs[, test_indexes] <- NA # set the test columns indexes to NA
data_tasmin_obs[, test_indexes] <- NA # set the test columns indexes to NA
data_tasmax_obs[, test_indexes] <- NA # set the test columns indexes to NA
# anomalies_pr_obs[, test_indexes] <- NA # set the test columns indexes to NA
# anomalies_tasmin_obs[, test_indexes] <- NA # set the test columns indexes to NA
# anomalies_tasmax_obs[, test_indexes] <- NA # set the test columns indexes to NA


### Now, apply SPR to data_..._val with different values of n_spatterns and
### return the number with the lower RMSE


## Precipitation
data_pr_test <- data_pr_obs[, - test_indexes]
n_pr <- min(ncol(data_pr_grid), ncol(data_pr_test))
prop_NA_pr <- 0.50
n_spatterns_pr <- custom_round(n_pr * prop_NA_pr)
spr_interp_pr_B1_30 <- foreach(i = 1:length(n_spatterns_pr), .combine = rbind,
                            .packages = c("bigmemory", "bigstatsr")) %dopar% {
  # Perform SPR for each n_spatterns
  number <- n_spatterns_pr
  interp_spr <- spr_interp(data_rcm_grid = data_pr_grid, n_spatterns = number,
                         data_obs_grid = data_pr_obs)
  # Back transformation for positivity
  interp_spr <- log(1 + exp(interp_spr)) # back transformation for positivity
  interp_spr[interp_spr < 0.5] <- 0 # lower bundary to remove noise
  return(interp_spr)
}

# Save the results in appropriate folder
save(spr_interp_pr_B1_30, file = "Projet_1/Interpolation_finale_3/Simulations_SPR/Region_B1/spr_interp_pr_B1_30.RData")


## Minimum temperature
data_tasmin_test <- data_tasmin_obs[, - test_indexes]
n_tasmin <- min(ncol(data_tasmin_grid), ncol(data_tasmin_test))
prop_NA_tasmin <- 0.30
n_spatterns_tasmin <- custom_round(n_tasmin * prop_NA_tasmin)
spr_interp_tasmin_B1_30 <- foreach(i = 1:length(n_spatterns_tasmin), .combine = rbind,
                            .packages = c("bigmemory", "bigstatsr")) %dopar% {
  # Perform SPR for each n_spatterns
  number <- n_spatterns_tasmin
  interp_spr <- spr_interp(data_rcm_grid = data_tasmin_grid, n_spatterns = number,
                         data_obs_grid = data_tasmin_obs)
  return(interp_spr)
}

# Save the results in appropriate folder
save(spr_interp_tasmin_B1_30, file = "Projet_1/Interpolation_finale_3/Simulations_SPR/Region_B1/spr_interp_tasmin_B1_30.RData")


## Maximum temperature
data_tasmax_test <- data_tasmax_obs[, - test_indexes]
n_tasmax <- min(ncol(data_tasmax_grid), ncol(data_tasmax_test))
prop_NA_tasmax <- 0.40
n_spatterns_tasmax <- custom_round(n_tasmax * prop_NA_tasmax)
spr_interp_tasmax_B1_30 <- foreach(i = 1:length(n_spatterns_tasmax), .combine = rbind,
                            .packages = c("bigmemory", "bigstatsr")) %dopar% {
  # Perform SPR for each n_spatterns
  number <- n_spatterns_tasmax
  interp_spr <- spr_interp(data_rcm_grid = data_tasmax_grid, n_spatterns = number,
                         data_obs_grid = data_tasmax_obs)
  return(interp_spr)
}

# Save the results in appropriate folder
save(spr_interp_tasmax_B1_30, file = "Projet_1/Interpolation_finale_3/Simulations_SPR/Region_B1/spr_interp_tasmax_B1_30.RData")

#######################################################################################
################################## 50% of missing data ################################


### RCM data
data_pr_grid <- data.pr[rcm_period, region]
data_tasmin_grid <- data.tasmin[rcm_period, region]
data_tasmax_grid <- data.tasmax[rcm_period, region]
# climatology.pr <- climatology_pr[, region]
# climatology.tasmin <- climatology_tasmin[, region]
# climatology.tasmax <- climatology_tasmax[, region]


### Observations data
data_pr_obs <- data.pr[obs_period, region]
### Different transformations for PRECIPITATION DATA
data_pr_obs[data_pr_obs <= 10^(-5)] <- 10^(-5) # lower bundary for positivity transformation
data_pr_obs <- log(exp(data_pr_obs) - 1) # ensuring positive values for precipitations
data_tasmin_obs <- data.tasmin[obs_period, region]
data_tasmax_obs <- data.tasmax[obs_period, region]


# # Extend the climatology over 10 years
# climatology_fun <- function(climatology, nyears) {
#   days_by_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
#   month_indexe <- rep(1:12, times = days_by_month)
#   climatology_year <- climatology[month_indexe, ]
#   climatology_nyears <- matrix(rep(climatology_year, times = nyears), nrow = 365 * nyears, byrow = TRUE)
#   return(climatology_nyears)
# }

# # Calculate anomalies for grid data
# anomalies_pr_grid <- data_pr_grid - climatology_fun(climatology.pr, nyears = 30)
# anomalies_tasmin_grid <- data_tasmin_grid - climatology_fun(climatology.tasmin, nyears = 30)
# anomalies_tasmax_grid <- data_tasmax_grid - climatology_fun(climatology.tasmax, nyears = 30)

# # Calculate anomalies for observations
# anomalies_pr_obs <- data_pr_obs - climatology_fun(climatology.pr, nyears = 10)
# anomalies_tasmin_obs <- data_tasmin_obs - climatology_fun(climatology.tasmin, nyears = 10)
# anomalies_tasmax_obs <- data_tasmax_obs - climatology_fun(climatology.tasmax, nyears = 10)

### Train and Test set
perc <- 0.5 # Change this percentage to represent each density
custom_round <- function(x) { # Custom rounding function
  ifelse(x - floor(x) == 0.5, ceiling(x), round(x))
}
n_test <- custom_round(perc * ncol(data_pr_obs)) # number of cells in test set
set.seed(seed) # set the seed
test_indexes <- sample(1:ncol(data_pr_obs), size = n_test, replace = FALSE)

### Extract data for test and train sets
#### data_..._obs contains the observed data with test indexes set to NA
data_pr_obs[, test_indexes] <- NA # set the test columns indexes to NA
data_tasmin_obs[, test_indexes] <- NA # set the test columns indexes to NA
data_tasmax_obs[, test_indexes] <- NA # set the test columns indexes to NA
# anomalies_pr_obs[, test_indexes] <- NA # set the test columns indexes to NA
# anomalies_tasmin_obs[, test_indexes] <- NA # set the test columns indexes to NA
# anomalies_tasmax_obs[, test_indexes] <- NA # set the test columns indexes to NA


### Now, apply SPR to data_..._val with different values of n_spatterns and
### return the number with the lower RMSE


## Precipitation
data_pr_test <- data_pr_obs[, - test_indexes]
n_pr <- min(ncol(data_pr_grid), ncol(data_pr_test))
prop_NA_pr <- 0.45
n_spatterns_pr <- custom_round(n_pr * prop_NA_pr)
spr_interp_pr_B1_50 <- foreach(i = 1:length(n_spatterns_pr), .combine = rbind,
                            .packages = c("bigmemory", "bigstatsr")) %dopar% {
  # Perform SPR for each n_spatterns
  number <- n_spatterns_pr
  interp_spr <- spr_interp(data_rcm_grid = data_pr_grid, n_spatterns = number,
                         data_obs_grid = data_pr_obs)
  # Back transformation for positivity
  interp_spr <- log(1 + exp(interp_spr)) # back transformation for positivity
  interp_spr[interp_spr < 0.5] <- 0 # lower bundary to remove noise
  return(interp_spr)
}

# Save the results in appropriate folder
save(spr_interp_pr_B1_50, file = "Projet_1/Interpolation_finale_3/Simulations_SPR/Region_B1/spr_interp_pr_B1_50.RData")


## Minimum temperature
data_tasmin_test <- data_tasmin_obs[, - test_indexes]
n_tasmin <- min(ncol(data_tasmin_grid), ncol(data_tasmin_test))
prop_NA_tasmin <- 0.20
n_spatterns_tasmin <- custom_round(n_tasmin * prop_NA_tasmin)
spr_interp_tasmin_B1_50 <- foreach(i = 1:length(n_spatterns_tasmin), .combine = rbind,
                            .packages = c("bigmemory", "bigstatsr")) %dopar% {
  # Perform SPR for each n_spatterns
  number <- n_spatterns_tasmin
  interp_spr <- spr_interp(data_rcm_grid = data_tasmin_grid, n_spatterns = number,
                         data_obs_grid = data_tasmin_obs)
  return(interp_spr)
}

# Save the results in appropriate folder
save(spr_interp_tasmin_B1_50, file = "Projet_1/Interpolation_finale_3/Simulations_SPR/Region_B1/spr_interp_tasmin_B1_50.RData")


## Maximum temperature
data_tasmax_test <- data_tasmax_obs[, - test_indexes]
n_tasmax <- min(ncol(data_tasmax_grid), ncol(data_tasmax_test))
prop_NA_tasmax <- 0.30
n_spatterns_tasmax <- custom_round(n_tasmax * prop_NA_tasmax)
spr_interp_tasmax_B1_50 <- foreach(i = 1:length(n_spatterns_tasmax), .combine = rbind,
                            .packages = c("bigmemory", "bigstatsr")) %dopar% {
  # Perform SPR for each n_spatterns
  number <- n_spatterns_tasmax
  interp_spr <- spr_interp(data_rcm_grid = data_tasmax_grid, n_spatterns = number,
                         data_obs_grid = data_tasmax_obs)
  return(interp_spr)
}

# Save the results in appropriate folder
save(spr_interp_tasmax_B1_50, file = "Projet_1/Interpolation_finale_3/Simulations_SPR/Region_B1/spr_interp_tasmax_B1_50.RData")

#######################################################################################
################################## 70% of missing data ################################


### RCM data
data_pr_grid <- data.pr[rcm_period, region]
data_tasmin_grid <- data.tasmin[rcm_period, region]
data_tasmax_grid <- data.tasmax[rcm_period, region]
# climatology.pr <- climatology_pr[, region]
# climatology.tasmin <- climatology_tasmin[, region]
# climatology.tasmax <- climatology_tasmax[, region]


### Observations data
data_pr_obs <- data.pr[obs_period, region]
### Different transformations for PRECIPITATION DATA
data_pr_obs[data_pr_obs <= 10^(-5)] <- 10^(-5) # lower bundary for positivity transformation
data_pr_obs <- log(exp(data_pr_obs) - 1) # ensuring positive values for precipitations
data_tasmin_obs <- data.tasmin[obs_period, region]
data_tasmax_obs <- data.tasmax[obs_period, region]


# # Extend the climatology over 10 years
# climatology_fun <- function(climatology, nyears) {
#   days_by_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
#   month_indexe <- rep(1:12, times = days_by_month)
#   climatology_year <- climatology[month_indexe, ]
#   climatology_nyears <- matrix(rep(climatology_year, times = nyears), nrow = 365 * nyears, byrow = TRUE)
#   return(climatology_nyears)
# }

# # Calculate anomalies for grid data
# anomalies_pr_grid <- data_pr_grid - climatology_fun(climatology.pr, nyears = 30)
# anomalies_tasmin_grid <- data_tasmin_grid - climatology_fun(climatology.tasmin, nyears = 30)
# anomalies_tasmax_grid <- data_tasmax_grid - climatology_fun(climatology.tasmax, nyears = 30)

# # Calculate anomalies for observations
# anomalies_pr_obs <- data_pr_obs - climatology_fun(climatology.pr, nyears = 10)
# anomalies_tasmin_obs <- data_tasmin_obs - climatology_fun(climatology.tasmin, nyears = 10)
# anomalies_tasmax_obs <- data_tasmax_obs - climatology_fun(climatology.tasmax, nyears = 10)

### Train and Test set
perc <- 0.7 # Change this percentage to represent each density
custom_round <- function(x) { # Custom rounding function
  ifelse(x - floor(x) == 0.5, ceiling(x), round(x))
}
n_test <- custom_round(perc * ncol(data_pr_obs)) # number of cells in test set
set.seed(seed) # set the seed
test_indexes <- sample(1:ncol(data_pr_obs), size = n_test, replace = FALSE)

### Extract data for test and train sets
#### data_..._obs contains the observed data with test indexes set to NA
data_pr_obs[, test_indexes] <- NA # set the test columns indexes to NA
data_tasmin_obs[, test_indexes] <- NA # set the test columns indexes to NA
data_tasmax_obs[, test_indexes] <- NA # set the test columns indexes to NA
# anomalies_pr_obs[, test_indexes] <- NA # set the test columns indexes to NA
# anomalies_tasmin_obs[, test_indexes] <- NA # set the test columns indexes to NA
# anomalies_tasmax_obs[, test_indexes] <- NA # set the test columns indexes to NA


### Now, apply SPR to data_..._val with different values of n_spatterns and
### return the number with the lower RMSE


## Precipitation
data_pr_test <- data_pr_obs[, - test_indexes]
n_pr <- min(ncol(data_pr_grid), ncol(data_pr_test))
prop_NA_pr <- 0.40
n_spatterns_pr <- custom_round(n_pr * prop_NA_pr)
spr_interp_pr_B1_70 <- foreach(i = 1:length(n_spatterns_pr), .combine = rbind,
                            .packages = c("bigmemory", "bigstatsr")) %dopar% {
  # Perform SPR for each n_spatterns
  number <- n_spatterns_pr
  interp_spr <- spr_interp(data_rcm_grid = data_pr_grid, n_spatterns = number,
                         data_obs_grid = data_pr_obs)
  # Back transformation for positivity
  interp_spr <- log(1 + exp(interp_spr)) # back transformation for positivity
  interp_spr[interp_spr < 0.5] <- 0 # lower bundary to remove noise
  return(interp_spr)
}

# Save the results in appropriate folder
save(spr_interp_pr_B1_70, file = "Projet_1/Interpolation_finale_3/Simulations_SPR/Region_B1/spr_interp_pr_B1_70.RData")


## Minimum temperature
data_tasmin_test <- data_tasmin_obs[, - test_indexes]
n_tasmin <- min(ncol(data_tasmin_grid), ncol(data_tasmin_test))
prop_NA_tasmin <- 0.25
n_spatterns_tasmin <- custom_round(n_tasmin * prop_NA_tasmin)
spr_interp_tasmin_B1_70 <- foreach(i = 1:length(n_spatterns_tasmin), .combine = rbind,
                            .packages = c("bigmemory", "bigstatsr")) %dopar% {
  # Perform SPR for each n_spatterns
  number <- n_spatterns_tasmin
  interp_spr <- spr_interp(data_rcm_grid = data_tasmin_grid, n_spatterns = number,
                         data_obs_grid = data_tasmin_obs)
  return(interp_spr)
}

# Save the results in appropriate folder
save(spr_interp_tasmin_B1_70, file = "Projet_1/Interpolation_finale_3/Simulations_SPR/Region_B1/spr_interp_tasmin_B1_70.RData")


## Maximum temperature
data_tasmax_test <- data_tasmax_obs[, - test_indexes]
n_tasmax <- min(ncol(data_tasmax_grid), ncol(data_tasmax_test))
prop_NA_tasmax <- 0.30
n_spatterns_tasmax <- custom_round(n_tasmax * prop_NA_tasmax)
spr_interp_tasmax_B1_70 <- foreach(i = 1:length(n_spatterns_tasmax), .combine = rbind,
                            .packages = c("bigmemory", "bigstatsr")) %dopar% {
  # Perform SPR for each n_spatterns
  number <- n_spatterns_tasmax
  interp_spr <- spr_interp(data_rcm_grid = data_tasmax_grid, n_spatterns = number,
                         data_obs_grid = data_tasmax_obs)
  return(interp_spr)
}

# Save the results in appropriate folder
save(spr_interp_tasmax_B1_70, file = "Projet_1/Interpolation_finale_3/Simulations_SPR/Region_B1/spr_interp_tasmax_B1_70.RData")

#######################################################################################
################################## 90% of missing data ################################


### RCM data
data_pr_grid <- data.pr[rcm_period, region]
data_tasmin_grid <- data.tasmin[rcm_period, region]
data_tasmax_grid <- data.tasmax[rcm_period, region]
# climatology.pr <- climatology_pr[, region]
# climatology.tasmin <- climatology_tasmin[, region]
# climatology.tasmax <- climatology_tasmax[, region]


### Observations data
data_pr_obs <- data.pr[obs_period, region]
### Different transformations for PRECIPITATION DATA
data_pr_obs[data_pr_obs <= 10^(-5)] <- 10^(-5) # lower bundary for positivity transformation
data_pr_obs <- log(exp(data_pr_obs) - 1) # ensuring positive values for precipitations
data_tasmin_obs <- data.tasmin[obs_period, region]
data_tasmax_obs <- data.tasmax[obs_period, region]


# # Extend the climatology over 10 years
# climatology_fun <- function(climatology, nyears) {
#   days_by_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
#   month_indexe <- rep(1:12, times = days_by_month)
#   climatology_year <- climatology[month_indexe, ]
#   climatology_nyears <- matrix(rep(climatology_year, times = nyears), nrow = 365 * nyears, byrow = TRUE)
#   return(climatology_nyears)
# }

# # Calculate anomalies for grid data
# anomalies_pr_grid <- data_pr_grid - climatology_fun(climatology.pr, nyears = 30)
# anomalies_tasmin_grid <- data_tasmin_grid - climatology_fun(climatology.tasmin, nyears = 30)
# anomalies_tasmax_grid <- data_tasmax_grid - climatology_fun(climatology.tasmax, nyears = 30)

# # Calculate anomalies for observations
# anomalies_pr_obs <- data_pr_obs - climatology_fun(climatology.pr, nyears = 10)
# anomalies_tasmin_obs <- data_tasmin_obs - climatology_fun(climatology.tasmin, nyears = 10)
# anomalies_tasmax_obs <- data_tasmax_obs - climatology_fun(climatology.tasmax, nyears = 10)

### Train and Test set
perc <- 0.9 # Change this percentage to represent each density
custom_round <- function(x) { # Custom rounding function
  ifelse(x - floor(x) == 0.5, ceiling(x), round(x))
}
n_test <- custom_round(perc * ncol(data_pr_obs)) # number of cells in test set
set.seed(seed) # set the seed
test_indexes <- sample(1:ncol(data_pr_obs), size = n_test, replace = FALSE)

### Extract data for test and train sets
#### data_..._obs contains the observed data with test indexes set to NA
data_pr_obs[, test_indexes] <- NA
data_tasmin_obs[, test_indexes] <- NA
data_tasmax_obs[, test_indexes] <- NA
# anomalies_pr_obs[, test_indexes] <- NA # set the test columns indexes to NA
# anomalies_tasmin_obs[, test_indexes] <- NA # set the test columns indexes to NA
# anomalies_tasmax_obs[, test_indexes] <- NA # set the test columns indexes to NA

### Now, apply SPR to data_..._val with different values of n_spatterns and
### return the number with the lower RMSE


## Precipitation
data_pr_test <- data_pr_obs[, - test_indexes]
n_pr <- min(ncol(data_pr_grid), ncol(data_pr_test))
prop_NA_pr <- 0.30
n_spatterns_pr <- custom_round(n_pr * prop_NA_pr)
spr_interp_pr_B1_90 <- foreach(i = 1:length(n_spatterns_pr), .combine = rbind,
                            .packages = c("bigmemory", "bigstatsr")) %dopar% {
  # Perform SPR for each n_spatterns
  number <- n_spatterns_pr
  interp_spr <- spr_interp(data_rcm_grid = data_pr_grid, n_spatterns = number,
                         data_obs_grid = data_pr_obs)
  # Back transformation for positivity
  interp_spr <- log(1 + exp(interp_spr)) # back transformation for positivity
  interp_spr[interp_spr < 0.5] <- 0 # lower bundary to remove noise
  return(interp_spr)
}

# Save the results in appropriate folder
save(spr_interp_pr_B1_90, file = "Projet_1/Interpolation_finale_3/Simulations_SPR/Region_B1/spr_interp_pr_B1_90.RData")


## Minimum temperature
data_tasmin_test <- data_tasmin_obs[, - test_indexes]
n_tasmin <- min(ncol(data_tasmin_grid), ncol(data_tasmin_test))
prop_NA_tasmin <- 0.30
n_spatterns_tasmin <- custom_round(n_tasmin * prop_NA_tasmin)
spr_interp_tasmin_B1_90 <- foreach(i = 1:length(n_spatterns_tasmin), .combine = rbind,
                            .packages = c("bigmemory", "bigstatsr")) %dopar% {
  # Perform SPR for each n_spatterns
  number <- n_spatterns_tasmin
  interp_spr <- spr_interp(data_rcm_grid = data_tasmin_grid, n_spatterns = number,
                         data_obs_grid = data_tasmin_obs)
  return(interp_spr)
}

# Save the results in appropriate folder
save(spr_interp_tasmin_B1_90, file = "Projet_1/Interpolation_finale_3/Simulations_SPR/Region_B1/spr_interp_tasmin_B1_90.RData")


## Maximum temperature
data_tasmax_test <- data_tasmax_obs[, - test_indexes]
n_tasmax <- min(ncol(data_tasmax_grid), ncol(data_tasmax_test))
prop_NA_tasmax <- 0.30
n_spatterns_tasmax <- custom_round(n_tasmax * prop_NA_tasmax)
spr_interp_tasmax_B1_90 <- foreach(i = 1:length(n_spatterns_tasmax), .combine = rbind,
                            .packages = c("bigmemory", "bigstatsr")) %dopar% {
  # Perform SPR for each n_spatterns
  number <- n_spatterns_tasmax
  interp_spr <- spr_interp(data_rcm_grid = data_tasmax_grid, n_spatterns = number,
                         data_obs_grid = data_tasmax_obs)
  return(interp_spr)
}

# Save the results in appropriate folder
save(spr_interp_tasmax_B1_90, file = "Projet_1/Interpolation_finale_3/Simulations_SPR/Region_B1/spr_interp_tasmax_B1_90.RData")


# Stop the cluster
stopCluster(cl)
