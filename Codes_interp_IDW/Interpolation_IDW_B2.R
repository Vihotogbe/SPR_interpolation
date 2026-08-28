########################### IDW interpolation function ############################
# Charging libraries
library(gstat)
library(sf)
library(sp)
library(terra)

# IDW interpolation function
idw_interp <- function(data_train, coord_train, coord_test, power){
  data_idw.pred <- NULL
  for (i in 1:nrow(data_train)){
    data_Y <- data_train[i, ]
    data_idw <- data.frame(coord_train, values = data_Y)
    coordinates(data_idw) <- ~ lon + lat
    # Prediction locations
    grille <- data.frame(coord_test)
    coordinates(grille) <- ~ lon + lat
    # Interpolation IDW
    idw_result <- gstat(formula = values ~ 1, locations = data_idw, 
                        nmax = nrow(data_idw), set = list(idp = power))
    idw_result_pred <- predict(idw_result, grille)
    idw_pred <- idw_result_pred$var1.pred
    data_idw.pred <- rbind(data_idw.pred, idw_pred)
  }
  return(data_idw.pred)
}
################################# END OF THE FUNCTION #################################


########################## SIMULATIONS ###########################################

# LOADING DATA
load(file = "Data/data.coord.sf.RData") # rotated coord
load(file = "Data/data.pr.RData") # precipitation
load(file = "Data/data.tasmin.RData") # tasmin
load(file = "Data/data.tasmax.RData") # tasmax

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

# Defining the rotation CRS for the rotated coordinates
rotation_crs <- paste0("+proj=ob_tran +o_proj=longlat +o_lon_p=", 0,
                       " +o_lat_p=", 42.5, " +lon_0=", 83-180, " +to_meter=0.01745329")
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
region <- region_B2 # choose the region

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
data_coord <- data.coord.sf[region,]

# Convert rotated coordinates fron degrees to meters
data_coord <- st_as_sf(as.data.frame(data_coord), coords = c("lon", "lat"), crs = rotation_crs)
data_coord <- st_transform(data_coord, crs = 4326)
data_coord <- as.data.frame(st_coordinates(data_coord))
colnames(data_coord) <- c("lon", "lat") # renaming the columns

### Observations data
data_pr_obs <- data.pr[obs_period, region]
### Different transformations for PRECIPITATION DATA
data_pr_obs[data_pr_obs <= 10^(-5)] <- 10^(-5) # lower bundary for positivity transformation
data_pr_obs <- log(exp(data_pr_obs) - 1) # ensuring positive values for precipitations
data_tasmin_obs <- data.tasmin[obs_period, region]
data_tasmax_obs <- data.tasmax[obs_period, region]


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

### Define train and test sets
#### train sets
data_pr_train <- data_pr_obs[, - test_indexes]
data_tasmin_train <- data_tasmin_obs[, - test_indexes]
data_tasmax_train <- data_tasmax_obs[, - test_indexes]
coord_train <- data_coord[- test_indexes, ]
#### test sets
data_pr_test <- data_pr_obs[, test_indexes]
data_tasmin_test <- data_tasmin_obs[, test_indexes]
data_tasmax_test <- data_tasmax_obs[, test_indexes]
coord_test <- data_coord[test_indexes, ]


### Data for comparison in validation : because of positivity transformation
data_pr_obs_0 <- data.pr[obs_period, region]
data_pr_obs_0[data_pr_obs_0 < 0.5] <- 0 # removing the lower values as noise
data_tasmin_obs_0 <- data.tasmin[obs_period, region]
data_tasmax_obs_0 <- data.tasmax[obs_period, region]


### Now, apply IDW to data_..._train with the value of idw_power choosen
### return the interpolated data


## Precipitation
powers <- 5
idw_interp_pr_B2_10 <- foreach(i = 1:length(powers), .combine = rbind,
                          .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  interp_idw_pr <- idw_interp(data_train = data_pr_train, coord_train = coord_train,
                              coord_test = coord_test, power = powers)
  interp_idw_pr <- log(1 + exp(interp_idw_pr))
  interp_idw_pr[interp_idw_pr < 0.5] <- 0
  return(interp_idw_pr)
}

# Save the results in appropriate folder
save(idw_interp_pr_B2_10, file = "Projet_1/Interpolation_finale/Simulations_IDW/Region_B2/idw_interp_pr_B2_10.RData")


## Minimum temperature
powers <- 5
idw_interp_tasmin_B2_10 <- foreach(i = 1:length(powers), .combine = rbind,
                          .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  interp_idw_tasmin <- idw_interp(data_train = data_tasmin_train, coord_train = coord_train,
                                  coord_test = coord_test, power = powers)
  return(interp_idw_tasmin)
}

# Save the results in appropriate folder
save(idw_interp_tasmin_B2_10, file = "Projet_1/Interpolation_finale/Simulations_IDW/Region_B2/idw_interp_tasmin_B2_10.RData")


## Maximum temperature
powers <- 5
idw_interp_tasmax_B2_10 <- foreach(i = 1:length(powers), .combine = rbind,
                          .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  interp_idw_tasmax <- idw_interp(data_train = data_tasmax_train, coord_train = coord_train,
                                  coord_test = coord_test, power = powers)
  return(interp_idw_tasmax)
}

# Save the results in appropriate folder
save(idw_interp_tasmax_B2_10, file = "Projet_1/Interpolation_finale/Simulations_IDW/Region_B2/idw_interp_tasmax_B2_10.RData")


#######################################################################################
################################## 30% of missing data ################################

### RCM data
data_pr_grid <- data.pr[rcm_period, region]
data_tasmin_grid <- data.tasmin[rcm_period, region]
data_tasmax_grid <- data.tasmax[rcm_period, region]
data_coord <- data.coord.sf[region,]

# Convert rotated coordinates fron degrees to meters
data_coord <- st_as_sf(as.data.frame(data_coord), coords = c("lon", "lat"), crs = rotation_crs)
data_coord <- st_transform(data_coord, crs = 4326)
data_coord <- as.data.frame(st_coordinates(data_coord))
colnames(data_coord) <- c("lon", "lat") # renaming the columns

### Observations data
data_pr_obs <- data.pr[obs_period, region]
### Different transformations for PRECIPITATION DATA
data_pr_obs[data_pr_obs <= 10^(-5)] <- 10^(-5) # lower bundary for positivity transformation
data_pr_obs <- log(exp(data_pr_obs) - 1) # ensuring positive values for precipitations
data_tasmin_obs <- data.tasmin[obs_period, region]
data_tasmax_obs <- data.tasmax[obs_period, region]


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

### Define train and test sets
#### train sets
data_pr_train <- data_pr_obs[, - test_indexes]
data_tasmin_train <- data_tasmin_obs[, - test_indexes]
data_tasmax_train <- data_tasmax_obs[, - test_indexes]
coord_train <- data_coord[- test_indexes, ]
#### test sets
data_pr_test <- data_pr_obs[, test_indexes]
data_tasmin_test <- data_tasmin_obs[, test_indexes]
data_tasmax_test <- data_tasmax_obs[, test_indexes]
coord_test <- data_coord[test_indexes, ]


### Data for comparison in validation : because of positivity transformation
data_pr_obs_0 <- data.pr[obs_period, region]
data_pr_obs_0[data_pr_obs_0 < 0.5] <- 0 # removing the lower values as noise
data_tasmin_obs_0 <- data.tasmin[obs_period, region]
data_tasmax_obs_0 <- data.tasmax[obs_period, region]


### Now, apply IDW to data_..._train with the value of idw_power choosen
### return the interpolated data


## Precipitation
powers <- 5
idw_interp_pr_B2_30 <- foreach(i = 1:length(powers), .combine = rbind,
                          .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  interp_idw_pr <- idw_interp(data_train = data_pr_train, coord_train = coord_train,
                              coord_test = coord_test, power = powers)
  interp_idw_pr <- log(1 + exp(interp_idw_pr))
  interp_idw_pr[interp_idw_pr < 0.5] <- 0
  return(interp_idw_pr)
}

# Save the results in appropriate folder
save(idw_interp_pr_B2_30, file = "Projet_1/Interpolation_finale/Simulations_IDW/Region_B2/idw_interp_pr_B2_30.RData")


## Minimum temperature
powers <- 4
idw_interp_tasmin_B2_30 <- foreach(i = 1:length(powers), .combine = rbind,
                          .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  interp_idw_tasmin <- idw_interp(data_train = data_tasmin_train, coord_train = coord_train,
                                  coord_test = coord_test, power = powers)
  return(interp_idw_tasmin)
}

# Save the results in appropriate folder
save(idw_interp_tasmin_B2_30, file = "Projet_1/Interpolation_finale/Simulations_IDW/Region_B2/idw_interp_tasmin_B2_30.RData")


## Maximum temperature
powers <- 5
idw_interp_tasmax_B2_30 <- foreach(i = 1:length(powers), .combine = rbind,
                          .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  interp_idw_tasmax <- idw_interp(data_train = data_tasmax_train, coord_train = coord_train,
                                  coord_test = coord_test, power = powers)
  return(interp_idw_tasmax)
}

# Save the results in appropriate folder
save(idw_interp_tasmax_B2_30, file = "Projet_1/Interpolation_finale/Simulations_IDW/Region_B2/idw_interp_tasmax_B2_30.RData")


#######################################################################################
################################## 50% of missing data ################################

### RCM data
data_pr_grid <- data.pr[rcm_period, region]
data_tasmin_grid <- data.tasmin[rcm_period, region]
data_tasmax_grid <- data.tasmax[rcm_period, region]
data_coord <- data.coord.sf[region,]

# Convert rotated coordinates fron degrees to meters
data_coord <- st_as_sf(as.data.frame(data_coord), coords = c("lon", "lat"), crs = rotation_crs)
data_coord <- st_transform(data_coord, crs = 4326)
data_coord <- as.data.frame(st_coordinates(data_coord))
colnames(data_coord) <- c("lon", "lat") # renaming the columns

### Observations data
data_pr_obs <- data.pr[obs_period, region]
### Different transformations for PRECIPITATION DATA
data_pr_obs[data_pr_obs <= 10^(-5)] <- 10^(-5) # lower bundary for positivity transformation
data_pr_obs <- log(exp(data_pr_obs) - 1) # ensuring positive values for precipitations
data_tasmin_obs <- data.tasmin[obs_period, region]
data_tasmax_obs <- data.tasmax[obs_period, region]


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

### Define train and test sets
#### train sets
data_pr_train <- data_pr_obs[, - test_indexes]
data_tasmin_train <- data_tasmin_obs[, - test_indexes]
data_tasmax_train <- data_tasmax_obs[, - test_indexes]
coord_train <- data_coord[- test_indexes, ]
#### test sets
data_pr_test <- data_pr_obs[, test_indexes]
data_tasmin_test <- data_tasmin_obs[, test_indexes]
data_tasmax_test <- data_tasmax_obs[, test_indexes]
coord_test <- data_coord[test_indexes, ]


### Data for comparison in validation : because of positivity transformation
data_pr_obs_0 <- data.pr[obs_period, region]
data_pr_obs_0[data_pr_obs_0 < 0.5] <- 0 # removing the lower values as noise
data_tasmin_obs_0 <- data.tasmin[obs_period, region]
data_tasmax_obs_0 <- data.tasmax[obs_period, region]


### Now, apply IDW to data_..._train with the value of idw_power choosen
### return the interpolated data


## Precipitation
powers <- 4
idw_interp_pr_B2_50 <- foreach(i = 1:length(powers), .combine = rbind,
                          .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  interp_idw_pr <- idw_interp(data_train = data_pr_train, coord_train = coord_train,
                              coord_test = coord_test, power = powers)
  interp_idw_pr <- log(1 + exp(interp_idw_pr))
  interp_idw_pr[interp_idw_pr < 0.5] <- 0
  return(interp_idw_pr)
}

# Save the results in appropriate folder
save(idw_interp_pr_B2_50, file = "Projet_1/Interpolation_finale/Simulations_IDW/Region_B2/idw_interp_pr_B2_50.RData")


## Minimum temperature
powers <- 4
idw_interp_tasmin_B2_50 <- foreach(i = 1:length(powers), .combine = rbind,
                          .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  interp_idw_tasmin <- idw_interp(data_train = data_tasmin_train, coord_train = coord_train,
                                  coord_test = coord_test, power = powers)
  return(interp_idw_tasmin)
}

# Save the results in appropriate folder
save(idw_interp_tasmin_B2_50, file = "Projet_1/Interpolation_finale/Simulations_IDW/Region_B2/idw_interp_tasmin_B2_50.RData")


## Maximum temperature
powers <- 4
idw_interp_tasmax_B2_50 <- foreach(i = 1:length(powers), .combine = rbind,
                          .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  interp_idw_tasmax <- idw_interp(data_train = data_tasmax_train, coord_train = coord_train,
                                  coord_test = coord_test, power = powers)
  return(interp_idw_tasmax)
}

# Save the results in appropriate folder
save(idw_interp_tasmax_B2_50, file = "Projet_1/Interpolation_finale/Simulations_IDW/Region_B2/idw_interp_tasmax_B2_50.RData")


#######################################################################################
################################## 70% of missing data ################################

### RCM data
data_pr_grid <- data.pr[rcm_period, region]
data_tasmin_grid <- data.tasmin[rcm_period, region]
data_tasmax_grid <- data.tasmax[rcm_period, region]
data_coord <- data.coord.sf[region,]

# Convert rotated coordinates fron degrees to meters
data_coord <- st_as_sf(as.data.frame(data_coord), coords = c("lon", "lat"), crs = rotation_crs)
data_coord <- st_transform(data_coord, crs = 4326)
data_coord <- as.data.frame(st_coordinates(data_coord))
colnames(data_coord) <- c("lon", "lat") # renaming the columns

### Observations data
data_pr_obs <- data.pr[obs_period, region]
### Different transformations for PRECIPITATION DATA
data_pr_obs[data_pr_obs <= 10^(-5)] <- 10^(-5) # lower bundary for positivity transformation
data_pr_obs <- log(exp(data_pr_obs) - 1) # ensuring positive values for precipitations
data_tasmin_obs <- data.tasmin[obs_period, region]
data_tasmax_obs <- data.tasmax[obs_period, region]


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

### Define train and test sets
#### train sets
data_pr_train <- data_pr_obs[, - test_indexes]
data_tasmin_train <- data_tasmin_obs[, - test_indexes]
data_tasmax_train <- data_tasmax_obs[, - test_indexes]
coord_train <- data_coord[- test_indexes, ]
#### test sets
data_pr_test <- data_pr_obs[, test_indexes]
data_tasmin_test <- data_tasmin_obs[, test_indexes]
data_tasmax_test <- data_tasmax_obs[, test_indexes]
coord_test <- data_coord[test_indexes, ]


### Data for comparison in validation : because of positivity transformation
data_pr_obs_0 <- data.pr[obs_period, region]
data_pr_obs_0[data_pr_obs_0 < 0.5] <- 0 # removing the lower values as noise
data_tasmin_obs_0 <- data.tasmin[obs_period, region]
data_tasmax_obs_0 <- data.tasmax[obs_period, region]


### Now, apply IDW to data_..._train with the value of idw_power choosen
### return the interpolated data


## Precipitation
powers <- 4
idw_interp_pr_B2_70 <- foreach(i = 1:length(powers), .combine = rbind,
                          .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  interp_idw_pr <- idw_interp(data_train = data_pr_train, coord_train = coord_train,
                              coord_test = coord_test, power = powers)
  interp_idw_pr <- log(1 + exp(interp_idw_pr))
  interp_idw_pr[interp_idw_pr < 0.5] <- 0
  return(interp_idw_pr)
}

# Save the results in appropriate folder
save(idw_interp_pr_B2_70, file = "Projet_1/Interpolation_finale/Simulations_IDW/Region_B2/idw_interp_pr_B2_70.RData")


## Minimum temperature
powers <- 3
idw_interp_tasmin_B2_70 <- foreach(i = 1:length(powers), .combine = rbind,
                          .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  interp_idw_tasmin <- idw_interp(data_train = data_tasmin_train, coord_train = coord_train,
                                  coord_test = coord_test, power = powers)
  return(interp_idw_tasmin)
}

# Save the results in appropriate folder
save(idw_interp_tasmin_B2_70, file = "Projet_1/Interpolation_finale/Simulations_IDW/Region_B2/idw_interp_tasmin_B2_70.RData")


## Maximum temperature
powers <- 4
idw_interp_tasmax_B2_70 <- foreach(i = 1:length(powers), .combine = rbind,
                          .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  interp_idw_tasmax <- idw_interp(data_train = data_tasmax_train, coord_train = coord_train,
                                  coord_test = coord_test, power = powers)
  return(interp_idw_tasmax)
}

# Save the results in appropriate folder
save(idw_interp_tasmax_B2_70, file = "Projet_1/Interpolation_finale/Simulations_IDW/Region_B2/idw_interp_tasmax_B2_70.RData")


#######################################################################################
################################## 90% of missing data ################################

### RCM data
data_pr_grid <- data.pr[rcm_period, region]
data_tasmin_grid <- data.tasmin[rcm_period, region]
data_tasmax_grid <- data.tasmax[rcm_period, region]
data_coord <- data.coord.sf[region,]

# Convert rotated coordinates fron degrees to meters
data_coord <- st_as_sf(as.data.frame(data_coord), coords = c("lon", "lat"), crs = rotation_crs)
data_coord <- st_transform(data_coord, crs = 4326)
data_coord <- as.data.frame(st_coordinates(data_coord))
colnames(data_coord) <- c("lon", "lat") # renaming the columns

### Observations data
data_pr_obs <- data.pr[obs_period, region]
### Different transformations for PRECIPITATION DATA
data_pr_obs[data_pr_obs <= 10^(-5)] <- 10^(-5) # lower bundary for positivity transformation
data_pr_obs <- log(exp(data_pr_obs) - 1) # ensuring positive values for precipitations
data_tasmin_obs <- data.tasmin[obs_period, region]
data_tasmax_obs <- data.tasmax[obs_period, region]


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
data_pr_obs[, test_indexes] <- NA # set the test columns indexes to NA
data_tasmin_obs[, test_indexes] <- NA # set the test columns indexes to NA
data_tasmax_obs[, test_indexes] <- NA # set the test columns indexes to NA

### Define train and test sets
#### train sets
data_pr_train <- data_pr_obs[, - test_indexes]
data_tasmin_train <- data_tasmin_obs[, - test_indexes]
data_tasmax_train <- data_tasmax_obs[, - test_indexes]
coord_train <- data_coord[- test_indexes, ]
#### test sets
data_pr_test <- data_pr_obs[, test_indexes]
data_tasmin_test <- data_tasmin_obs[, test_indexes]
data_tasmax_test <- data_tasmax_obs[, test_indexes]
coord_test <- data_coord[test_indexes, ]


### Data for comparison in validation : because of positivity transformation
data_pr_obs_0 <- data.pr[obs_period, region]
data_pr_obs_0[data_pr_obs_0 < 0.5] <- 0 # removing the lower values as noise
data_tasmin_obs_0 <- data.tasmin[obs_period, region]
data_tasmax_obs_0 <- data.tasmax[obs_period, region]


### Now, apply IDW to data_..._train with the value of idw_power choosen
### return the interpolated data


## Precipitation
powers <- 3
idw_interp_pr_B2_90 <- foreach(i = 1:length(powers), .combine = rbind,
                          .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  interp_idw_pr <- idw_interp(data_train = data_pr_train, coord_train = coord_train,
                              coord_test = coord_test, power = powers)
  interp_idw_pr <- log(1 + exp(interp_idw_pr))
  interp_idw_pr[interp_idw_pr < 0.5] <- 0
  return(interp_idw_pr)
}

# Save the results in appropriate folder
save(idw_interp_pr_B2_90, file = "Projet_1/Interpolation_finale/Simulations_IDW/Region_B2/idw_interp_pr_B2_90.RData")


## Minimum temperature
powers <- 2
idw_interp_tasmin_B2_90 <- foreach(i = 1:length(powers), .combine = rbind,
                          .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  interp_idw_tasmin <- idw_interp(data_train = data_tasmin_train, coord_train = coord_train,
                                  coord_test = coord_test, power = powers)
  return(interp_idw_tasmin)
}

# Save the results in appropriate folder
save(idw_interp_tasmin_B2_90, file = "Projet_1/Interpolation_finale/Simulations_IDW/Region_B2/idw_interp_tasmin_B2_90.RData")


## Maximum temperature
powers <- 3
idw_interp_tasmax_B2_90 <- foreach(i = 1:length(powers), .combine = rbind,
                          .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  interp_idw_tasmax <- idw_interp(data_train = data_tasmax_train, coord_train = coord_train,
                                  coord_test = coord_test, power = powers)
  return(interp_idw_tasmax)
}

# Save the results in appropriate folder
save(idw_interp_tasmax_B2_90, file = "Projet_1/Interpolation_finale/Simulations_IDW/Region_B2/idw_interp_tasmax_B2_90.RData")



# Stop the cluster
stopCluster(cl)
