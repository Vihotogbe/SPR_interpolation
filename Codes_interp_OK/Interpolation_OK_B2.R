########################### OK interpolation function ############################
# Charging libraries
library(gstat)
library(sf)
library(sp)
library(terra)


# OK interpolation function
ok_interp <- function(data_train, coord_train, coord_test, model) {
  data_ok.pred <- NULL
  # Loop for each location
  for (i in 1:nrow(data_train)) {
    data_Y <- data_train[i, ]
    data_ok <- data.frame(coord_train, values = data_Y)
    coordinates(data_ok) <- ~ lon + lat
    
    # Define test grid
    grille <- data.frame(coord_test)
    coordinates(grille) <- ~ lon + lat
    
    # Compute empirical variogram
    vgm_emp <- variogram(values ~ 1, data_ok, cutoff = max(dist(coord_train)) * 0.75, width = 0.2)
    
    # Initialize variogram model
    vgm_init <- vgm(model, 
                     psill = max(median(vgm_emp$gamma, na.rm = TRUE), 1e-6), 
                     range = max(vgm_emp$dist) / 3, 
                     nugget = max(min(vgm_emp$gamma, na.rm = TRUE), 1e-6))
    
    # Fit variogram model with error handling
    vgm_fit <- tryCatch({
      fit.variogram(vgm_emp, vgm_init)
    }, error = function(e) {
      vgm_init # Use initial values if fitting fails
    })
    
    # Ensure range is positive
    if (vgm_fit$range[2] <= 0) {
      vgm_fit$range[2] <- max(vgm_emp$dist) / 3
    }
    
    # Ensure psill is positive
    if (vgm_fit$psill[2] <= 0) {
      vgm_fit$psill[2] <- 1e-6
    }
    
    # Ordinary Kriging interpolation
    krig_result <- gstat(formula = values ~ 1, locations = data_ok, model = vgm_fit)
    krig_result_pred <- predict(krig_result, grille)
    
    # Extract predictions
    ok_pred <- krig_result_pred$var1.pred
    data_ok.pred <- rbind(data_ok.pred, ok_pred)
  }
  
  return(data_ok.pred)
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
rotation_crs <- paste0("+proj=ob_tran +o_proj=longlat +o_lon_p=", 0,
                       " +o_lat_p=", 42.5, " +lon_0=", 83-180, " +to_meter=0.01745329")
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



### Now, apply KED to data_..._train with the value of ok_power choosen
### return the interpolated data


## Precipitation
model <- "Sph"
ok_interp_pr_B2_10 <- foreach(i = 1:length(model), .combine = rbind,
                            .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  # Perform IDW for each model
  model <- model
  interp_ok_pr <- ok_interp(data_train = data_pr_train, coord_train = coord_train,
                              coord_test = coord_test, model = model)
  interp_ok_pr <- log(1 + exp(interp_ok_pr))
  interp_ok_pr[interp_ok_pr < 0.5] <- 0
  return(interp_ok_pr)
}

# Save the results in appropriate folder
save(ok_interp_pr_B2_10, file = "Projet_1/Interpolation_finale/Simulations_OK/Region_B2/ok_interp_pr_B2_10.RData")


## Minimum temperature
model <- "Sph"
ok_interp_tasmin_B2_10 <- foreach(i = 1:length(model), .combine = rbind,
                            .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  # Perform IDW for each model
  model <- model
  interp_ok_tasmin <- ok_interp(data_train = data_tasmin_train, coord_train = coord_train,
                              coord_test = coord_test, model = model)
  return(interp_ok_tasmin)
}

# Save the results in appropriate folder
save(ok_interp_tasmin_B2_10, file = "Projet_1/Interpolation_finale/Simulations_OK/Region_B2/ok_interp_tasmin_B2_10.RData")


## Maximum temperature
model <- "Sph"
ok_interp_tasmax_B2_10 <- foreach(i = 1:length(model), .combine = rbind,
                            .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  # Perform IDW for each model
  model <- model
  interp_ok_tasmax <- ok_interp(data_train = data_tasmax_train, coord_train = coord_train,
                              coord_test = coord_test, model = model)
  return(interp_ok_tasmax)
}

# Save the results in appropriate folder
save(ok_interp_tasmax_B2_10, file = "Projet_1/Interpolation_finale/Simulations_OK/Region_B2/ok_interp_tasmax_B2_10.RData")


#######################################################################################
################################## 30% of missing data ################################


### RCM data
data_pr_grid <- data.pr[rcm_period, region]
data_tasmin_grid <- data.tasmin[rcm_period, region]
data_tasmax_grid <- data.tasmax[rcm_period, region]
data_coord <- data.coord.sf[region,]

# Convert rotated coordinates fron degrees to meters
rotation_crs <- paste0("+proj=ob_tran +o_proj=longlat +o_lon_p=", 0,
                       " +o_lat_p=", 42.5, " +lon_0=", 83-180, " +to_meter=0.01745329")
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



### Now, apply KED to data_..._train with the value of ok_power choosen
### return the interpolated data


## Precipitation
model <- "Sph"
ok_interp_pr_B2_30 <- foreach(i = 1:length(model), .combine = rbind,
                            .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  # Perform IDW for each model
  model <- model
  interp_ok_pr <- ok_interp(data_train = data_pr_train, coord_train = coord_train,
                              coord_test = coord_test, model = model)
  interp_ok_pr <- log(1 + exp(interp_ok_pr))
  interp_ok_pr[interp_ok_pr < 0.5] <- 0
  return(interp_ok_pr)
}

# Save the results in appropriate folder
save(ok_interp_pr_B2_30, file = "Projet_1/Interpolation_finale/Simulations_OK/Region_B2/ok_interp_pr_B2_30.RData")


## Minimum temperature
model <- "Sph"
ok_interp_tasmin_B2_30 <- foreach(i = 1:length(model), .combine = rbind,
                            .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  # Perform IDW for each model
  model <- model
  interp_ok_tasmin <- ok_interp(data_train = data_tasmin_train, coord_train = coord_train,
                              coord_test = coord_test, model = model)
  return(interp_ok_tasmin)
}

# Save the results in appropriate folder
save(ok_interp_tasmin_B2_30, file = "Projet_1/Interpolation_finale/Simulations_OK/Region_B2/ok_interp_tasmin_B2_30.RData")


## Maximum temperature
model <- "Sph"
ok_interp_tasmax_B2_30 <- foreach(i = 1:length(model), .combine = rbind,
                            .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  # Perform IDW for each model
  model <- model
  interp_ok_tasmax <- ok_interp(data_train = data_tasmax_train, coord_train = coord_train,
                              coord_test = coord_test, model = model)
  return(interp_ok_tasmax)
}

# Save the results in appropriate folder
save(ok_interp_tasmax_B2_30, file = "Projet_1/Interpolation_finale/Simulations_OK/Region_B2/ok_interp_tasmax_B2_30.RData")

#######################################################################################
################################## 50% of missing data ################################


### RCM data
data_pr_grid <- data.pr[rcm_period, region]
data_tasmin_grid <- data.tasmin[rcm_period, region]
data_tasmax_grid <- data.tasmax[rcm_period, region]
data_coord <- data.coord.sf[region,]

# Convert rotated coordinates fron degrees to meters
rotation_crs <- paste0("+proj=ob_tran +o_proj=longlat +o_lon_p=", 0,
                       " +o_lat_p=", 42.5, " +lon_0=", 83-180, " +to_meter=0.01745329")
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


### Now, apply KED to data_..._train with the value of ok_power choosen
### return the interpolated data


## Precipitation
model <- "Sph"
ok_interp_pr_B2_50 <- foreach(i = 1:length(model), .combine = rbind,
                            .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  # Perform IDW for each model
  model <- model
  interp_ok_pr <- ok_interp(data_train = data_pr_train, coord_train = coord_train,
                              coord_test = coord_test, model = model)
  interp_ok_pr <- log(1 + exp(interp_ok_pr))
  interp_ok_pr[interp_ok_pr < 0.5] <- 0
  return(interp_ok_pr)
}

# Save the results in appropriate folder
save(ok_interp_pr_B2_50, file = "Projet_1/Interpolation_finale/Simulations_OK/Region_B2/ok_interp_pr_B2_50.RData")


## Minimum temperature
model <- "Sph"
ok_interp_tasmin_B2_50 <- foreach(i = 1:length(model), .combine = rbind,
                            .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  # Perform IDW for each model
  model <- model
  interp_ok_tasmin <- ok_interp(data_train = data_tasmin_train, coord_train = coord_train,
                              coord_test = coord_test, model = model)
  return(interp_ok_tasmin)
}

# Save the results in appropriate folder
save(ok_interp_tasmin_B2_50, file = "Projet_1/Interpolation_finale/Simulations_OK/Region_B2/ok_interp_tasmin_B2_50.RData")


## Maximum temperature
model <- "Sph"
ok_interp_tasmax_B2_50 <- foreach(i = 1:length(model), .combine = rbind,
                            .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  # Perform IDW for each model
  model <- model
  interp_ok_tasmax <- ok_interp(data_train = data_tasmax_train, coord_train = coord_train,
                              coord_test = coord_test, model = model)
  return(interp_ok_tasmax)
}

# Save the results in appropriate folder
save(ok_interp_tasmax_B2_50, file = "Projet_1/Interpolation_finale/Simulations_OK/Region_B2/ok_interp_tasmax_B2_50.RData")

#######################################################################################
################################## 70% of missing data ################################


### RCM data
data_pr_grid <- data.pr[rcm_period, region]
data_tasmin_grid <- data.tasmin[rcm_period, region]
data_tasmax_grid <- data.tasmax[rcm_period, region]
data_coord <- data.coord.sf[region,]

# Convert rotated coordinates fron degrees to meters
rotation_crs <- paste0("+proj=ob_tran +o_proj=longlat +o_lon_p=", 0,
                       " +o_lat_p=", 42.5, " +lon_0=", 83-180, " +to_meter=0.01745329")
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


### Now, apply KED to data_..._train with the value of ok_power choosen
### return the interpolated data


## Precipitation
model <- "Sph"
ok_interp_pr_B2_70 <- foreach(i = 1:length(model), .combine = rbind,
                            .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  # Perform IDW for each model
  model <- model
  interp_ok_pr <- ok_interp(data_train = data_pr_train, coord_train = coord_train,
                              coord_test = coord_test, model = model)
  interp_ok_pr <- log(1 + exp(interp_ok_pr))
  interp_ok_pr[interp_ok_pr < 0.5] <- 0
  return(interp_ok_pr)
}

# Save the results in appropriate folder
save(ok_interp_pr_B2_70, file = "Projet_1/Interpolation_finale/Simulations_OK/Region_B2/ok_interp_pr_B2_70.RData")


## Minimum temperature
model <- "Sph"
ok_interp_tasmin_B2_70 <- foreach(i = 1:length(model), .combine = rbind,
                            .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  # Perform IDW for each model
  model <- model
  interp_ok_tasmin <- ok_interp(data_train = data_tasmin_train, coord_train = coord_train,
                              coord_test = coord_test, model = model)
  return(interp_ok_tasmin)
}

# Save the results in appropriate folder
save(ok_interp_tasmin_B2_70, file = "Projet_1/Interpolation_finale/Simulations_OK/Region_B2/ok_interp_tasmin_B2_70.RData")


## Maximum temperature
model <- "Sph"
ok_interp_tasmax_B2_70 <- foreach(i = 1:length(model), .combine = rbind,
                            .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  # Perform IDW for each model
  model <- model
  interp_ok_tasmax <- ok_interp(data_train = data_tasmax_train, coord_train = coord_train,
                              coord_test = coord_test, model = model)
  return(interp_ok_tasmax)
}

# Save the results in appropriate folder
save(ok_interp_tasmax_B2_70, file = "Projet_1/Interpolation_finale/Simulations_OK/Region_B2/ok_interp_tasmax_B2_70.RData")


#######################################################################################
################################## 90% of missing data ################################


### RCM data
data_pr_grid <- data.pr[rcm_period, region]
data_tasmin_grid <- data.tasmin[rcm_period, region]
data_tasmax_grid <- data.tasmax[rcm_period, region]
data_coord <- data.coord.sf[region,]

# Convert rotated coordinates fron degrees to meters
rotation_crs <- paste0("+proj=ob_tran +o_proj=longlat +o_lon_p=", 0,
                       " +o_lat_p=", 42.5, " +lon_0=", 83-180, " +to_meter=0.01745329")
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


### Now, apply KED to data_..._train with the value of ok_power choosen
### return the interpolated data


## Precipitation
model <- "Exp"
ok_interp_pr_B2_90 <- foreach(i = 1:length(model), .combine = rbind,
                            .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  # Perform IDW for each model
  model <- model
  interp_ok_pr <- ok_interp(data_train = data_pr_train, coord_train = coord_train,
                              coord_test = coord_test, model = model)
  interp_ok_pr <- log(1 + exp(interp_ok_pr))
  interp_ok_pr[interp_ok_pr < 0.5] <- 0
  return(interp_ok_pr)
}

# Save the results in appropriate folder
save(ok_interp_pr_B2_90, file = "Projet_1/Interpolation_finale/Simulations_OK/Region_B2/ok_interp_pr_B2_90.RData")


## Minimum temperature
model <- "Exp"
ok_interp_tasmin_B2_90 <- foreach(i = 1:length(model), .combine = rbind,
                            .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  # Perform IDW for each model
  model <- model
  interp_ok_tasmin <- ok_interp(data_train = data_tasmin_train, coord_train = coord_train,
                              coord_test = coord_test, model = model)
  return(interp_ok_tasmin)
}

# Save the results in appropriate folder
save(ok_interp_tasmin_B2_90, file = "Projet_1/Interpolation_finale/Simulations_OK/Region_B2/ok_interp_tasmin_B2_90.RData")


## Maximum temperature
model <- "Exp"
ok_interp_tasmax_B2_90 <- foreach(i = 1:length(model), .combine = rbind,
                            .packages = c("gstat", "sf", "sp", "terra")) %dopar% {
  # Perform IDW for each model
  model <- model
  interp_ok_tasmax <- ok_interp(data_train = data_tasmax_train, coord_train = coord_train,
                              coord_test = coord_test, model = model)
  return(interp_ok_tasmax)
}

# Save the results in appropriate folder
save(ok_interp_tasmax_B2_90, file = "Projet_1/Interpolation_finale/Simulations_OK/Region_B2/ok_interp_tasmax_B2_90.RData")


# Stop the cluster
stopCluster(cl)
