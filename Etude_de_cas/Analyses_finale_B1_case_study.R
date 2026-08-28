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


# Defining the per:iod of RCM and observations data
rcm_period <- 1:10950 # 30 years of data (1980-2009)
obs_period_1 <- 7301:10950 # 10 years of data (2000-2009)
obs_period_2 <- 10951:14600 # 10 years of data (2010-2019)
obs_period_3 <- 14601:18250 # 10 years of data (2020-2029)


# Choose the period and the region
obs_period <- obs_period_1 # choose the period
region <- region_B1 # choose the region


seed <- 2244677 # Define the seed


####################################################################################################
########################################### CHOOSE A REGION ########################################
####################################################################################################

# Load the interpolated data
# SPR interpolated data
## REGION : B1
load("Projet_1/Interpolation_finale/Etude_de_cas/spr_interp_pr_B1.RData")
load("Projet_1/Interpolation_finale/Etude_de_cas/spr_interp_tasmin_B1.RData")
load("Projet_1/Interpolation_finale/Etude_de_cas/spr_interp_tasmax_B1.RData")

# KED interpolated data
## REGION : B1
load("Projet_1/Interpolation_finale/Etude_de_cas/ked_interp_pr_B1.RData")
load("Projet_1/Interpolation_finale/Etude_de_cas/ked_interp_tasmin_B1.RData")
load("Projet_1/Interpolation_finale/Etude_de_cas/ked_interp_tasmax_B1.RData")
######################################### CHOOSE THE REGION SIZE ####################################

# Taille de région
taille_region_A1 <- 3025
taille_region_A2 <- 1332
taille_region_A3 <- 342
taille_region_B1 <- 2970
taille_region_B2 <- 1332
taille_region_B3 <- 342

# Choix de la taille de la région
taille_region <- taille_region_B1 # choose the region size


####################################################################################################
############################# Function to calculate daily statistics ####################################

# Daily RMSE
rmse_day <- function(data_interp, data_obs) {
  rmse <- NULL
  for (i in 1:nrow(data_interp)) {
    rmse <- c(rmse, sqrt(mean((data_interp[i, ] - data_obs[i, ])^2)))
  }
  return(rmse)
}
# rmseee <- rmse_day(spr_interp_pr_B1_10[, test_indexes], data_pr_obs_0)

# Daily SSIM
ssim_global <- function(Z, Z_hat, K1 = 0.01, K2 = 0.03, L = NULL) {
  # Determine the dynamic range (L)
    L <- max(Z, Z_hat) - min(Z, Z_hat)  # Default to range of data
  
  # Compute means
  mu_x <- mean(Z)
  mu_y <- mean(Z_hat)
  
  # Compute variances
  sigma_x2 <- var(as.vector(Z))
  sigma_y2 <- var(as.vector(Z_hat))
  
  # Compute covariance
  sigma_xy <- cov(as.vector(Z), as.vector(Z_hat))
  
  # Compute constants
  C1 <- (K1 * L)^2
  C2 <- (K2 * L)^2
  
  # Compute SSIM
  ssim_value <- ((2 * mu_x * mu_y + C1) * (2 * sigma_xy + C2)) /
    ((mu_x^2 + mu_y^2 + C1) * (sigma_x2 + sigma_y2 + C2))
  
  return(ssim_value)
}


# Custom rounding
custom_round <- function(x) { # Custom rounding function
  ifelse(x - floor(x) == 0.5, ceiling(x), round(x))
}

# Calculate statistics
# Calculate mean, median, standard deviation, and confidence interval
calculate_stats <- function(rmse_values) {
  mean_val <- mean(rmse_values)
  median_val <- median(rmse_values)
  sd_val <- sd(rmse_values)
  ci_val <- quantile(rmse_values, probs = c(0.025, 0.975))
  return(round(c(mean_val, median_val, sd_val, ci_val), 2))
}

####################################################################################################
######################################### ONLY 3 GAUGED STATIONS ###################################
####################################################################################################

### Train and Test set indexes
perc <- 0.999 # Change this percentage to represent each density
custom_round <- function(x) { # Custom rounding function
  ifelse(x - floor(x) == 0.5, ceiling(x), round(x))
}
n_test <- custom_round(perc * taille_region) # number of cells in test set
set.seed(seed) # set the seed
test_indexes <- sample(1:taille_region, size = n_test, replace = FALSE)


# Data for comparison
data_pr_obs_0 <- data.pr[obs_period, region]
data_pr_obs_0 <- data_pr_obs_0[, test_indexes]
data_pr_obs_0[data_pr_obs_0 < 0.5] <- 0 # removing the lower values (noise) of pr 
data_tasmin_obs_0 <- data.tasmin[obs_period, region]
data_tasmin_obs_0 <- data_tasmin_obs_0[, test_indexes]
data_tasmax_obs_0 <- data.tasmax[obs_period, region]
data_tasmax_obs_0 <- data_tasmax_obs_0[, test_indexes]

### RMSE values
# SPR
v1 <- calculate_stats(rmse_day(spr_interp_pr_B1[, test_indexes], data_pr_obs_0))
v01 <- mean(ssim_global(spr_interp_pr_B1[, test_indexes], data_pr_obs_0))
v2 <- calculate_stats(rmse_day(spr_interp_tasmin_B1[, test_indexes], data_tasmin_obs_0))
v02 <- mean(ssim_global(spr_interp_tasmin_B1[, test_indexes], data_tasmin_obs_0))
v3 <- calculate_stats(rmse_day(spr_interp_tasmax_B1[, test_indexes], data_tasmax_obs_0))
v03 <- mean(ssim_global(spr_interp_tasmax_B1[, test_indexes], data_tasmax_obs_0))

# KED
v10 <- calculate_stats(rmse_day(ked_interp_pr_B1, data_pr_obs_0))
v010 <- mean(ssim_global(ked_interp_pr_B1, data_pr_obs_0))
v11 <- calculate_stats(rmse_day(ked_interp_tasmin_B1, data_tasmin_obs_0))
v011 <- mean(ssim_global(ked_interp_tasmin_B1, data_tasmin_obs_0))
v12 <- calculate_stats(rmse_day(ked_interp_tasmax_B1, data_tasmax_obs_0))
v012 <- mean(ssim_global(ked_interp_tasmax_B1, data_tasmax_obs_0))

######################################################################################################
################################### TABLES OF RMSE statistics ####################################

# pr
spr_ked_rmse_stats_pr <- as.matrix(rbind(v1, v10))
spr_ked_rmse_stats_pr <- cbind(spr_ked_rmse_stats_pr, round(c(v01, v010), 4))
colnames(spr_ked_rmse_stats_pr) <- c("Mean", "Median", "Sd", "2.5%", "97.5%", "SSIM")
rownames(spr_ked_rmse_stats_pr) <- c("SPR", "KED")

# tasmin
spr_ked_rmse_stats_tasmin <- as.matrix(rbind(v2, v11))
spr_ked_rmse_stats_tasmin <- cbind(spr_ked_rmse_stats_tasmin, round(c(v02, v011), 4))
colnames(spr_ked_rmse_stats_tasmin) <- c("Mean", "Median", "Sd", "2.5%", "97.5%", "SSIM")
rownames(spr_ked_rmse_stats_tasmin) <- c("SPR", "KED")

# tasmax
spr_ked_rmse_stats_tasmax <- as.matrix(rbind(v3, v12))
spr_ked_rmse_stats_tasmax <- cbind(spr_ked_rmse_stats_tasmax, round(c(v03, v012), 4))
colnames(spr_ked_rmse_stats_tasmax) <- c("Mean", "Median", "Sd", "2.5%", "97.5%", "SSIM")
rownames(spr_ked_rmse_stats_tasmax) <- c("SPR", "KED")

cat("#################### RMSE STATISTICS ####################\n")
cat("########################## PR ##########################\n")
print(spr_ked_rmse_stats_pr)
cat("########################## TASMIN ##########################\n")
print(spr_ked_rmse_stats_tasmin)
cat("########################## TASMAX ##########################\n")
print(spr_ked_rmse_stats_tasmax)
cat("########################## END ##########################\n")