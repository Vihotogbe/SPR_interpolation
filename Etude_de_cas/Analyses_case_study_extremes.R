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
load("Projet_1/Analyses_finales/Etude_de_cas/spr_interp_pr_B1.RData")
load("Projet_1/Analyses_finales/Etude_de_cas/spr_interp_tasmin_B1.RData")
load("Projet_1/Analyses_finales/Etude_de_cas/spr_interp_tasmax_B1.RData")


# KED interpolated data
## REGION : B1
load("Projet_1/Analyses_finales/Etude_de_cas/ked_interp_pr_B1.RData")
load("Projet_1/Analyses_finales/Etude_de_cas/ked_interp_tasmin_B1.RData")
load("Projet_1/Analyses_finales/Etude_de_cas/ked_interp_tasmax_B1.RData")
ked_interp_pr_B1_10 <- as.matrix(ked_interp_pr_B1_10)
ked_interp_tasmin_B1_10 <- as.matrix(ked_interp_tasmin_B1_10)
ked_interp_tasmax_B1_10 <- as.matrix(ked_interp_tasmax_B1_10)
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

######################################### OBSERVATION ##############################################
#################################### ONLY 3 GAUGED STATIONS ########################################
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
###################################################################################

####### QUANTILE LOSS #######

# Étape 1 : fonction quantile loss
quantile_loss <- function(y_true, y_pred, tau = 0.90) {
  error <- y_true - y_pred
  loss <- ifelse(error >= 0, tau * error, (1 - tau) * -error)
  return(mean(loss, na.rm = TRUE))
}


### PRECIPITATION ###

# Étape 2 : Calcul journalier de la quantile loss
n_days <- nrow(data_pr_obs_0)
ql_spr <- numeric(n_days)
ql_ked <- numeric(n_days)

for (i in 1:n_days) {
  ql_spr[i] <- quantile_loss(data_pr_obs_0[i, ], spr_interp_pr_B1[i, test_indexes], tau = 0.95)
  ql_ked[i] <- quantile_loss(data_pr_obs_0[i, ], ked_interp_pr_B1_10[i, ], tau = 0.95)
}

# Étape 3 : Regrouper pour graphique
df_plot <- data.frame(
  #day = 1:n_days,
  SPR = ql_spr,
  KED = ql_ked
)

### Calculer les moyennes journalières
jour_annuel <- rep(1:365, times = 10)
df_plot$jour_annuel <- jour_annuel
# Calculer les moyennes journalières
moyennes_journalieres <- aggregate(. ~ jour_annuel, data = df_plot, FUN = mean)
head(moyennes_journalieres)
tail(moyennes_journalieres)

# Transformer pour ggplot
library(tidyr)
df_long <- pivot_longer(df_plot, cols = c("SPR", "KED"), names_to = "Model", values_to = "QuantileLoss")

# Étape 4 : Graphique
library(ggplot2)
ggplot(df_long, aes(x = jour_annuel, y = QuantileLoss, color = Model)) +
  geom_line(alpha = 0.7) +
  labs(title = "Quantile Loss (τ = 0.90) par jour pour chaque modèle",
       x = "Jour (sur 10 ans)",
       y = "Quantile Loss",
       color = "Modèle") +
  theme_minimal() +
  geom_smooth(se = FALSE, method = "loess") +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    strip.text = element_blank(),
    panel.grid = element_blank()
  )

##### MINIMUM TEMPERATURES
# Étape 2 : Calcul journalier de la quantile loss
n_days <- nrow(data_tasmin_obs_0)
ql_spr <- numeric(n_days)
ql_ked <- numeric(n_days)

for (i in 1:n_days) {
  ql_spr[i] <- quantile_loss(data_tasmin_obs_0[i, ], spr_interp_tasmin_B1[i, test_indexes], tau = 0.95)
  ql_ked[i] <- quantile_loss(data_tasmin_obs_0[i, ], ked_interp_tasmin_B1_10[i, ], tau = 0.95)
}

# Étape 3 : Regrouper pour graphique
df_plot <- data.frame(
  #day = 1:n_days,
  SPR = ql_spr,
  KED = ql_ked
)

### Calculer les moyennes journalières
jour_annuel <- rep(1:365, times = 10)
df_plot$jour_annuel <- jour_annuel
moyennes_journalieres <- aggregate(. ~ jour_annuel, data = df_plot, FUN = mean)
head(moyennes_journalieres)
tail(moyennes_journalieres)

# Transformer pour ggplot
library(tidyr)
df_long <- pivot_longer(df_plot, cols = c("SPR", "KED"), names_to = "Model", values_to = "QuantileLoss")
# Étape 4 : Graphique
library(ggplot2)
ggplot(df_long, aes(x = jour_annuel, y = QuantileLoss, color = Model)) +
  geom_line(alpha = 0.7) +
  labs(title = "Quantile Loss (τ = 0.90) par jour pour chaque modèle",
       x = "Jour (sur 10 ans)",
       y = "Quantile Loss",
       color = "Modèle") +
  theme_minimal() +
  geom_smooth(se = FALSE, method = "loess") +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    strip.text = element_blank(),
    panel.grid = element_blank()
  )


##### MAXIMUM TEMPERATURES
# Étape 2 : Calcul journalier de la quantile loss
n_days <- nrow(data_tasmax_obs_0)
ql_spr <- numeric(n_days)
ql_ked <- numeric(n_days)

for (i in 1:n_days) {
  ql_spr[i] <- quantile_loss(data_tasmax_obs_0[i, ], spr_interp_tasmax_B1[i, test_indexes], tau = 0.95)
  ql_ked[i] <- quantile_loss(data_tasmax_obs_0[i, ], ked_interp_tasmax_B1_10[i, ], tau = 0.95)
}

# Étape 3 : Regrouper pour graphique
df_plot <- data.frame(
  #day = 1:n_days,
  SPR = ql_spr,
  KED = ql_ked
)

### Calculer les moyennes journalières
jour_annuel <- rep(1:365, times = 10)
df_plot$jour_annuel <- jour_annuel
moyennes_journalieres <- aggregate(. ~ jour_annuel, data = df_plot, FUN = mean)
head(moyennes_journalieres)
tail(moyennes_journalieres)

# Transformer pour ggplot

library(tidyr)
df_long <- pivot_longer(df_plot, cols = c("SPR", "KED"), names_to = "Model", values_to = "QuantileLoss")

# Étape 4 : Graphique

library(ggplot2)
ggplot(df_long, aes(x = jour_annuel, y = QuantileLoss, color = Model)) +
  geom_line(alpha = 0.7) +
  labs(title = "Quantile Loss (τ = 0.90) par jour pour chaque modèle",
       x = "Jour (sur 10 ans)",
       y = "Quantile Loss",
       color = "Modèle") +
  theme_minimal() +
  geom_smooth(se = FALSE, method = "loess") +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    strip.text = element_blank(),
    panel.grid = element_blank()
  )


############################################################################
sum(moyennes_journalieres$SPR > moyennes_journalieres$KED) # 0
sum(moyennes_journalieres$SPR < moyennes_journalieres$KED) # 3650
sum(moyennes_journalieres$SPR == moyennes_journalieres$KED) # 0
###########################################################################
############### % of data > 95% percentile ################################
# Calculer le 95e percentile
percentile_95_obs_pr <- apply(data_pr_obs_0, 1, function(x) quantile(x, 0.95, na.rm = TRUE))
percentile_95_spr_pr <- apply(spr_interp_pr_B1[, test_indexes], 1, function(x) quantile(x, 0.95, na.rm = TRUE))
percentile_95_ked_pr <- apply(ked_interp_pr_B1_10, 1, function(x) quantile(x, 0.95, na.rm = TRUE))
head(percentile_95_obs_pr)
head(percentile_95_spr_pr)
head(percentile_95_ked_pr)
jour_annuel <- rep(1:365, times = 10)
data_percentile_pr <- data.frame(
  jour_annuel = jour_annuel,
  obs_pr = percentile_95_obs_pr,
  spr_pr = percentile_95_spr_pr,
  ked_pr = percentile_95_ked_pr
)
mean_day_pr <- aggregate(. ~ jour_annuel, data = data_percentile_pr, FUN = mean)
head(mean_day_pr)

# pr
range(data_pr_obs_0[1,])
range(spr_interp_pr_B1[, test_indexes][1,])
range(ked_interp_pr_B1_10[1,])
# tasmin
range(data_tasmin_obs_0[1,])
range(spr_interp_tasmin_B1[, test_indexes][1,])
range(ked_interp_tasmin_B1_10[1,])
# tasmax
range(data_tasmax_obs_0[1,])
range(spr_interp_tasmax_B1[, test_indexes][1,])
range(ked_interp_tasmax_B1_10[1,])
#######################################################################
############ AUTRE ANALYSE
i <- 1
threshold <- quantile(data_pr_obs_0[i,], 0.95)
obs_bin <- as.integer(data_pr_obs_0[i,] > threshold)
pred_bin_spr <- as.integer(spr_interp_pr_B1[i, test_indexes] > threshold)
pred_bin_ked <- as.integer(ked_interp_pr_B1_10[i,] > threshold)
confusionMatrix_spr <- table(obs_bin, pred_bin_spr)
confusionMatrix_spr
# Accuracy
accuracy_spr <- (confusionMatrix_spr[1, 1] + confusionMatrix_spr[2, 2]) / sum(confusionMatrix_spr) # accuracy
accuracy_spr
# Sensitivity
sensitivity_spr <- confusionMatrix_spr[2, 2] / (confusionMatrix_spr[2, 1] + confusionMatrix_spr[2, 2]) # true positive rate
sensitivity_spr
# Specificity
specificity_spr <- confusionMatrix_spr[1, 1] / (confusionMatrix_spr[1, 1] + confusionMatrix_spr[1, 2]) # true negative rate
specificity_spr

### KED
confusionMatrix_ked <- table(obs_bin, pred_bin_ked)
confusionMatrix_ked
# Accuracy
accuracy_ked <- (confusionMatrix_ked[1, 1] + confusionMatrix_ked[2, 2]) / sum(confusionMatrix_ked) # accuracy
accuracy_ked
# Sensitivity
sensitivity_ked <- confusionMatrix_ked[2, 2] / (confusionMatrix_ked[2, 1] + confusionMatrix_ked[2, 2]) # true positive rate
sensitivity_ked
# Specificity
specificity_ked <- confusionMatrix_ked[1, 1] / (confusionMatrix_ked[1, 1] + confusionMatrix_ked[1, 2]) # true negative rate
specificity_ked


spr_interp_pr_B1[1,1:10]
ked_interp_pr_B1_10[1,1:10]
data_pr_obs_0[1,1:10]
range(spr_interp_pr_B1[1,test_indexes])
range(ked_interp_pr_B1_10[1,])
range(data_pr_obs_0[1,])
#######################
threshold <- quantile(data_pr_obs_0, 0.95, na.rm = TRUE)
indices_extremes <- which(data_pr_obs_0 > threshold)
# MAE et RMSE sur les extrêmes
mae_extreme <- mean(abs(spr_interp_pr_B1[indices_extremes] - data_pr_obs_0[indices_extremes]), na.rm = TRUE)
rmse_extreme <- sqrt(mean((spr_interp_pr_B1[indices_extremes] - data_pr_obs_0[indices_extremes])^2, na.rm = TRUE))

# Pour comparaison : erreurs globales
mae_global <- mean(abs(spr_interp_pr_B1 - data_pr_obs_0), na.rm = TRUE)
rmse_global <- sqrt(mean((spr_interp_pr_B1 - data_pr_obs_0)^2, na.rm = TRUE))

# Affichage
cat("MAE global :", mae_global, "\n")
cat("MAE extrêmes :", mae_extreme, "\n")
cat("RMSE global :", rmse_global, "\n")
cat("RMSE extrêmes :", rmse_extreme, "\n")
pred_extreme <- which(spr_interp_pr_B1 > threshold)

# Vrai positif = intersection
true_extreme_detected <- length(intersect(indices_extremes, pred_extreme))
precision <- true_extreme_detected / length(pred_extreme)
recall <- true_extreme_detected / length(indices_extremes)

cat("Précision (prédits extrêmes qui sont vraiment extrêmes) :", precision, "\n")
cat("Rappel (vrais extrêmes bien prédits comme extrêmes) :", recall, "\n")
plot(data_pr_obs_0, spr_interp_pr_B1, pch = 16, col = "grey80",
     xlab = "Observé", ylab = "Prévu", main = "Valeurs extrêmes en rouge")
points(data_pr_obs_0[indices_extremes], spr_interp_pr_B1[indices_extremes], col = "red", pch = 16)
abline(a = 0, b = 1, col = "blue", lwd = 2)
qqplot(data_pr_obs_0, spr_interp_pr_B1,
       main = "QQ-plot Observé vs Prévu (queues)",
       xlab = "Quantiles observés", ylab = "Quantiles prédits")
abline(0, 1, col = "blue")
#########################################################################################
ked_interp_pr_B1_10[1,1:10]
ked_interp_pr_B1_10[2,1:10]
spr_interp_pr_B1[1,1:10]
spr_interp_pr_B1[2,1:10]
