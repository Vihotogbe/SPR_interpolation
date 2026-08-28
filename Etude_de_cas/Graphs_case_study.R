########################## ANALYSIS
library(ggplot2)
library(ggspatial)
library(maps)
library(tmap)
library(sf)
library(sp)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(dplyr)
theme_set(theme_light(18))
library(RColorBrewer)
library(gridExtra)
library(scales)
library(cowplot)

# LOADING DATA
load(file = "Data/data.coord.sf.RData") # rotated coord
load(file = "Data/data.pr.RData") # precipitation
load(file = "Data/data.tasmin.RData") # tasmin
load(file = "Data/data.tasmax.RData") # tasmax
load("Data/climatology_pr.RData")
load("Data/climatology_tasmin.RData")
load("Data/climatology_tasmax.RData")

# REGIONS A and B DETERMINATION : Showing the 3 sizes of region
rlon <- data.coord.sf[, "lon"]
rlat <- data.coord.sf[, "lat"]
region_B1 <- which(rlat >= 7 & rlat <= 13 & rlon >= 14 & rlon <= 20, arr.ind = TRUE)


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


# Choix de la taille de la région
taille_region <- 2970 # choose the region size

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
# data_pr_obs_0 <- data_pr_obs_0[, test_indexes]
data_pr_obs_0[data_pr_obs_0 < 0.5] <- 0 # removing the lower values (noise) of pr 
data_tasmin_obs_0 <- data.tasmin[obs_period, region]
# data_tasmin_obs_0 <- data_tasmin_obs_0[, test_indexes]
data_tasmax_obs_0 <- data.tasmax[obs_period, region]
# data_tasmax_obs_0 <- data_tasmax_obs_0[, test_indexes]
###################################################################################

############ SPR data
data_pr_interp_spr <- spr_interp_pr_B1#[, test_indexes]
data_tasmin_interp_spr <- spr_interp_tasmin_B1#[, test_indexes]
data_tasmax_interp_spr <- spr_interp_tasmax_B1#[, test_indexes]

############ KED data
data_pr_interp_ked <- data_pr_obs_0 # Use the observed data for KED
data_pr_interp_ked[, test_indexes] <- ked_interp_pr_B1 # Replace the observed data with KED interpolated data
# data_pr_interp_ked <- ked_interp_pr_B1_10
data_tasmin_interp_ked <- data_tasmin_obs_0 # Use the observed data for KED
data_tasmin_interp_ked[, test_indexes] <- ked_interp_tasmin_B1 # Replace the observed data with KED interpolated data
data_tasmax_interp_ked <- data_tasmax_obs_0 # Use the observed data for KED
data_tasmax_interp_ked[, test_indexes] <- ked_interp_tasmax_B1 # Replace the observed data with KED interpolated data

# Give a day and get the month of the corresponding climatology
day_to_month <- function(day) { # Supposing a 10 years period
  # Monthly boundaries for a 365-day year
  debut_mois <- c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)
  fin_mois <- c(30, 58, 89, 119, 150, 180, 211, 242, 272, 303, 333, 364)
  # The day
  i <- day
  nyears <- (i-1) %/% 365 # Get the year, 0 is the first year, 1 is the second year, etc.
  aa <- nyears * 365 + 1 # Get the first day of the year
  indice_mois <- which((i - aa) >= debut_mois & (i - aa) <= fin_mois) # Get the month
  return(indice_mois)
}
####################################################################################################
##################################### PLOTING DAILY VALUES FOR SOME DAYS ##########################
####################################################################################################
####################################################### MAPS

## Define colors palettes
### Temperature
temp_colors <- c(
  rgb(5, 48, 97, maxColorValue = 255),   # dark blue
  rgb(33, 102, 172, maxColorValue = 255),
  rgb(67, 147, 195, maxColorValue = 255),
  rgb(146, 197, 222, maxColorValue = 255),
  rgb(209, 229, 240, maxColorValue = 255),
  rgb(247, 247, 247, maxColorValue = 255),
  rgb(253, 219, 199, maxColorValue = 255),
  rgb(244, 165, 130, maxColorValue = 255),
  rgb(214, 96, 77, maxColorValue = 255),
  rgb(178, 24, 43, maxColorValue = 255),
  rgb(103, 0, 31, maxColorValue = 255)    # dark red
)
palette_temp <- colorRampPalette(temp_colors) # Create a color palette function
num_colors_temp <- length(temp_colors) # Number of discrete colors
discrete_palette_temp <- palette_temp(num_colors_temp) # Generate the discrete color palette
### Precipitation
palette_precip <- brewer.pal(9, "Blues")

#### The crs for rotated coordinate system
rotation_crs <- paste0("+proj=ob_tran +o_proj=longlat +o_lon_p=", 0,
                " +o_lat_p=", 42.5, " +lon_0=", 83-180, " +to_meter=0.01745329")

#### Charge the world, rivers and lakes shapefiles
world <- ne_countries(scale = 50, returnclass = "sf")
rivers <- ne_load(scale = 10, type = "rivers_lake_centerlines", destdir = "maps", returnclass = "sf")
lakes <- ne_load(scale = 10, type = "lakes", destdir = "maps", returnclass = "sf")

### Coordinates for the map
data_coord <- data.coord.sf[region, ]
coord_train <- data_coord#[-test_indexes, ] # training coordinates
# Reprojeter dans le même système que la carte principale
points_train_sf <- data_coord[-test_indexes, ]
points_train_sf <- st_as_sf(points_train_sf, coords = c("lon", "lat"), crs = rotation_crs)
######################################################################################
# PRECIPITATION
### Maps for the 3 variables for a given day + for climatology
### Choose a wet day for precipitation and dry day for temperature
xx <- rowSums(data_pr_obs_0)
day_pr <- 3519 #1883 # which.max(xx) + 500

## PRECIPITATION
rcoord.sf <- data_coord#[test_indexes, ]
data.graph <- data.frame(rcoord.sf, values_obs = data_pr_obs_0[day_pr, ],
                          values_spr = data_pr_interp_spr[day_pr, ],
                          values_ked = data_pr_interp_ked[day_pr, ]) # add climatology values
### Tranforming to match the crs of the rotated coordinate system
data.graph.sf <- st_as_sf(data.graph, coords = c("lon", "lat"), crs = rotation_crs)

# Values for the color scale for precipitation
min_break <- round(min(data.graph[, c(3, 4, 5)]))
max_break <- round(max(data.graph[, c(3, 4, 5)]))
mean_break <- round(0.5 * (min_break + max_break))

#### LEGEND
p <- ggplot() +
  geom_sf(data = data.graph.sf, aes(color = values_obs), size = 6, shape = 15) +  # Apply color to points
  scale_colour_gradientn(colors = palette_precip,
          limits = c(min_break, max_break+0.5),
          breaks = c(min_break, mean_break, max_break),
          guide = guide_colourbar(barheight = 8, barwidth = 0.5, title.position = "right"),
                         name = "pr (mm/day)") +
  theme(legend.title = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 15),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "right",
        legend.margin = margin(t = -30),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  geom_sf(data = world, fill = "transparent", color = "black", linewidth = 0.2) +
  coord_sf(crs = rotation_crs, xlim = c(14.25, 19.75), ylim = c(7.38, 12.7)) +
  labs(x = "", y = "")
# Extract legend
legend <- get_legend(p)

ggsave("Projet_1/Interpolation_finale/Etude_de_cas/Graphs/legend_pr_B1.png", width = 1,
       height = 4, bg = "white", plot = cowplot::plot_grid(legend))
#############################################################################################

### OBSERVED
ggplot() +
  geom_sf(data = data.graph.sf, aes(color = values_obs), size = 4, shape = 15) +  # Apply color to points
  scale_colour_gradientn(
    colors = palette_precip,
    limits = c(min_break, max_break),
    breaks = c(min_break, mean_break, max_break),
    guide = guide_colourbar(barheight = 0.5, barwidth = 8, title.position = "bottom", 
                            title.hjust = 0.5, direction = "horizontal"), name = "pr (mm/day)") +
  geom_sf(data = points_train_sf, shape = 22, color = "black", fill = "green", size = 4, stroke = 1) +
  theme(
    legend.title = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 15),
    legend.text = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.margin = margin(t = -30),
    plot.margin = unit(c(0.5, 0, 0, 0), "cm")
  ) +
  geom_sf(data = world, fill = "transparent", color = "black", linewidth = 0.2) +
  coord_sf(crs = rotation_crs, xlim = c(14.25, 19.75), ylim = c(7.38, 12.7)) +
  labs(x = "", y = "")

### Saving the plot
ggsave("Projet_1/Interpolation_finale/Etude_de_cas/Graphs/day_pr_B1_obs.png", width = 4,
       height = 4, bg = "white")

### SPR
ggplot() +
  geom_sf(data = data.graph.sf, aes(color = values_spr), size = 4, shape = 15) +  # Apply color to points
  scale_colour_gradientn(colors = palette_precip,
          limits = c(min_break, max_break),
          breaks = c(min_break, mean_break, max_break),
          guide = guide_colourbar(barheight = 8, barwidth = 0.5, title.position = "right"),
                         name = "pr (mm/day)") +
  geom_sf(data = points_train_sf, shape = 22, color = "black", fill = "green", size = 4, stroke = 1) +
  theme(legend.title = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 10),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "none",
        legend.margin = margin(t = -30),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  geom_sf(data = world, fill = "transparent", color = "black", linewidth = 0.2) +
  coord_sf(crs = rotation_crs, xlim = c(14.25, 19.75), ylim = c(7.38, 12.7)) +
  labs(x = "", y = "")

### Saving the plot
ggsave("Projet_1/Interpolation_finale/Etude_de_cas/Graphs/day_pr_B1_spr.png", width = 4,
       height = 4, bg = "white")

### KED
ggplot() +
  geom_sf(data = data.graph.sf, aes(color = values_ked), size = 4, shape = 15) +  # Apply color to points
  scale_colour_gradientn(colors = palette_precip,
          limits = c(min_break, max_break+0.5),
          breaks = c(min_break, mean_break, max_break),
          guide = guide_colourbar(barheight = 8, barwidth = 0.5, title.position = "right"),
                         name = "pr (mm/day)") +
  geom_sf(data = points_train_sf, shape = 22, color = "black", fill = "green", size = 4, stroke = 1) +
  theme(legend.title = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 10),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "none",
        legend.margin = margin(t = -30),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  geom_sf(data = world, fill = "transparent", color = "black", linewidth = 0.2) +
  coord_sf(crs = rotation_crs, xlim = c(14.25, 19.75), ylim = c(7.38, 12.7)) +
  labs(x = "", y = "")

### Saving the plot
ggsave("Projet_1/Interpolation_finale/Etude_de_cas/Graphs/day_pr_B1_ked.png", width = 4,
       height = 4, bg = "white")


## MINIMUM TEMPERATURE
### Choose cold day for tasmin
xx <- rowSums(data_tasmin_obs_0)
day_tasmin <- 100 # which.min(xx) - 1500

# Data for the map of tasmin
rcoord.sf <- data_coord#[test_indexes, ]
data.graph_1 <- data.frame(rcoord.sf, values_obs = data_tasmin_obs_0[day_tasmin, ],
                          values_spr = data_tasmin_interp_spr[day_tasmin, ],
                          values_ked = data_tasmin_interp_ked[day_tasmin, ]) # change the day
### Choose a dry day for tasmax
xx <- rowSums(data_tasmax_obs_0)
day_tasmax <- 100 #2045 #which.max(xx) - 500 #3491 # 2800 and 2050

# Data for the map of tasmax
data.graph_2 <- data.frame(rcoord.sf, values_obs = data_tasmax_obs_0[day_tasmax, ],
                          values_spr = data_tasmax_interp_spr[day_tasmax, ],
                          values_ked = data_tasmax_interp_ked[day_tasmax, ]) # change the day

### Tranforming to match the crs of the rotated coordinate system
data.graph.sf <- st_as_sf(data.graph_1, coords = c("lon", "lat"), crs = rotation_crs)
# Values for the color scale for temperature
min_break <- round(min(data.graph_1[, c(3, 4, 5)]))
max_break <- round(max(data.graph_2[, c(3, 4, 5)]))
mean_break <- round(0.5 * (min_break + max_break))
###
#### LEGEND
p <- ggplot() +
  geom_sf(data = data.graph.sf, aes(color = values_obs), size = 6, shape = 15) +  # Apply color to points
  scale_colour_gradientn(colors = discrete_palette_temp,
          limits = c(min_break, max_break+0.5),
           breaks = c(min_break, mean_break, max_break),
          guide = guide_colourbar(barheight = 8, barwidth = 0.5, title.position = "right"),
                         name = "tmin (°C)") +
  theme(legend.title = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 15),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "right",
        legend.margin = margin(t = -30),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  geom_sf(data = world, fill = "transparent", color = "black", linewidth = 0.2) +
  coord_sf(crs = rotation_crs, xlim = c(14.25, 19.75), ylim = c(7.38, 12.7)) +
  labs(x = "", y = "")

# Extract legend
legend <- get_legend(p)

ggsave("Projet_1/Interpolation_finale/Etude_de_cas/Graphs/legend_temp_tasmin_B1.png", width = 1,
       height = 4, bg = "white", plot = cowplot::plot_grid(legend))
#############################################################################################

### OBSERVED
ggplot() +
  geom_sf(data = data.graph.sf, aes(color = values_obs), size = 4, shape = 15) +  # Apply color to points
  scale_colour_gradientn(colors = discrete_palette_temp,
          limits = c(min_break, max_break),
           breaks = c(min_break, mean_break, max_break),
          guide = guide_colourbar(barheight = 8, barwidth = 0.5, title.position = "right"),
                         name = "tasmin (°C)") +
  geom_sf(data = points_train_sf, shape = 22, color = "black", fill = "green", size = 4, stroke = 1) +
  theme(legend.title = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 15),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "none",
        legend.margin = margin(t = -30),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  geom_sf(data = world, fill = "transparent", color = "black", linewidth = 0.2) +
  coord_sf(crs = rotation_crs, xlim = c(14.25, 19.75), ylim = c(7.38, 12.7)) +
  labs(x = "", y = "")
### Saving the plot
ggsave("Projet_1/Interpolation_finale/Etude_de_cas/Graphs/day_tasmin_B1_obs.png", width = 4,
       height = 4, bg = "white")

### SPR
ggplot() +
  geom_sf(data = data.graph.sf, aes(color = values_spr), size = 4, shape = 15) +  # Apply color to points
  scale_colour_gradientn(colors = discrete_palette_temp,
          limits = c(min_break, max_break+0.5),
           breaks = c(min_break, mean_break, max_break),
          guide = guide_colourbar(barheight = 8, barwidth = 0.5, title.position = "right"),
                         name = "tasmin (°C)") +
  geom_sf(data = points_train_sf, shape = 22, color = "black", fill = "green", size = 4, stroke = 1) +  # train points
  theme(legend.title = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 10),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "none",
        legend.margin = margin(t = -30),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  geom_sf(data = world, fill = "transparent", color = "black", linewidth = 0.2) +
  coord_sf(crs = rotation_crs, xlim = c(14.25, 19.75), ylim = c(7.38, 12.7)) +
  labs(x = "", y = "")

### Saving the plot 
ggsave("Projet_1/Interpolation_finale/Etude_de_cas/Graphs/day_tasmin_B1_spr.png", width = 4,
       height = 4, bg = "white")

### KED
ggplot() +
  geom_sf(data = data.graph.sf, aes(color = values_ked), size = 4, shape = 15) +  # Apply color to points
  scale_colour_gradientn(colors = discrete_palette_temp,
          limits = c(min_break, max_break+0.5),
           breaks = c(min_break, mean_break, max_break),
          guide = guide_colourbar(barheight = 8, barwidth = 0.5, title.position = "right"),
                         name = "tasmin (°C)") +
  geom_sf(data = points_train_sf, shape = 22, color = "black", fill = "green", size = 4, stroke = 1) +  # train points
  theme(legend.title = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 10),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "none",
        legend.margin = margin(t = -30),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  geom_sf(data = world, fill = "transparent", color = "black", linewidth = 0.2) +
  coord_sf(crs = rotation_crs, xlim = c(14.25, 19.75), ylim = c(7.38, 12.7)) +
  labs(x = "", y = "")
### Saving the plot
ggsave("Projet_1/Interpolation_finale/Etude_de_cas/Graphs/day_tasmin_B1_ked.png", width = 4,
       height = 4, bg = "white")


## MAXIMUM TEMPERATURE

#### LEGEND
p <- ggplot() +
  geom_sf(data = data.graph.sf, aes(color = values_obs), size = 6, shape = 15) +  # Apply color to points
  scale_colour_gradientn(colors = discrete_palette_temp,
          limits = c(min_break, max_break+0.5),
           breaks = c(min_break, mean_break, max_break),
          guide = guide_colourbar(barheight = 8, barwidth = 0.5, title.position = "right"),
                         name = "tmax (°C)") +
  theme(legend.title = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 15),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "right",
        legend.margin = margin(t = -30),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  geom_sf(data = world, fill = "transparent", color = "black", linewidth = 0.2) +
  coord_sf(crs = rotation_crs, xlim = c(14.25, 19.75), ylim = c(7.38, 12.7)) +
  labs(x = "", y = "")
  
# Extract legend
legend <- get_legend(p)

ggsave("Projet_1/Interpolation_finale/Etude_de_cas/Graphs/legend_temp_tasmax_B1.png", width = 1,
       height = 4, bg = "white", plot = cowplot::plot_grid(legend))

### Tranforming to match the crs of the rotated coordinate system
data.graph.sf <- st_as_sf(data.graph_2, coords = c("lon", "lat"), crs = rotation_crs)
# Values for the color scale for temperature
min_break <- round(min(data.graph_1[, c(3, 4, 5)]))
max_break <- round(max(data.graph_2[, c(3, 4, 5)]))
mean_break <- round(0.5 * (min_break + max_break))

### OBSERVED
ggplot() +
  geom_sf(data = data.graph.sf, aes(color = values_obs), size = 4, shape = 15) +  # Apply color to points
  scale_colour_gradientn(colors = discrete_palette_temp,
          limits = c(min_break, max_break),
           breaks = c(min_break, mean_break, max_break),
          guide = guide_colourbar(barheight = 8, barwidth = 0.5, title.position = "right"),
                         name = "tasmax (°C)") +
  geom_sf(data = points_train_sf, shape = 22, color = "black", fill = "green", size = 4, stroke = 1) +  # train points
  theme(legend.title = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 15),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "none",
        legend.margin = margin(t = -30),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  geom_sf(data = world, fill = "transparent", color = "black", linewidth = 0.2) +
  coord_sf(crs = rotation_crs, xlim = c(14.25, 19.75), ylim = c(7.38, 12.7)) +
  labs(x = "", y = "")

### Saving the plot
ggsave("Projet_1/Interpolation_finale/Etude_de_cas/Graphs/day_tasmax_B1_obs.png", width = 4,
       height = 4, bg = "white")

### SPR
ggplot() +
  geom_sf(data = data.graph.sf, aes(color = values_spr), size = 4, shape = 15) +  # Apply color to points
  scale_colour_gradientn(colors = discrete_palette_temp,
          limits = c(min_break, max_break+0.5),
           breaks = c(min_break, mean_break, max_break),
          guide = guide_colourbar(barheight = 8, barwidth = 0.5, title.position = "right"),
                         name = "tasmax (°C)") +
  geom_sf(data = points_train_sf, shape = 22, color = "black", fill = "green", size = 4, stroke = 1) +  # train points
  theme(legend.title = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 10),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "none",
        legend.margin = margin(t = -30),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  geom_sf(data = world, fill = "transparent", color = "black", linewidth = 0.2) +
  coord_sf(crs = rotation_crs, xlim = c(14.25, 19.75), ylim = c(7.38, 12.7)) +
  labs(x = "", y = "")

### Saving the plot
ggsave("Projet_1/Interpolation_finale/Etude_de_cas/Graphs/day_tasmax_B1_spr.png", width = 4,
       height = 4, bg = "white")

### KED
ggplot() +
  geom_sf(data = data.graph.sf, aes(color = values_ked), size = 4, shape = 15) +  # Apply color to points
  scale_colour_gradientn(colors = discrete_palette_temp,
          limits = c(min_break, max_break+0.5),
           breaks = c(min_break, mean_break, max_break),
          guide = guide_colourbar(barheight = 8, barwidth = 0.5, title.position = "right"),
                         name = "tasmax (°C)") +
  geom_sf(data = points_train_sf, shape = 22, color = "black", fill = "green", size = 4, stroke = 1) +  # train points
  theme(legend.title = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 10),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "none",
        legend.margin = margin(t = -30),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  geom_sf(data = world, fill = "transparent", color = "black", linewidth = 0.2) +
  coord_sf(crs = rotation_crs, xlim = c(14.25, 19.75), ylim = c(7.38, 12.7)) +
  labs(x = "", y = "")

### Saving the plot
ggsave("Projet_1/Interpolation_finale/Etude_de_cas/Graphs/day_tasmax_B1_ked.png", width = 4,
       height = 4, bg = "white")

##########################################################################################
##########################################################################################

##########################################################################################
######################### MEAN ABSOLUTE ERROR ON THE MAP FOR THE DAY #####################
##########################################################################################
## PLOT for each variable for the choosed day
### Choose a day for precipitation, tasmin and tasmax in order to have a good representation
### of SPR over KED
day_pr <- 3519 # which.max(xx_pr) + 500
day_tasmin <- 100 # which.min(xx_tasmin) - 1500
day_tasmax <-  100 # which.max(xx_tasmax) - 1500

## Calculating RMSE FOR EACH VARIABLE AND THE CHOSED DAY
### RMSE for precipitation
rmse_pr_spr <- sqrt((data_pr_obs_0[day_pr, test_indexes] - data_pr_interp_spr[day_pr, test_indexes])^2)
rmse_pr_ked <- sqrt((data_pr_obs_0[day_pr, test_indexes] - data_pr_interp_ked[day_pr, test_indexes])^2)
rmse_pr_spr <- round(rmse_pr_spr)
rmse_pr_ked <- round(rmse_pr_ked)

### RMSE for tasmin
rmse_tasmin_spr <- sqrt((data_tasmin_obs_0[day_tasmin, test_indexes] - data_tasmin_interp_spr[day_tasmin, test_indexes])^2)
rmse_tasmin_ked <- sqrt((data_tasmin_obs_0[day_tasmin, test_indexes] - data_tasmin_interp_ked[day_tasmin, test_indexes])^2)
rmse_tasmin_spr <- round(rmse_tasmin_spr)
rmse_tasmin_ked <- round(rmse_tasmin_ked)

### RMSE for tasmax
rmse_tasmax_spr <- sqrt((data_tasmax_obs_0[day_tasmax, test_indexes] - data_tasmax_interp_spr[day_tasmax, test_indexes])^2)
rmse_tasmax_ked <- sqrt((data_tasmax_obs_0[day_tasmax, test_indexes] - data_tasmax_interp_ked[day_tasmax, test_indexes])^2)
rmse_tasmax_spr <- round(rmse_tasmax_spr)
rmse_tasmax_ked <- round(rmse_tasmax_ked)

# Colors palette for RMSE
rmse_colors <- c("white","orange", "red", "darkred")
palette_rmse <- colorRampPalette(rmse_colors)
num_colors_rmse <- length(rmse_colors)
discrete_palette_rmse <- palette_rmse(num_colors_rmse) # Generate the discrete color palette

## PRECIPITATION
### Breaks for the color scale
min_break <- round(min(rmse_pr_spr, rmse_pr_ked))
max_break <- round(max(rmse_pr_spr, rmse_pr_ked))
mean_break <- round(0.5 * (min_break + max_break))

# Data for the map
rcoord.sf <- data_coord[test_indexes, ]
data.graph <- data.frame(rcoord.sf, values_spr = rmse_pr_spr, values_ked = rmse_pr_ked)
### Tranforming to match the crs of the rotated coordinate system
data.graph.sf <- st_as_sf(data.graph, coords = c("lon", "lat"), crs = rotation_crs)

###################################################### LEGEND
p <- ggplot() +
  geom_sf(data = data.graph.sf, aes(color = values_spr), size = 6, shape = 15) +  # Apply color to points
  scale_colour_gradientn(colors = discrete_palette_rmse,
          limits = c(min_break, max_break+0.5),
           breaks = c(min_break, mean_break, max_break),
          guide = guide_colourbar(barheight = 8, barwidth = 0.5, title.position = "right"),
                         name = "pr (mm/day)") +
  theme(legend.title = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 15),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "right",
        legend.margin = margin(t = -30),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  geom_sf(data = world, fill = "transparent", color = "black", linewidth = 0.2) +
  coord_sf(crs = rotation_crs, xlim = c(14.25, 19.75), ylim = c(7.38, 12.7)) +
  labs(x = "", y = "")
# Extract legend
legend <- cowplot::get_legend(p)
# Plot the legend
plot_grid(legend, ncol = 1, rel_heights = c(1))
# Saving the legend
ggsave("Projet_1/Interpolation_finale/Etude_de_cas/Graphs/legend_rmse_pr_B1.png", width = 1,
       height = 4, bg = "white", plot = cowplot::plot_grid(legend))
#############################################################################################

### SPR
ggplot() +
  geom_sf(data = data.graph.sf, aes(color = values_spr), size = 4, shape = 15) +  # Apply color to points
  scale_colour_gradientn(colors = discrete_palette_rmse,
          limits = c(min_break, max_break+0.5),
           breaks = c(min_break, mean_break, max_break),
          guide = guide_colourbar(barheight = 8, barwidth = 0.5, title.position = "right"),
                         name = "RMSE (mm/day)") +
  geom_sf(data = points_train_sf, shape = 22, color = "black", fill = "green", size = 4, stroke = 1) +  # train points
  theme(legend.title = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 15),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "none",
        legend.margin = margin(t = -30),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  geom_sf(data = world, fill = "transparent", color = "black", linewidth = 0.2) +
  coord_sf(crs = rotation_crs, xlim = c(14.25, 19.75), ylim = c(7.38, 12.7)) +
  labs(x = "", y = "")

### Saving the plot
ggsave("Projet_1/Interpolation_finale/Etude_de_cas/Graphs/rmse_pr_B1_spr.png", width = 4,
       height = 4, bg = "white")

### KED
ggplot() +
  geom_sf(data = data.graph.sf, aes(color = values_ked), size = 4, shape = 15) +  # Apply color to points
  scale_colour_gradientn(colors = discrete_palette_rmse,
          limits = c(min_break, max_break+0.5),
           breaks = c(min_break, mean_break, max_break),
          guide = guide_colourbar(barheight = 8, barwidth = 0.5, title.position = "right"),
                         name = "RMSE (mm/day)") +
   geom_sf(data = points_train_sf, shape = 22, color = "black", fill = "green", size = 4, stroke = 1) +  # train points
  theme(legend.title = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 10),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "none",
        legend.margin = margin(t = -30),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  geom_sf(data = world, fill = "transparent", color = "black", linewidth = 0.2) +
  coord_sf(crs = rotation_crs, xlim = c(14.25, 19.75), ylim = c(7.38, 12.7)) +
  labs(x = "", y = "")

### Saving the plot
ggsave("Projet_1/Interpolation_finale/Etude_de_cas/Graphs/rmse_pr_B1_ked.png", width = 4,
       height = 4, bg = "white")


## MINIMUM TEMPERATURE
### Breaks for the color scale
min_break <- round(min(rmse_tasmin_spr, rmse_tasmin_ked))
max_break <- round(max(rmse_tasmin_spr, rmse_tasmin_ked))
mean_break <- round(0.5 * (min_break + max_break))

# Data for the map
rcoord.sf <- data_coord[test_indexes, ]
data.graph <- data.frame(rcoord.sf, values_spr = rmse_tasmin_spr, values_ked = rmse_tasmin_ked)
### Tranforming to match the crs of the rotated coordinate system
data.graph.sf <- st_as_sf(data.graph, coords = c("lon", "lat"), crs = rotation_crs)

########################### LEGEND
p <- ggplot() +
  geom_sf(data = data.graph.sf, aes(color = values_spr), size = 6, shape = 15) +  # Apply color to points
  scale_colour_gradientn(colors = discrete_palette_rmse,
          limits = c(min_break, max_break+0.5),
           breaks = c(min_break, mean_break, max_break),
          guide = guide_colourbar(barheight = 8, barwidth = 0.5, title.position = "right"),
                         name = "tmin (°C)") +
  theme(legend.title = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 15),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "right",
        legend.margin = margin(t = -30),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  geom_sf(data = world, fill = "transparent", color = "black", linewidth = 0.2) +
  coord_sf(crs = rotation_crs, xlim = c(14.25, 19.75), ylim = c(7.38, 12.7)) +
  labs(x = "", y = "")
# Extract legend
legend <- get_legend(p)

ggsave("Projet_1/Interpolation_finale/Etude_de_cas/Graphs/legend_rmse_tasmin_B1.png", width = 1,
       height = 4, bg = "white", plot = cowplot::plot_grid(legend))
#############################################################################################

### SPR
ggplot() +
  geom_sf(data = data.graph.sf, aes(color = values_spr), size = 4, shape = 15) +  # Apply color to points
  scale_colour_gradientn(colors = discrete_palette_rmse,
          limits = c(min_break, max_break+0.5),
           breaks = c(min_break, mean_break, max_break),
          guide = guide_colourbar(barheight = 8, barwidth = 0.5, title.position = "right"),
                         name = "RMSE (°C)") +
  geom_sf(data = points_train_sf, shape = 22, color = "black", fill = "green", size = 4, stroke = 1) +  # train points
  theme(legend.title = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 15),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "none",
        legend.margin = margin(t = -30),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  geom_sf(data = world, fill = "transparent", color = "black", linewidth = 0.2) +
  coord_sf(crs = rotation_crs, xlim = c(14.25, 19.75), ylim = c(7.38, 12.7)) +
  labs(x = "", y = "")

### Saving the plot
ggsave("Projet_1/Interpolation_finale/Etude_de_cas/Graphs/rmse_tasmin_B1_spr.png", width = 4,
       height = 4, bg = "white")

### KED
ggplot() +
  geom_sf(data = data.graph.sf, aes(color = values_ked), size = 4, shape = 15) +  # Apply color to points
  scale_colour_gradientn(colors = discrete_palette_rmse,
          limits = c(min_break, max_break+0.5),
           breaks = c(min_break, mean_break, max_break),
          guide = guide_colourbar(barheight = 8, barwidth = 0.5, title.position = "right"),
                         name = "RMSE (°C)") +
  geom_sf(data = points_train_sf, shape = 22, color = "black", fill = "green", size = 4, stroke = 1) +  # train points
  theme(legend.title = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 10),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "none",
        legend.margin = margin(t = -30),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  geom_sf(data = world, fill = "transparent", color = "black", linewidth = 0.2) +
  coord_sf(crs = rotation_crs, xlim = c(14.25, 19.75), ylim = c(7.38, 12.7)) +
  labs(x = "", y = "")

### Saving the plot
ggsave("Projet_1/Interpolation_finale/Etude_de_cas/Graphs/rmse_tasmin_B1_ked.png", width = 4,
       height = 4, bg = "white")

## MAXIMUM TEMPERATURE
### Breaks for the color scale
min_break <- round(min(rmse_tasmax_spr, rmse_tasmax_ked))
max_break <- round(max(rmse_tasmax_spr, rmse_tasmax_ked))
mean_break <- round(0.5 * (min_break + max_break))

# Data for the map
rcoord.sf <- data_coord[test_indexes, ]
data.graph <- data.frame(rcoord.sf, values_spr = rmse_tasmax_spr, values_ked = rmse_tasmax_ked)

data.graph.sf <- st_as_sf(data.graph, coords = c("lon", "lat"), crs = rotation_crs)

##################################################### LEGEND
p <- ggplot() +
  geom_sf(data = data.graph.sf, aes(color = values_spr), size = 6, shape = 15) +  # Apply color to points
  scale_colour_gradientn(colors = discrete_palette_rmse,
          limits = c(min_break, max_break+0.5),
           breaks = c(min_break, mean_break, max_break),
          guide = guide_colourbar(barheight = 8, barwidth = 0.5, title.position = "right"),
                         name = "tmax (°C)") +
  theme(legend.title = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 15),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "right",
        legend.margin = margin(t = -30),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  geom_sf(data = world, fill = "transparent", color = "black", linewidth = 0.2) +
  coord_sf(crs = rotation_crs, xlim = c(14.25, 19.75), ylim = c(7.38, 12.7)) +
  labs(x = "", y = "")
# Extract legend
legend <- get_legend(p)

ggsave("Projet_1/Interpolation_finale/Etude_de_cas/Graphs/legend_rmse_tasmax_B1.png", width = 1,
       height = 4, bg = "white", plot = cowplot::plot_grid(legend))
#############################################################################################

### SPR
ggplot() +
  geom_sf(data = data.graph.sf, aes(color = values_spr), size = 6, shape = 15) +  # Apply color to points
  scale_colour_gradientn(colors = discrete_palette_rmse,
          limits = c(min_break, max_break+0.5),
           breaks = c(min_break, mean_break, max_break),
          guide = guide_colourbar(barheight = 8, barwidth = 0.5, title.position = "right"),
                         name = "RMSE (°C)") +
  geom_sf(data = points_train_sf, shape = 22, color = "black", fill = "green", size = 4, stroke = 1) +  # train points
  theme(legend.title = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 10),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "none",
        legend.margin = margin(t = -30),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  geom_sf(data = world, fill = "transparent", color = "black", linewidth = 0.2) +
  coord_sf(crs = rotation_crs, xlim = c(14.25, 19.75), ylim = c(7.38, 12.7)) +
  labs(x = "", y = "")

### Saving the plot
ggsave("Projet_1/Interpolation_finale/Etude_de_cas/Graphs/rmse_tasmax_B1_spr.png", width = 4,
       height = 4, bg = "white")

### KED
ggplot() +
  geom_sf(data = data.graph.sf, aes(color = values_ked), size = 6, shape = 15) +  # Apply color to points
  scale_colour_gradientn(colors = discrete_palette_rmse,
          limits = c(min_break, max_break+0.5),
           breaks = c(min_break, mean_break, max_break),
          guide = guide_colourbar(barheight = 8, barwidth = 0.5, title.position = "right"),
                         name = "RMSE (°C)") +
  geom_sf(data = points_train_sf, shape = 22, color = "black", fill = "green", size = 4, stroke = 1) +  # train points
  theme(legend.title = element_text(angle = 90, vjust = 1, hjust = 0.5, size = 15),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.position = "none",
        legend.margin = margin(t = -30),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  geom_sf(data = world, fill = "transparent", color = "black", linewidth = 0.2) +
  coord_sf(crs = rotation_crs, xlim = c(14.25, 19.75), ylim = c(7.38, 12.7)) +
  labs(x = "", y = "")

### Saving the plot
ggsave("Projet_1/Interpolation_finale/Etude_de_cas/Graphs/rmse_tasmax_B1_ked.png", width = 4,
       height = 4, bg = "white")
##########################################################################################

########################### SPEARMAN CORRELATION FOR THE 3 DAYS ##########################
ssim_day <- function(Z, Z_hat, K1 = 0.01, K2 = 0.03, L = NULL) {
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

# PRECIPITATION
## SPR
day_pr <- 3519
month_pr <- day_to_month(day_pr)
climatology_pr_B1 <- climatology_pr[, region]
climatology_pr_B1 <- climatology_pr_B1[, test_indexes]
climatology_pr_B1_day <- climatology_pr_B1[month_pr, ]
corr_spr_pr <- cor(climatology_pr_B1_day, data_pr_interp_spr[day_pr, test_indexes], method = "spearman")
corr_spr_pr
ssim_spr_pr <- ssim_day(data_pr_obs_0[day_pr, test_indexes], data_pr_interp_spr[day_pr, test_indexes])
ssim_spr_pr

## KED
corr_ked_pr <- cor(climatology_pr_B1_day, data_pr_interp_ked[day_pr, test_indexes], method = "spearman")
corr_ked_pr
ssim_ked_pr <- ssim_day(data_pr_obs_0[day_pr, test_indexes], data_pr_interp_ked[day_pr, test_indexes])
ssim_ked_pr

# TASMIN
## SPR
day_tasmin <- 100
month_tasmin <- day_to_month(day_tasmin)
climatology_tasmin_B1 <- climatology_tasmin[, region]
climatology_tasmin_B1 <- climatology_tasmin_B1[, test_indexes]
climatology_tasmin_B1_day <- climatology_tasmin_B1[month_tasmin, ]
corr_spr_tasmin <- cor(climatology_tasmin_B1_day, data_tasmin_interp_spr[day_tasmin, test_indexes], method = "spearman")
corr_spr_tasmin
ssim_spr_tasmin <- ssim_day(data_tasmin_obs_0[day_tasmin, test_indexes], data_tasmin_interp_spr[day_tasmin, test_indexes])
ssim_spr_tasmin

## KED
corr_ked_tasmin <- cor(climatology_tasmin_B1_day, data_tasmin_interp_ked[day_tasmin, test_indexes], method = "spearman")
corr_ked_tasmin
ssim_ked_tasmin <- ssim_day(data_tasmin_obs_0[day_tasmin, test_indexes], data_tasmin_interp_ked[day_tasmin, test_indexes])
ssim_ked_tasmin

# TASMAX
## SPR
day_tasmax <-  100
month_tasmax <- day_to_month(day_tasmax)
climatology_tasmax_B1 <- climatology_tasmax[, region]
climatology_tasmax_B1 <- climatology_tasmax_B1[, test_indexes]
climatology_tasmax_B1_day <- climatology_tasmax_B1[month_tasmax, ]
corr_spr_tasmax <- cor(climatology_tasmax_B1_day, data_tasmax_interp_spr[day_tasmax, test_indexes], method = "spearman")
corr_spr_tasmax
ssim_spr_tasmax <- ssim_day(data_tasmax_obs_0[day_tasmax, test_indexes], data_tasmax_interp_spr[day_tasmax, test_indexes])
ssim_spr_tasmax
## KED
corr_ked_tasmax <- cor(climatology_tasmax_B1_day, data_tasmax_interp_ked[day_tasmax, test_indexes], method = "spearman")
corr_ked_tasmax
ssim_ked_tasmax <- ssim_day(data_tasmax_obs_0[day_tasmax, test_indexes], data_tasmax_interp_ked[day_tasmax, test_indexes])
ssim_ked_tasmax
##########################################################################################
sum(corr_spr_pr, corr_spr_tasmin, corr_spr_tasmax)/3
sum(corr_ked_pr, corr_ked_tasmin, corr_ked_tasmax)/3
sum(ssim_spr_pr, ssim_spr_tasmin, ssim_spr_tasmax)/3
sum(ssim_ked_pr, ssim_ked_tasmin, ssim_ked_tasmax)/3
