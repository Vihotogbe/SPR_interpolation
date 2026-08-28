########################## DATA PROCESSING
# Packages loading
install.packages <- c("bigstatsr", "bigmemory", "gstat", "fields", "RNetCDF", "ggplot2",
                    "ggspatial", "maps", "tmap", "sf", "sp", "rnaturalearth",
                    "rnaturalearthdata", "rnaturalearthhires", "dplyr",
                    "RColorBrewer", "gridExtra", "scales")
new.packages <- install.packages[!(install.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)


library(bigstatsr)
library(bigmemory)
library(gstat)
library(fields)
library(RNetCDF)
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
### Load the necessary packages

# Function to process data
data_processing <- function(dir_seq, var) {
  dir_seq <- dir_seq
  ## Function do process data for 1 year
  data_process <- function(data) {
    ### Re-organizing matrix of pixels into a vector
    data.mat <- matrix(ncol = 280 * 280,
                       nrow = 365) # ndays
    ### Data loading
    data_nc <- data
    ### Data extraction
    data_extract <- var.get.nc(data_nc, var)
    ### Function to fill data.mat et coord matrix for one year
    k <- 1
    for (i in 1:280)  {
      for (j in 1:280) {
        data.mat[, k] <- data_extract[i, j, ]
        k <- k + 1
      }
    }
    return(data.mat)
  }
  ## Function to process data for all years for
  data_all <- NULL
  nbre_elements <- length(dir_seq)
  for (yy in 1:nbre_elements) {
    data_nc <- open.nc(dir_seq[yy])
    results_data <- data_process(data = data_nc)
    data_all <- rbind(data_all, results_data)
  }
  return(data_all)
}

# Period of extraction
years <- seq(1980, 2029)
dir_seq.pr <- paste("/home/julie/Data/Climex/day/kdj/day_pr_historical-r1-r10i1p1-rcp85_kdj",
                    years, "se.nc", sep = "_")
dir_seq.tasmin <- paste("/home/julie/Data/Climex/day/kdj/day_tasmin_historical-r1-r10i1p1-rcp85_kdj",
                        years, "se.nc", sep = "_")
dir_seq.tasmax <- paste("/home/julie/Data/Climex/day/kdj/day_tasmax_historical-r1-r10i1p1-rcp85_kdj",
                        years, "se.nc", sep = "_")

# Extract coordinates
coord <- matrix(ncol = 2, nrow = 280 * 280)
colnames(coord) <- c("lon", "lat")
data_nc <- open.nc(dir_seq.pr[1])
lat <- var.get.nc(data_nc, "lat")[, , 1] # same for 1:365
lon <- var.get.nc(data_nc, "lon")[, , 1] # same for 1:365
k <- 1
for (i in 1:280){
  for (j in 1:280){
    coord[k, ] <- c(lon[i, j], lat[i, j])
    k <- k + 1
  }
}
dim(coord) # 78400 pixels

# save the coordinates as .csv file
# write.csv(coord, "Data/coord.csv", row.names = FALSE)

# Create an sf object from coordinates (lon, lat)
# Create a data frame with the coordinates from rorated pole
# The crs for rotated coordinate system
rotation_crs <- paste0("+proj=ob_tran +o_proj=longlat +o_lon_p=", 0,
                       " +o_lat_p=", 42.5, " +lon_0=", 83-180, " +to_meter=0.01745329")

## Transform data.coord to a spatial object
rcoord <- st_as_sf(as.data.frame(coord), coords = c("lon", "lat"), crs = 4326)
## Convert data.coord to the rotated pole
rcoord.sf <- st_transform(rcoord, crs = rotation_crs)
## Extract the coordinates into a data.frame
rcoord.sf.df <- as.data.frame(st_coordinates(rcoord.sf))
#colnames(rcoord.sf.df) <- c("lon", "lat")


# Climex coverage and the region of interest
rlon <- rcoord.sf.df[, "X"]
rlat <- rcoord.sf.df[, "Y"]

### GLOBAL VISUALIZATION FOR THE SUBREGION
### Rect of the region of interest
### Climex coverage and the region of interest
#### climex coverage
lat_min <- min(rlat)
lat_max <- max(rlat)
lon_min <- min(rlon)
lon_max <- max(rlon)
rect_climex <- data.frame(
  lon = c(lon_min, lon_max, lon_max, lon_min, lon_min),
  lat = c(lat_min, lat_min, lat_max, lat_max, lat_min)
)
#### region of interest
lat_min <- min(rlat) + 7
lat_max <- max(rlat) - 2
lon_min <- min(rlon) + 5
lon_max <- max(rlon) - 8
rect_QC_region <- data.frame(
  lon = c(lon_min, lon_max, lon_max, lon_min, lon_min),
  lat = c(lat_min, lat_min, lat_max, lat_max, lat_min)
)
### Tranforming to match the crs of the rotated coordinate system
rect_climex_sf <- st_as_sf(rect_climex, coords = c("lon", "lat"), crs = rotation_crs) %>%
    st_combine() %>% st_cast("POLYGON")
rect_QC_region_sf <- st_as_sf(rect_QC_region, coords = c("lon", "lat"), crs = rotation_crs) %>%
    st_combine() %>% st_cast("POLYGON")


# Plot of the climex and interest region
# Charge the world rivers and lakes data
if (!file.exists("maps/ne_10m_rivers.shp")) {
  ne_download(scale = 10, type = "rivers_lake_centerlines", category = "physical",
              destdir = "maps/", load = FALSE) # major rivers
  ne_download(scale = 10, type = "lakes", category = "physical",
              destdir = "maps/", load = FALSE) # major lakes
}
world <- ne_countries(scale = 50, returnclass = "sf")
rivers <- ne_load(scale = 10, type = "rivers_lake_centerlines", destdir = "maps", returnclass = "sf")
lakes <- ne_load(scale = 10, type = "lakes", destdir = "maps", returnclass = "sf")

# Plot the rotated coordinate system in the entire region
ggplot() +
geom_sf(data = world, fill = "transparent") +
geom_sf(data = rivers, colour = "blue", linewidth = 0.2) +
geom_sf(data = lakes, fill = "lightblue") +
geom_sf(data = rect_climex_sf, color = "darkred", linewidth = 1, fill = "transparent") +
# geom_sf(data = rect_QC_region_sf, color = "black", linewidth = 1, fill = "transparent") +
coord_sf(crs = st_crs(rcoord.sf), xlim = c(-5, 35), ylim = c(-15, 25)) +
scale_x_continuous(breaks = seq(-140, -20, by = 10)) +
scale_y_continuous(breaks = seq(35, 70, by = 5)) +
theme(axis.text = element_text(size = 10)) +
labs(x = "", y = "")
# Saving the plot
# ggsave("Projet_1/Graphs/world_climex_QC_plot_2_1.png", width = 4.5, height = 4, dpi = 500)
ggsave("/home/vihoua@labos.polymtl.ca/Projet_2/Graphs/world_climex_QC_plot_2_1.png", width = 4.5, height = 4, dpi = 500)

# Plot the rotated coordinate system in the region of interest
ggplot() +
geom_sf(data = world, fill = "transparent") +
geom_sf(data = rivers, colour = "blue", linewidth = 0.2) +
geom_sf(data = lakes, fill = "lightblue") +
geom_sf(data = rect_climex_sf, color = "red", linewidth = 1, fill = "transparent") +
geom_sf(data = rect_QC_region_sf, color = "green", linewidth = 1, fill = "transparent") +
coord_sf(crs = st_crs(rcoord.sf), xlim = c(8, 25), ylim = c(-5, 15.5)) +
scale_x_continuous(breaks = seq(-140, -20, by = 5)) +
scale_y_continuous(breaks = seq(40, 60, by = 5)) +
theme(axis.text = element_text(size = 10)) +
labs(x = "", y = "")
# Saving the plot
ggsave("/home/vihoua@labos.polymtl.ca/Projet_2/Graphs/world_climex_QC_zum_plot.png", width = 4.5, height = 4, dpi = 500)



# Extact data for precipitation and temperature
# PRECIPITATIONS
data.pr <- data_processing(dir_seq = dir_seq.pr, var = "pr")
### Different transformations for PRECIPITATION DATA
data.pr <- data.pr * 86400 # Transformer pr en mm/day (from kg/m2/s)

# MINIMUM TEMPERATURES
data.tasmin <- data_processing(dir_seq = dir_seq.tasmin, var = "tasmin")
# Transforming MINIMUM TEMPERATURE DATA to celsius (from kevin)
data.tasmin <- data.tasmin - 273.15

# MAXIMUM TEMPERATURES
data.tasmax <- data_processing(dir_seq = dir_seq.tasmax, var = "tasmax")
### Transforming MAXIMUM TEMPERATURE DATA to celsius (from kevin)
data.tasmax <- data.tasmax - 273.15


# Select Quebec subregion with rotated coordinates and corresponding data
quebec_region <- which(rlat >= -5.615007 & rlat <= 16.075004 & rlon >= 7.694992 & rlon <= 25.385023,
                       arr.ind = TRUE)
data.coord <- coord[quebec_region, ]
data.coord.sf <- rcoord.sf.df[quebec_region, ]
colnames(data.coord.sf) <- c("lon", "lat")
data.pr <- data.pr[, quebec_region]
data.tasmin <- data.tasmin[, quebec_region]
data.tasmax <- data.tasmax[, quebec_region]

# Save the data
save(coord, file = "Data/coord.RData")
save(data.coord, file = "Data/data.coord.RData")
save(data.coord.sf, file = "Data/data.coord.sf.RData")
save(data.pr, file = "Data/data.pr.RData")
save(data.tasmin, file = "Data/data.tasmin.RData")
save(data.tasmax, file = "Data/data.tasmax.RData")
#########################################################


# Defining the period of RCM and observations data
rcm_period <- 1:10950 # 30 years of data (1980-2009)
obs_period_1 <- 7300:10950 # 10 years of data (2000-2009)
obs_period_2 <- 10951:14600 # 10 years of data (2010-2019)
obs_period_3 <- 14601:18250 # 10 years of data (2020-2029)


## Function to calculate the climatology

### Function to calculate mensual climatology for temperature data
climatology_mensual_temp <- function(data, nyears = 30) {
  # Function to calculate monthly means
  monthly_means <- function(data) {
    # Define the days corresponding to each month
    months <- list(
      January = 1:31,
      February = 32:59,
      March = 60:90,
      April = 91:120,
      May = 121:151,
      June = 152:181,
      July = 182:212,
      August = 213:243,
      September = 244:273,
      October = 274:304,
      November = 305:334,
      December = 335:365
    )
    # Initialize a matrix to store monthly means
    monthly_means <- matrix(0, nrow = 12, ncol = ncol(data))
    rownames(monthly_means) <- names(months)
    # Calculate the mean for each month
    for (i in 1:12) {
      month_days <- months[[i]]
      monthly_means[i, ] <- colMeans(data[month_days, ])
    }
    return(monthly_means)
  }
  # Create a matrix with 12 rows to store annual monthly means
  ndays <- 365
  days <- seq(1:ndays)
  moyennes_months <- matrix(0, nrow = 12, ncol = ncol(data))
  # Loop to calculate the monthly mean for each year and sum them over years
  for (i in 1:nyears) {
    data_year <- data[days, ]
    moyennes_months_year <- monthly_means(data_year)
    moyennes_months <- moyennes_months + moyennes_months_year
    days <- days + ndays
  }
  # Calculate the climatology by dividing the sum of monthly means over years
  # by the number of years
  climatology_months <- moyennes_months / nyears
  return(climatology_months)
}

### Mensual climatology for tasmin data
data.tasmin.rcm <- data.tasmin[rcm_period, ]
climatology_tasmin <- climatology_mensual_temp(data = data.tasmin.rcm, nyears = 30)


### Mensual climatology for tasmax data
data.tasmax.rcm <- data.tasmax[rcm_period, ]
climatology_tasmax <- climatology_mensual_temp(data = data.tasmax.rcm, nyears = 30)


### Function to calculate mensual climatology for precipitation data
climatology_mensual_pr <- function(data, nyears = 30) {
  # Function to calculate monthly total
  monthly_total <- function(data) {
    # Define the days corresponding to each month
    months <- list(
      January = 1:31,
      February = 32:59,
      March = 60:90,
      April = 91:120,
      May = 121:151,
      June = 152:181,
      July = 182:212,
      August = 213:243,
      September = 244:273,
      October = 274:304,
      November = 305:334,
      December = 335:365
    )
    # Initialize a matrix to store monthly means
    monthly_sum <- matrix(0, nrow = 12, ncol = ncol(data))
    rownames(monthly_sum) <- names(months)
    # Calculate the sum for each month
    for (i in 1:12) {
      month_days <- months[[i]]
      monthly_sum[i, ] <- colSums(data[month_days, ])
    }
    return(monthly_sum)
  }
  # Create a matrix with 12 rows to store annual monthly means
  ndays <- 365
  days <- seq(1:ndays)
  total_months <- matrix(0, nrow = 12, ncol = ncol(data))
  # Loop to calculate the monthly total for each year and sum them over years
  for (i in 1:nyears) {
    data_year <- data[days, ]
    total_months_year <- monthly_total(data_year)
    total_months <- total_months + total_months_year
    days <- days + ndays
  }
  # Calculate the climatology by dividing the sum of monthly means over years
  # by the number of years
  climatology_months <- total_months / nyears
  return(climatology_months)
}


data.pr.rcm <- data.pr[rcm_period, ]
climatology_pr <- climatology_mensual_pr(data = data.pr.rcm, nyears = 30)


# Save the climatology
save(climatology_tasmin, file = "Data/climatology_tasmin.RData")
save(climatology_tasmax, file = "Data/climatology_tasmax.RData")
save(climatology_pr, file = "Data/climatology_pr.RData")



# REGIONS A and B DETERMINATION : Showing the 3 sizes of region
# Delimitation of the two regions of comparison
## Region A : around Montreal (south QC)
lat_min.1 <- -2
lat_max.1 <- 4
lon_min.1 <- 13
lon_max.1 <- 19
rect_A1 <- data.frame(
  lon = c(lon_min.1, lon_max.1, lon_max.1, lon_min.1, lon_min.1),
  lat = c(lat_min.1, lat_min.1, lat_max.1, lat_max.1, lat_min.1)
)

rect_A1_sf <- st_as_sf(rect_A1, coords = c("lon", "lat"), crs = rotation_crs) %>% 
    st_combine() %>% st_cast("POLYGON")

sz <- 1
rect_A2 <- data.frame(
  lon = c(lon_min.1 + sz, lon_max.1 - sz, lon_max.1 - sz, lon_min.1 + sz, lon_min.1 + sz),
  lat = c(lat_min.1 + sz, lat_min.1 + sz, lat_max.1 - sz, lat_max.1 - sz, lat_min.1 + sz)
)
rect_A2_sf <- st_as_sf(rect_A2, coords = c("lon", "lat"), crs = rotation_crs) %>%
                        st_combine() %>% st_cast("POLYGON")

sz <- 2
rect_A3 <- data.frame(
  lon = c(lon_min.1 + sz, lon_max.1 - sz, lon_max.1 - sz, lon_min.1 + sz, lon_min.1 + sz),
  lat = c(lat_min.1 + sz, lat_min.1 + sz, lat_max.1 - sz, lat_max.1 - sz, lat_min.1 + sz)
)
rect_A3_sf <- st_as_sf(rect_A3, coords = c("lon", "lat"), crs = rotation_crs) %>%
                        st_combine() %>% st_cast("POLYGON")

## Region B : around north QC
lat_min.2 <- 7
lat_max.2 <- 13
lon_min.2 <- 14
lon_max.2 <- 20
rect_B1 <- data.frame(
  lon = c(lon_min.2, lon_max.2, lon_max.2, lon_min.2, lon_min.2),
  lat = c(lat_min.2, lat_min.2, lat_max.2, lat_max.2, lat_min.2)
)
## Tranforming to match the crs of the rotated coordinate system
rect_B1_sf <- st_as_sf(rect_B1, coords = c("lon", "lat"), crs = rotation_crs) %>%
    st_combine() %>% st_cast("POLYGON")

sz <- 1
rect_B2 <- data.frame(
  lon = c(lon_min.2 + sz, lon_max.2 - sz, lon_max.2 - sz, lon_min.2 + sz, lon_min.2 + sz),
  lat = c(lat_min.2 + sz, lat_min.2 + sz, lat_max.2 - sz, lat_max.2 - sz, lat_min.2 + sz)
)
rect_B2_sf <- st_as_sf(rect_B2, coords = c("lon", "lat"), crs = rotation_crs) %>%
                        st_combine() %>% st_cast("POLYGON")

sz <- 2
rect_B3 <- data.frame(
  lon = c(lon_min.2 + sz, lon_max.2 - sz, lon_max.2 - sz, lon_min.2 + sz, lon_min.2 + sz),
  lat = c(lat_min.2 + sz, lat_min.2 + sz, lat_max.2 - sz, lat_max.2 - sz, lat_min.2 + sz)
)
rect_B3_sf <- st_as_sf(rect_B3, coords = c("lon", "lat"), crs = rotation_crs) %>%
                        st_combine() %>% st_cast("POLYGON")


## Plot the regions with CRS that preserves right angles
ggplot() +
  geom_sf(data = world, fill = "transparent", color = "black", linewidth = 0.2) +
  geom_sf(data = rect_A1_sf, color = "blue", linewidth = 1, fill = "transparent") +
  # geom_sf(data = rect_A2_sf, color = "red", linewidth = 1, fill = "transparent") +
  # geom_sf(data = rect_A3_sf, color = "green", linewidth = 1, fill = "transparent") +
  # geom_sf(data = rect_B1_sf, color = "blue", linewidth = 1, fill = "transparent") +
  # geom_sf(data = rect_B2_sf, color = "red", linewidth = 1, fill = "transparent") +
  # geom_sf(data = rect_B3_sf, color = "green", linewidth = 1, fill = "transparent") +
  coord_sf(crs = rotation_crs, xlim = c(8, 25), ylim = c(-5, 15.5)) +
  scale_x_continuous(breaks = seq(-140, -20, by = 5)) +
  scale_y_continuous(breaks = seq(40, 60, by = 5)) +
  theme(axis.text = element_text(size = 10)) +
  labs(x = "", y = "")
## Saving the plot
# ggsave("Projet_1/Graphs/plot_2_regions.png", width = 4.5, height = 4)
ggsave("/home/vihoua@labos.polymtl.ca/Projet_2/Graphs/plot_2_regions.png", width = 4.5, height = 4, dpi = 500)

# Regions of interest
region_A1 <- which(rlat >= -2 & rlat <= 4 & rlon >= 13 & rlon <= 19, arr.ind = TRUE)
region_A2 <- which(rlat >= -1 & rlat <= 3 & rlon >= 14 & rlon <= 18, arr.ind = TRUE)
region_A3 <- which(rlat >= 0 & rlat <= 2 & rlon >= 15 & rlon <= 17, arr.ind = TRUE)
region_B1 <- which(rlat >= 7.5 & rlat <= 13.5 & rlon >= 14 & rlon <= 20, arr.ind = TRUE)
region_B2 <- which(rlat >= 8 & rlat <= 12 & rlon >= 15 & rlon <= 19, arr.ind = TRUE)
region_B3 <- which(rlat >= 9 & rlat <= 11 & rlon >= 16 & rlon <= 18, arr.ind = TRUE)

length(region_B1)
