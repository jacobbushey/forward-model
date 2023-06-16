library(tidyverse)
library(ggmap)
library(ncdf4)
library(terra)
library(viridis)
library(ggplot2)

library(pracma)
library(imager)
#library(magick)
library(edges)

setwd("~/forward-model/code")
source("ssec.R")  

setwd("~/forward-model/data/Mosaic")

file <- 'MAIR_2021_08_11_RF08_gim-delivery_e9f5868_2023-06-01_mosaic_30m_MethaneAIR_L3_mosaic_20210811T174014_20210811T193548_dpp.nc'

nc_data <- nc_open(file)
# Save the print(nc) dump to a text file
{
  sink('RF08_file_contents.txt')
  print(nc_data)
  sink()
}

mair.lon <- ncvar_get(nc_data, "lon")
mair.lat <- ncvar_get(nc_data, "lat")
mair.p0 <- ncvar_get(nc_data, "apriori_data/surface_pressure")

mair.tau_start <- ncatt_get(nc_data, varid = 0, attname = "time_coverage_start")$value
mair.tau_end <- ncatt_get(nc_data, varid = 0, attname = "time_coverage_end")$value

mair.xch4 <- ncvar_get(nc_data, "xch4")
mair.xch4[mair.xch4 > (10^35)] = NA
#mair.xch4 <- mair.xch4[4536:1, 1:3492]
mair.xch4 <- t(mair.xch4)
mair.xch4 <- mair.xch4[3492:1, 1:4536]

#extent <- ext(c(min(mair.lat), max(mair.lat), min(mair.lon), max(mair.lon)))
extent <- ext(c(min(mair.lon), max(mair.lon), min(mair.lat), max(mair.lat)))
segment_rast <- rast(mair.xch4, extent = extent)
plot(segment_rast)

# MAKE A MASK A THAT HAS A VALUE OF 1 IF THERE ARE POINTS IN THE GRID CELL, AND A VALUE OF 0 IF THERE ARE NOT. IT ALSO SMOOTHS OUT THE DATA POINTS. CREATING A BUFFER REGION SO THAT YOU DON'T MISS ANYTHING
#se <- strel("disk",200) #structuring element used to dilate the image
#A <- imdilate(edge(isnan(mair.xch4), "roberts"), se);
#A <- imager::dilate(edge(isnan(mair.xch4), "roberts"), se) #or is the package called erode?
#mair.xch4.df <- as.data.frame(mair.xch4)
#mair.xch4.t <- t(mair.xch4)[3492:1, 1:4536]
mair.xch4.cimg <- imager::as.cimg(!is.na(mair.xch4[3492:1, 1:4536]))
#temp <- imager::cannyEdges(mair.xch4.cimg)
#A <- imager::dilate(magick::edges())
A <- imager::isoblur(mair.xch4.cimg, sigma = 10)
A <- imager::imrotate(A, -90)
A.mat <- t(as.matrix(A))
A.rast <- rast(A.mat, extent = extent) #DOING A LOT OF TRANSPOSING HERE... IT SEEMS LIKE EACH OF THESE TRANSFORMATIONS EITHER TRANSPOSES OR FLIPS MY IAMGE
A.df <- as.data.frame(A.rast, xy = TRUE)
colnames(A.df) <- c('lon', 'lat', 'xch4')

segment_df <- as.data.frame(segment_rast, xy = TRUE)
colnames(segment_df) <- c('lon', 'lat', 'xch4')

colorbar_min <- 1800
colorbar_max <- 1930

ggplot(data = segment_df) +
  geom_raster(mapping = aes(x = lon, y = lat, fill = xch4)) +
  scale_fill_gradientn(colours = ssec(100), limits = c(colorbar_min, colorbar_max))

ggplot(data = segment_df) +
  geom_histogram(mapping = aes(x = xch4), binwidth = 5, colour = 'blue', fill = 'white') +
  xlim(c(colorbar_min, colorbar_max)) +
  ggtitle(file) +
  theme(plot.title = element_text(hjust = 0.5))

## Setting up the limits for the grid -----------------------

lat <- array()
lat[1] <- min(mair.lat) - 0.05
lat[2] <- max(mair.lat) + 0.05

lon <- array()
lon[1] <- min(mair.lon) - 0.05
lon[2] <- max(mair.lon) + 0.05

res <- 17

dx <- (lon[2] - lon[1])/res
dy <- (lat[2] - lat[1])/res

# generate grid ----------------
x <- pracma::linspace(lon[1], lon[2], res) #NEED TO ISNTALL PACKAGE PRACMA
y <- pracma::linspace(lat[1], lat[2], res)

# select points inside the mask
X <- matrix()
Y <- matrix()
tau <- matrix()
N <- matrix()

#for (i in c(1, res)){
#  for (j in c(1, res)){
#    cx <- min(abs(mair.lon-x[i]));
#    cy <- min(abs(mair.lat-y[j]));
#    if (A[cx,cy] > 0){
#      X <- c(X, x[i]);
#      Y <- c(Y, y[j]);
#      #tau <- c(tau mair.tau(cx,cy));
#    end
#    }
#  }
#}

for (i in c(1, res-1)){
  for (j in c(1, res-1)){
    mair.xch4.inbox <- mair.xch4 %>%
      dplyr::mutate()
    
    if (A[cx,cy] > 0){
      X <- c(X, x[i]);
      Y <- c(Y, y[j]);
      #tau <- c(tau mair.tau(cx,cy));
      end
    }
  }
}