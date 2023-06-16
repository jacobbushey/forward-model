library(tidyverse)
library(ggmap)
library(ncdf4)
library(terra)
library(viridis)
library(ggplot2)

library(pracma)
library(imager)
library(lubridate)

setwd("~/forward-model/data/Mosaic")

# Load the RF08 data file
file <- 'MAIR_2021_08_11_RF08_gim-delivery_e9f5868_2023-06-01_mosaic_30m_MethaneAIR_L3_mosaic_20210811T174014_20210811T193548_dpp.nc'

# Open the data file and save the contents to a text file
nc_data <- nc_open(file)
# Save the print(nc) dump to a text file
{
  sink('RF08_file_contents.txt')
  print(nc_data)
  sink()
}

# Retrieve the variables of interest
# JULIAN USED THIS DOT NOTATION TO BUILD A STRUCT, SHOULD I BE DOING SOMETHING SIMILAR, INSTEAD OF VARIABLES?
# LIKELY DOESN'T MATTER, AS LONG AS I GET A GOOD RECEPTOR LIST
mair.lon <- ncvar_get(nc_data, "lon")
mair.lat <- ncvar_get(nc_data, "lat")
mair.p0 <- ncvar_get(nc_data, "apriori_data/surface_pressure")
mair.tau.start <- lubridate::as_datetime(ncatt_get(nc_data, varid = 0, attname = "time_coverage_start")$value, tz = 'UTC') #DO I HAVE THE RIGHT TIMEZONE?
mair.tau.end <- lubridate::as_datetime(ncatt_get(nc_data, varid = 0, attname = "time_coverage_end")$value, tz = 'UTC') # DO I HAVE THE RIGHT TIMEZONE?

# Retrieve the xch4 data and reorient it to match the proper lat/lon
mair.xch4.mat <- ncvar_get(nc_data, "xch4")
mair.xch4.mat[mair.xch4.mat > (10^35)] = NA
mair.xch4.mat <- t(mair.xch4.mat)
mair.xch4.mat <- mair.xch4.mat[3492:1, 1:4536]

# Create a raster object and plot to check the lat/lon
extent <- ext(c(min(mair.lon), max(mair.lon), min(mair.lat), max(mair.lat)))
mair.xch4.rast <- rast(mair.xch4.mat, extent = extent)
plot(mair.xch4.rast)

# Convert the matrix to a cimg so it can be blurred, creating a mask of regions where data exists
# the .cimg will appear rotate, but note the strange axes. What matters is that the dimensions match with the original.
mair.xch4.cimg <- imager::as.cimg(!is.na(mair.xch4.mat[3492:1, 1:4536]))
# blur the image to create a 'buffer' region, extending the range of your region of interest beyond where data was collected.
# this is necessary as the flight pattern is not aligned with the mosaic (personal correspondance with Julian)
A <- imager::isoblur(mair.xch4.cimg, sigma = 200) # change sigma to change how much the image is blurred, changing your region of interest
x11() # plot A to see how much the image has been blurred
plot(A)

# set up the limits for the grid
lat <- array()
lat[1] <- min(mair.lat) - 0.05
lat[2] <- max(mair.lat) + 0.05

lon <- array()
lon[1] <- min(mair.lon) - 0.05
lon[2] <- max(mair.lon) + 0.05

# set the resolution of the grid. THIS IS THE NUMBER THAT JULIAN USED
res <- 17

# set the step size for the space between your grid points 
dx <- (lon[2] - lon[1])/res
dy <- (lat[2] - lat[1])/res

# generate the grid
x <- pracma::linspace(lon[1], lon[2], res) #NEED TO ISNTALL PACKAGE PRACMA
y <- pracma::linspace(lat[1], lat[2], res)

## Retrieve time information from each individual segment ------------------------------------

setwd('~/forward-model/data/RF08')
files <- list.files(path = '~/forward-model/data/RF08', pattern = NULL)

tick <- 1
for (name in files){
  nc_data <- nc_open(name)
  varnames <- list(paste0('segment', tick, '.lon'),
                   paste0('segment', tick, '.lat'),
                   paste0('segment', tick, '.xch4.mat'),
                   paste0('segment', tick, '.tau.start'),
                   paste0('segment', tick, 'tau.end'))
  assign(as.character(varnames[1]), ncvar_get(nc_data, "lon"))
  assign(as.character(varnames[2]), ncvar_get(nc_data, "lat"))
  assign(as.character(varnames[3]), ncvar_get(nc_data, "xch4"))
  assign(as.character(varnames[4]), as_datetime(ncatt_get(nc_data, varid = 0, attname = "time_coverage_start")$value, tz = 'UTC'))
  assign(as.character(varnames[5]), as_datetime(ncatt_get(nc_data, varid = 0, attname = "time_coverage_end")$value, tz = 'UTC'))
  tick <- tick + 1
  
  # CREATE A LIST THAT KEEPS TRACK OF THESE VARIABLE NAMES SO I CAN GO THROUGH THEM ALL ITERATIVELY LATER
  time.int <- interval(start = mair.tau.start, end = mair.tau.end)
  mid.time <- mair.tau.start + (time.int / 2)
  
  # ITERATE THROUGH THIS LIST OF FILES, AND IF THE LAT AND LON OF THE BIG DATASET ARE FOUND WITHIN A SPECIFIC
  # SEGMENT, THEN ASSIGN THAT TIME VALUE
  # WHAT TO DO ABOUT OVERLAPS?
}


## -------------------------------------------------------------------------------------------

# select points inside the mask
# make arrays to store their x and y coordinates, timestamp, and the point ID number
X <- array()
Y <- array()
tau <- array()
N <- array()

tick <- 1 #counter to move to the next index in the arrays
for (i in c(1:res)){
  print(i) #printing to make sure that the loop is progressing properly
 for (j in c(1:res)){
   print(j)
   cx <- match(min(abs(mair.lon-x[i])), 
               abs(mair.lon-x[i]))
   cy <- match(min(abs(mair.lat-y[j])), 
               abs(mair.lat-y[j]))
   # finding the lon and lat indeces of the nearest data point to the current grid point by minimizing the distance between the two
   if (A[cy, cx, 1, 1] > 0.01){
     # this conditional statement means that the grid cell falls within our blurred mask. Note cy and cx look flipped because cy is latitude, which
     # will be in the row index, while cx is longitude, which will be in the column index
     X[tick] <- x[i] # if there is data near that grid point, then add the gridpoint to the receptor list
     Y[tick] <- y[j]
     N[tick] <- tick
     tick <- tick + 1
     #tau <- c(tau mair.tau(cx,cy)) # THIS IS WHERE I NEED TO ADD THE TIMESTAMP, BUT UNSURE HOW TO DO THIS
     # Either need to come up with some way to step from the tau max to tau min given in the file, or (more likely) use the timestamps from the individual legs that Josh gave me
   end
    }
  }
}

# Create a data frame holding the coordinates of the receptors
receptor.df <- data.frame(X, Y, N)

## Plotting the data ---------------------------------------------------
# Load the ssec color scheme
setwd("~/forward-model/code")
source("ssec.R") 

setwd("~/forward-model/code")
load('rf08_grid.RData')

# Create a dataframe in x, y, z format specifically for plotting
plot.df <- as.data.frame(mair.xch4.rast, xy = TRUE)
colnames(plot.df) <- c('lon', 'lat', 'xch4')

# Plot a histogram of the original data to help inform the colorbar
x11()
ggplot(data = plot.df) +
  geom_histogram(mapping = aes(x = xch4), binwidth = 5, colour = 'blue', fill = 'white') +
  xlim(c(colorbar_min, colorbar_max)) +
  ggtitle(file) +
  theme(plot.title = element_text(hjust = 0.5))

# Set the colorbar for the ssec color scale
colorbar_min <- 1800
colorbar_max <- 1930

# Plot the original data
x11()
ggplot(data = plot.df) +
  geom_raster(mapping = aes(x = lon, y = lat, fill = xch4)) +
  scale_fill_gradientn(colours = ssec(100), limits = c(colorbar_min, colorbar_max))

# Plot the original data with the receptor grid superimposed
x11()
ggplot() +
  geom_raster(data = plot.df, mapping = aes(x = lon, y = lat, fill = xch4)) +
  scale_fill_gradientn(colours = ssec(100), limits = c(colorbar_min, colorbar_max)) +
  geom_point(data = receptor.df, mapping = aes(x = X, y = Y)) +
  geom_text(data = receptor.df, mapping = aes(x = X, y = Y, label = N), vjust = -1)

# ONCE VARIABLES/STRUCTURES ARE FINALIZED, SAVE WORKSPACE SO IT CAN BE EASILY LOADED IN THE FUTURE

# NEED TO CHECK TO MAKE SURE THAT LAT AND LON WERE ASSIGNED CORRECTLY, AND THAT THINGS DIDN'T SHIFT AMIDST THE IMAGE PROCESSING AND OTHER TRANSFORMATIONS





