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

mair.time.int <- interval(start = mair.tau.start, end = mair.tau.end)
mair.time.length <- lubridate::time_length(time.int, unit = "second")
mair.mid.time <- mair.tau.start + (time.length / 2)

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
# NEED TO DISCUSS WHETHER THIS IS AN APPROPRIATE EXTENT OF BLURRING AS WELL AS METHOD

x11() # plot A to see how much the image has been blurred
plot(A)
ggsave(filename = "~/forward-model/figs/rf08_blurred.png", device = png, width = 8, height = 8, units = "in")

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










## Select gridded points inside the mask to be receptors for STILT -------------------------------------------------------------------------------------------

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
               abs(mair.lat-y[j])) # THIS ISN'T EUCLIDEAN DISTANCE, IT'S THE SMALLEST DISTANCE IN TWO DIFFERENT DIMENSIONS. WHAT DID JULIAN DO?
   # COULD USE STATS::DIST INSTEAD
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












## Retrieve time information from each individual segment, used to assign times to the receptors ------------------------------------

# make a list that has the names of the data files for each segment of the flight
setwd('~/forward-model/data/RF08')
files <- list.files(path = '~/forward-model/data/RF08', pattern="(_local.nc)$")

# for loop that assigns values to variables for each segment.
# also finds the 'middle' time for each segment, which will be used when calculating the back trajectories.
# each segment was only about 5 mins long, so this is an appropriate approximation
tick <- 1
#segment.list <- array()
for (name in files){
  nc_data <- nc_open(name)
  varnames <- list(paste0('segment', tick, '.lon'),
                   paste0('segment', tick, '.lat'),
                   paste0('segment', tick, '.xch4.mat'),
                   paste0('segment', tick, '.tau.start'),
                   paste0('segment', tick, '.tau.end'),
                   paste0('segment', tick, '.mid.time'),
                   paste0('segment', tick, '.reps'))
  assign(as.character(varnames[1]), ncvar_get(nc_data, "lon"))
  assign(as.character(varnames[2]), ncvar_get(nc_data, "lat"))
  assign(as.character(varnames[3]), ncvar_get(nc_data, "xch4"))
  assign(as.character(varnames[4]), as_datetime(ncatt_get(nc_data, varid = 0, attname = "time_coverage_start")$value, tz = 'UTC'))
  assign(as.character(varnames[5]), as_datetime(ncatt_get(nc_data, varid = 0, attname = "time_coverage_end")$value, tz = 'UTC'))
  
  # find the time interval, use it to calculate the mid time
  time.int <- interval(start = eval(parse(text = varnames[4])), end = eval(parse(text = varnames[5])))
  time.length <- lubridate::time_length(time.int, unit = "second")
  mid.time <- eval(parse(text = varnames[4])) + (time.length / 2)
  assign(as.character(varnames[6]), mid.time)
  
  # Apply the same transformations to the .xch4.mat for each of the segment matrices that you did for the total matrix.
  # this is necessary in order to get the lat and lon for each segment in dataframe form
  temp.mat <- eval(parse(text = varnames[3]))
  temp.mat[temp.mat > (10^35)] = NA
  temp.mat <- t(temp.mat)
  temp.mat <- temp.mat[dim(temp.mat)[1]:1, 1:dim(temp.mat)[2]]
  extent <- ext(min(eval(parse(text = varnames[1]))), max(eval(parse(text = varnames[1]))), min(eval(parse(text = varnames[2]))), max(eval(parse(text = varnames[2]))))
  temp.rast <- rast(temp.mat, extent = extent)
  temp.df <- terra::as.data.frame(temp.rast, xy = TRUE)
  
  # randomly sample 100 points from each segment as 'representatives' to use for clustering.
  reps.temp <- list()
  sample.idx <- sample(length(temp.df$x), size = 100, replace = FALSE, prob = NULL)
  reps.temp[[1]] <- temp.df$x[sample.idx]
  reps.temp[[2]] <- temp.df$y[sample.idx]
  reps.temp[[3]] <- rep(tick, 100)
  assign(as.character(varnames[7]), reps.temp)
  
  # create a list of segment names, potentially for looping through later
  #segment.list[tick] <- paste0('segment', tick)
  tick <- tick + 1
}

# create a list holding all representative coordinates and their ID's. Check by plotting
all.reps <- vector(mode = 'list', length = 3)
for (i in c(1:23)){
  temp.data <- eval(parse(text = paste0('segment', i, '.reps')))
  all.reps[[1]] <- c(all.reps[[1]], temp.data[[1]])
  all.reps[[2]] <- c(all.reps[[2]], temp.data[[2]])
  all.reps[[3]] <- c(all.reps[[3]], temp.data[[3]])
}

# convert the list to a dataframe and name the columns
all.reps.df <- data.frame(all.reps[[1]], all.reps[[2]], all.reps[[3]])
colnames(all.reps.df) <- c('lon', 'lat', 'id')

# plot to check that the ID's are distributed amongst the representatives appropriately
x11(width = 8, height = 8)
ggplot() +
  geom_point(data = all.reps.df, mapping = aes(x = lon, y = lat, color = as.character(id))) +
  scale_color_manual(values = rainbow(23))
#ggsave(filename = "~/forward-model/figs/rf08_rainbow.png", device = png, width = 8, height = 8, units = "in")

# x11()
# all.reps.df %>%
#   dplyr::filter(id < 6) %>%
#   ggplot() +
#   geom_point(mapping = aes(x = lon, y = lat, color = as.character(id))) +
#   scale_color_manual(values = rainbow(5))

# loop through every receptor. For each receptor, find its nearest neighbors from among the cluster representatives.
# store the index of the representative in receptor.id
# assign the id of the representative to the receptor
receptor.id <- array() # initialize the array for data storage
for (i in c(1:length(receptor.df$X))) {
  temp.list <- vector(mode = 'list', length = 2)
  temp.list[[1]] <- c(receptor.df$X[i], all.reps.df$lon)
  temp.list[[2]] <- c(receptor.df$Y[i], all.reps.df$lat)
  temp.df <- data.frame(temp.list[[1]], temp.list[[2]])
  temp.mat <- as.matrix(temp.df) # construct a matrix where the first row is the receptor lat/lon, and every other row contains a representative lat/lon
  
  # construct a distance matrix to find distance between each set of points
  dist.mat <- as.matrix(dist(temp.df, method = 'euclidean', upper = FALSE))
  
  dist.mat[1, ] <- 1000 # set the entire first row to 1000 so it will never be mistaken as the smallest distance
  min.idx <- which.min(dist.mat[ , 1]) - 1 # find the index of the smallest distance (subtract 1 b/c the first row was for the receptor itself, all set to 1000)
  
  receptor.id[i] <- all.reps.df$id[min.idx] # store this index
  
}

# add the receptor.id list to receptor.df in order to plot the receptors with their ID's and pair them with timestamps for STILT
receptor.df <- receptor.df %>%
  dplyr::mutate(id = receptor.id)

# Assign times according to the id
id.time <- array()
for (i in c(1:length(receptor.df$X))) {
  id <- receptor.df$id[i]
  id.time[i] <- eval(parse(text = paste0('segment', id, '.mid.time')))
  
}
id.time <- as.POSIXct(id.time, tz = 'UTC')

# add times to receptor.df
receptor.df <- receptor.df %>%
  dplyr::mutate(id.time = id.time)










## Plotting the data ---------------------------------------------------
# Load the ssec color scheme
setwd("~/forward-model/code")
source("ssec.R") 

setwd("~/forward-model/data/RF08")
#load('rf08_grid.RData')
load('rf08_grid_2.RData')

# Create a dataframe in x, y, z format specifically for plotting
plot.df <- as.data.frame(mair.xch4.rast, xy = TRUE)
colnames(plot.df) <- c('lon', 'lat', 'xch4')

# Plot a histogram of the original data to help inform the colorbar
x11()
ggplot(data = plot.df) +
  geom_histogram(mapping = aes(x = xch4), binwidth = 5, colour = 'blue', fill = 'white') +
  xlim(c(colorbar_min, 2000)) +
  ggtitle(file) +
  theme(plot.title = element_text(hjust = 0.5))
#ggsave(filename = "~/forward-model/figs/rf08_histogram.png", device = png, width = 8, height = 8, units = "in")

# Set the colorbar for the ssec color scale
colorbar_min <- 1800
colorbar_max <- 1930

# Plot the original data
x11()
ggplot(data = plot.df) +
  geom_raster(mapping = aes(x = lon, y = lat, fill = xch4)) +
  scale_fill_gradientn(colours = ssec(100), limits = c(colorbar_min, colorbar_max))
#ggsave(filename = "~/forward-model/figs/rf08_mosaic.png", device = png, width = 8, height = 8, units = "in")

# Plot the original data with the receptor grid superimposed
x11()
ggplot() +
  geom_raster(data = plot.df, mapping = aes(x = lon, y = lat, fill = xch4)) +
  scale_fill_gradientn(colours = ssec(100), limits = c(colorbar_min, colorbar_max)) +
  geom_point(data = receptor.df, mapping = aes(x = X, y = Y)) +
  geom_text(data = receptor.df, mapping = aes(x = X, y = Y, label = N), vjust = -1)
#ggsave(filename = "~/forward-model/figs/rf08_mosaic_with_receptors.png", device = png, width = 8, height = 8, units = "in")

# Plot the ID labeled data with the receptor grid superimposed
x11()
ggplot() +
  geom_point(data = all.reps.df, mapping = aes(x = lon, y = lat, color = as.character(id))) +
  scale_color_manual(values = rainbow(23)) +
  geom_point(data = receptor.df, mapping = aes(x = X, y = Y)) +
  geom_text(data = receptor.df, mapping = aes(x = X, y = Y, label = N), vjust = -1)
#ggsave(filename = "~/forward-model/figs/rf08_rainbow_with_receptors.png", device = png, width = 8, height = 8, units = "in")

# Plot just the receptors, labeled according to their nearest representative point
x11()
ggplot() +
  geom_point(data = receptor.df, mapping = aes(x = X, y = Y, color = as.character(id))) +
  scale_color_manual(values = rainbow(23)) +
  geom_text(data = receptor.df, mapping = aes(x = X, y = Y, label = N), vjust = -1)
#ggsave(filename = "~/forward-model/figs/rf08_receptors_id.png", device = png, width = 8, height = 8, units = "in")

# x11()
# ggplot() +
#   geom_point(data = receptor.df, mapping = aes(x = X, y = Y, color = as.character(id.time))) +
#   scale_color_manual(values = rainbow(23)) +
#   geom_text(data = receptor.df, mapping = aes(x = X, y = Y, label = N), vjust = -1)
# ggsave(filename = "~/forward-model/figs/rf08_rainbow_receptors.png", device = png, width = 8, height = 8, units = "in")

# x11()
# ggplot() +
#   geom_point(data = receptor.df, mapping = aes(x = X, y = Y, color = id.time)) +
#   scale_color_distiller(palette = "Spectral",
#                         breaks = as.numeric(receptor.df$id.time[c(1,50,100)]),
#                         labels = paste0(hour(receptor.df$id.time[c(1,50,100)]), ":", minute(receptor.df$id.time[c(1,50,100)])),
#                         name = "time")

# Plot the receptors, colored according to their timestamp with a continuous viridis colorscale.
x11()
ggplot() +
  geom_point(data = receptor.df, mapping = aes(x = X, y = Y, color = id.time), size = 3) +
  scale_colour_viridis(discrete = FALSE, 
                       option = "D", 
                       breaks = as.numeric(receptor.df$id.time[c(1, 100, 150, 200, 225)]), 
                       labels = paste0(hour(receptor.df$id.time[c(1, 100, 150, 200, 225)]), 
                                       ":", minute(receptor.df$id.time[c(1, 100, 150, 200, 225)])),
                       name = "time")
#ggsave(filename = "~/forward-model/figs/rf08_receptors_gradient.png", device = png, width = 8, height = 8, units = "in")

# ONCE VARIABLES/STRUCTURES ARE FINALIZED, SAVE WORKSPACE SO IT CAN BE EASILY LOADED IN THE FUTURE

# NEED TO CHECK TO MAKE SURE THAT LAT AND LON WERE ASSIGNED CORRECTLY, AND THAT THINGS DIDN'T SHIFT AMIDST THE IMAGE PROCESSING AND OTHER TRANSFORMATIONS





