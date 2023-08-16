###########################################################################!
# Extracting environmental data from Copernicus, not water column #########!
###########################################################################!
# Load the libraries
library(dplyr); library(raster); library(seegSDM)

#Load the database with the coordinates
data <- read.csv("./Data/Di_metrics_station.txt",header=T,dec=".",sep="\t", check.names = FALSE)

season<-read.csv("./Data/Raw_data/species_site_AFDW_2005.txt",header=T,dec=".",sep="\t")
#sample.mat <- dcast(season,station+year+ month +lon+lat~taxon, mean, value.var='wet_weight') #with this function we can create a matrix with the rows on the left side
season <- season %>% dplyr::select(station, year, month,lon, lat)

season$stat_year <- paste(season$station, season$year, sep = "_")
season <- season[which(!duplicated(season$stat_year)),]; season <- season %>% dplyr::select(station, year, month)
env <- data %>% dplyr::select(station, year, lon, lat)
env <- merge(env, season, by = c("station", "year")); env$month <- ifelse(env$month > 9, env$month, sprintf("%02d", env$month))
env$year_month <- paste(env$year,env$month, sep = ".")

#Prepare the dataset and extract the coordinates for which we want information
env <- env[which(env$year >= "2005"),] #select the years of interest
datasp <- env[,3:4]; coordinates(datasp)=~lon+lat #Select the coordinates to extract

# Example: Salinity at sea bottom - 2005-2020
salinity <- env
setwd("./env_data")
sal05 <- 'bot_sal_2005-2010.nc'; sal05 <- brick(sal05, var = "sob")
sal10 <- 'bot_sal_2010-2015.nc'; sal10 <- brick(sal10, var = "sob")
sal15 <- 'bot_sal_2015-2020.nc'; sal15 <- brick(sal15, var = "sob")
sal <- raster::stack(sal05,sal10,sal15) #Join all the files in one

var<- raster::extract(sal,datasp,method = 'bilinear', small = T, na.rm = T, df = T, exact = T, cellnumbers = F) #Extract the data for the selected coordinates

colnames(var) <- substr(colnames(var), 2,8) #Extract only the year and month for the name, expand to more characters if the day is also needed
var[1:5,1:5]

for (i in 1:nrow(salinity)){ #Obtain the median value for each month
  
  #salinity$bot_sal[i] <- var[i,which(colnames(var) == salinity$year_day[i])] #Environmental data for a specific day
  
  month <-var[i,which(colnames(var) == salinity$year_month[i])]
  salinity$bot_sal[i] <- apply(month, 1, median) #Median for a single month in a specific year
  
  # year <- var[i,which(substr(colnames(var), 1,4) == salinity$Year[i])]
  # salinity$bot_sal_var[i] <- apply(year, 1, sd) #Variability in a specific year
  
  month <- var[i,which(substr(colnames(var), 1,7) == salinity$year_month[i])]
  salinity$bot_sal_var[i] <- apply(month, 1, sd) #Variability in a specific month for a certain year
  
  
}

#Salinity: deal with NAs; 
salNA <- salinity[which(is.na(salinity$bot_sal)),];

#find the coordinates of the nearest nonNA cell in the raster
nas <- salNA[,3:4]
q = 25000 #we set a maximum distance for the function to look at coordinates, this is in meters
repeat{
  print(q)
  #With this nearestLand function, we calculate the nearest point in a raster file which is not NA
  #To deal with possible NA's in our data due to mistmatch between in-situ points and the satellite coordinates.
  land <- as.data.frame(nearestLand(nas, sal@layers[[1]], q)); colnames(land) <- c("lon", "lat")
  dup <- c(which(is.na(land)))
  q = q + 1000
  if(length(dup)==0){break}
} #We repeat it until there are not more NAs increasing the maximum distance by 1km each iteration

coordinates(land)= ~lon + lat
var<- raster::extract(sal,land,method = 'bilinear', small = T, na.rm = T, df = T, exact = T, cellnumbers = F);

colnames(var) <- substr(colnames(var), 2,8); 
for (i in 1:nrow(salNA)){
  
  #salNA$bot_sal[i] <- var[i,which(colnames(var) == salNA$year_day[i])]
  
  month <-var[i,which(colnames(var) == salNA$year_month[i])]
  salNA$bot_sal[i] <- apply(month, 1, median)
  
  # year <- var[i,which(substr(colnames(var), 1,4) == salNA$Year[i])]
  # salNA$bot_sal_var[i] <- apply(year, 1, sd)
  
  month <- var[i,which(substr(colnames(var), 1,7) == salNA$year_month[i])]
  salNA$bot_sal_var[i] <- apply(month, 1, sd)
  
} 
colnames(salNA)[7:8] <- c("sal","var_sal"); #Now we merge the original data with the nearest non-NA values
salinity <- salinity %>% full_join(salNA) %>% mutate(bot_sal = coalesce(bot_sal,sal), 
                                                     bot_sal_var = coalesce(bot_sal_var,var_sal)) %>% dplyr::select(-sal, -var_sal); rm(sal, sal05, sal10, sal15)
