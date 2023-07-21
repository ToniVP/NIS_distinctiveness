############################################################################!
# 3: Compute functional distinctiveness and other metrics at local scale, ##! 
#             with the corresponding local species pools                  ##!
############################################################################!
remove(list=ls())
#Load libraries and functions
source("./Code/Functions.R")
source("./Code/Libraries.R")

#Set the number of cores for parallel tasks
numcores <- detectCores()
registerDoParallel(numcores -1)

# 3.1: Functional distinctiveness and metrics for each single sampling event ----
#Load all the data sets
data<-read.csv("./Data/Raw_data/species_site_AFDW_2005.txt",header=T,dec=".",sep="\t", check.names = FALSE)
dist_matrix <- read.csv("./Data/dist_matrix_ovrll.txt",header=T,dec=".",sep="\t", check.names = FALSE)
traitraw <- read.csv("./Data/sp_traits.txt",header=T,dec=".",sep="\t", check.names = FALSE)
status <- read.csv("./Data/Raw_data/sp_status.txt",header=T,dec=".",sep="\t", check.names = FALSE) #Needed to identify the NIS
rownames(status) = status$taxon

#Be sure that all species are in the occurrences database
common_sp <- get_common_species(data, traitraw, 5); rownames(common_sp) <- gsub("\\.", " ", rownames(common_sp))
cols <- c(rownames(common_sp)); data <- cbind(data[1:5],data[,c(cols)])

#Calculate all the metrics for each single sampling event
#data <-  data[sample(nrow(data), 100),] #sample a subset to try it
#data = data[which(data$station == 911 & data$year == 2009),] #select a single station in case something does not work

#IMPORTANT: This is a loop that runs in parallel, within several cores in your computer, in order to reduce the computing time. 
# The number of cores of your computer destined to run this piece of code are set in the beginning; in case that you would not like to run this in parallel
# then "%dopar%" needs to be changed to "%do%". 
# For more information about parallel tasks in R, check this link: https://nceas.github.io/oss-lessons/parallel-computing-in-r/parallel-computing-in-r.html

system.time(stationDi <- foreach(t = 1:nrow(data), 
                                 .packages = c("dplyr", "vegan","data.table","raster",
                                                 "doParallel","foreach", "FD"), .errorhandling = "stop")%dopar%{ #start loop inside each station
   
   tryCatch({
     sample.mat <- data[t,]
     cols <- sample.mat[,6:ncol(sample.mat)]; num = 1
     
     #calculate NIS mean distinctiveness, distinctiveness x cell and distinctiveness x cell excluding NIS
     if (length(cols[,which(colSums(cols)> 0)]) > 1){cols <- cols[,which(colSums(cols)> 0)]} else {name = colnames(cols)[which(colSums(cols)> 0)]; 
     cols <- as.data.frame(cols[,which(colSums(cols)> 0)]); colnames(cols) <- name}
     
     #Di not weighted by local abundance/biomass
     results_notW <- mean_diss(dist_matrix, cols, status, num)
     Di_nW <- as.data.frame(results_notW[1]); Di_NAT_nW <- data.frame(non_NIS_Di =results_notW[[2]]); Di_ovrll_nW <- data.frame(overall_Di = results_notW[[3]])
     dist_nW <- cbind(Di_nW, Di_NAT_nW, Di_ovrll_nW)
     colnames(dist_nW) <- paste0(colnames(dist_nW), '_nw')
     
     #Di weighted by local abundance/biomass
     results <- mean_diss_rel_ab(dist_matrix, cols, status, num)
     Di <- as.data.frame(results[1]); Di_NAT <- data.frame(non_NIS_Di =results[[2]]); 
     Di_ovrll <- data.frame(overall_Di = results[[3]]); ratio_sp <- data.frame(ratio_sp = results[[4]])
     ratio <- as.data.frame(Di %>% dplyr::select(contains("_relab")) %>% rowSums(.)); colnames(ratio) <- "ratio_NISvsNAT"
     distances <- cbind(Di,Di_NAT,Di_ovrll,ratio, ratio_sp, dist_nW)
     
     #Functional diversity metrics, FRichness, FEveness, FDivergence, Rao's Q...
     traitcomm <- dist_matrix[c(colnames(cols)),c(colnames(cols))]
     if (ncol(cols) > 2){
       Fmetrics <- dbFD(traitcomm, messages = T, m = 10,
                        calc.FRic = TRUE, calc.FDiv = TRUE);
       
       FD <- lapply(Fmetrics[c(1:3,5:8)], mean, na.rm = TRUE);
       FDiv <- as.data.frame(do.call(cbind,FD))} else {FDiv <- NA}
     
     #Diversity metrics: Shannon, eveness...
     Sh<-diversity(cols, index = "shannon", MARGIN = 1, base =exp(1))
     S<-diversity(cols, index = "simpson", MARGIN = 1, base =exp(1))
     Richness <- ncol(cols) #the columns are all the species present in the location
     J<-Sh/log(Richness)
     divers <- data.frame(richness = Richness, Shannon = Sh, Simpson = S, Eveness = J)
     
     #Bind all the computed metrics
     metrics <- cbind(distances, FDiv)
     metrics <- cbind(metrics, divers)
     metrics <- metrics %>% select_if(~!all(is.na(.)))
     
   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
   
   
   data.frame(cbind (station = as.factor(data$station[t]), year = data$year[t], 
                     lon = data$lon[t],lat = data$lat[t], metrics))
 })

stationdi <- do.call(bind_rows, stationDi);
stationdi <- stationdi[colSums(!is.na(na_if(stationdi, 0))) > 0]
stationdi <- stationdi[which(stationdi$richness != 1),] #delete those stations that only have one species

#write.table(stationdi,file="./Data/Di_metrics_station.txt",sep="\t", row.names = TRUE)

# 3.2 Check the spatial autocorrelation for all the sampled locations; Moran's I test ----
Mor = data.frame()
autocorr <- c("non_NIS_Di", "mean_NIS_Di", "Marenzelleria", "Mya.arenaria", "Potamopyrgus.antipodarum","Streblospio.benedicti", "Polydora.cornuta")
for (i in 1:length(autocorr)){
  temp <- datat; colnames(temp)[which(colnames(datat) == autocorr[i])] <- "var"
  tempdata<- temp[which(!is.na(temp$var)),]
  data.distsXN <- as.matrix(dist(cbind(tempdata$lon,tempdata$lat)))
  data.dists.invXN <- 1/data.distsXN
  diag(data.dists.invXN) <- 0
  data.dists.invXN[is.infinite(data.dists.invXN)] <- 0
  tMor<-data.frame(Moran.I(tempdata$var, data.dists.invXN),Variable=autocorr[i])#data is auto correlated
  Mor<-rbind(tMor,Mor)}
#write.table(Mor,file="Moran_index.txt",sep="\t", row.names = TRUE)

# Fig. 5A-B; distribution of NIS and natives distinctiveness across the whole region ----
#Draw maps of distinctiveness to see how it's distributed. Linear interpolation of results
data<-read.csv("./Data/Di_metrics_station.txt",header=T,dec=".",sep="\t", check.names = FALSE)
datat <- data

metrics <- c("mean_NIS_Di","non_NIS_Di"); plots = list()

#NOTE: This process might not work some times, if it does not work, run the iterations manually
for (i in 1:length(metrics)){
  temp <- datat; colnames(temp)[which(colnames(temp) == metrics[i])] <- "var"
  temp <- temp[which(!is.na(temp$var)),]
  coords <- cbind(lon = temp$lon, lat = temp$lat)
  
  #Cut a polygon to crop the interpolated data
  myShape <- getDynamicAlphaHull(coords, partCount = 1, buff=500, clipToCoast = FALSE) #create a polygon to define the training area
  # plot(myShape[[1]])
  try <- temp
  
  #Krigging: linear interpolation
  coordinates(temp) <- ~ lon+lat
  kriging_result= autoKrige(var~1, temp)
  Kriging=as.data.frame(kriging_result$krige_output); Kriging<-fortify(Kriging)
  
  #Rasterize the linear interpolation to crop it for our interest area
  coords <- cbind(Kriging$x1, Kriging$x2)
  crs<- CRS("+proj=longlat")# proj4string of coords
  kd<-SpatialPointsDataFrame(coords = coords, data = Kriging,proj4string = crs)
  r<-raster(ncol=90, nrow=90, 
            xmn = min(Kriging$x1), ymn = min(Kriging$x2), xmx = max(Kriging$x1), ymx = max(Kriging$x2))
  
  krig <- rasterize(kd,r,field = kd$var1.pred, fun = mean, update=TRUE, updatevalue = "NA"); #plot(krig)
  shape_Di <- raster::crop(krig, myShape[[1]], snap = "out") #get the polygon to crop in, crop the rectangle
  train <- mask(x=krig, mask=myShape[[1]]); #plot(train)
  Kriging <- data.frame(rasterToPoints(train)); data <-try
  
  #Plot the results
  colpal<- c("midnightblue", "ghostwhite", "coral4") # Color palette for the map
  map <- ggplot(data, aes(x= lon, y = lat)) + 
    geom_tile(data=Kriging,aes(x,y,fill=layer))+
    scale_fill_gradientn(name = metrics[i], colours=colpal, na.value = 'white')+
    borders(fill="grey55", colour = "grey5") +
    theme(panel.background = element_rect(fill='lightsteelblue1', colour = 'black')) +
    coord_quickmap(xlim=c(min(data$lon)-1, max(data$lon)+1),ylim= c(min(data$lat)-1,max(data$lat)+1))+
    labs(x = "Lon",y="Lat")+
    theme(legend.position = c(1.2,0.2))+
    theme(legend.title = element_text(size = 11),legend.text = element_text(size = 10))+
    guides(shape = guide_legend(override.aes = list(size = 1))); map
  
  p <- list(map); names(p) <- metrics[i]
  name <- metrics[i]
  plots <- append(plots, p); 
  #ggsave(file=paste(name, ".svg", sep = ""), plot = p[[1]], width=10, height=10)
  
}
#plots

crossPlot <- ggarrange(plots[[1]], plots[[2]], 
                       ncol=1, nrow=2, common.legend = F, legend = "right"); crossPlot
#ggsave(file="description_metrics.svg", plot = crossPlot, width=20, height=20)


# Fig. S3; Local distinctiveness and distribution of NIS ----
NIS <- c("Marenzelleria", "Mya.arenaria", "Potamopyrgus.antipodarum","Streblospio.benedicti", "Polydora.cornuta"); plots = list()

#NOTE: This process might not work some times, if it does not work, run the iterations manually
for (i in 1:length(NIS)){
  temp <- datat; colnames(temp)[which(colnames(temp) == NIS[i])] <- "var"
  temp <- temp[which(!is.na(temp$var)),]
  coords <- cbind(lon = temp$lon, lat = temp$lat)
  
  #Cut a polygon to crop the interpolated data
  myShape <- getDynamicAlphaHull(coords, partCount = 1, buff=500, clipToCoast = FALSE) #create a polygon to define the training area
  # plot(myShape[[1]])
  try <- temp
  
  #Krigging: linear interpolation
  coordinates(temp) <- ~ lon+lat
  kriging_result= autoKrige(var~1, temp)
  Kriging=as.data.frame(kriging_result$krige_output); Kriging<-fortify(Kriging)
  
  if(NIS[i] == "Mya.arenaria"){ncol = 50} else {ncol = 90}
  
  #Rasterize the linear interpolation to crop it for our interest area
  coords <- cbind(Kriging$x1, Kriging$x2)
  crs<- CRS("+proj=longlat")# proj4string of coords
  kd<-SpatialPointsDataFrame(coords = coords, data = Kriging,proj4string = crs)
  r<-raster(ncol= ncol, nrow= ncol, 
            xmn = min(Kriging$x1), ymn = min(Kriging$x2), xmx = max(Kriging$x1), ymx = max(Kriging$x2))
  
  krig <- rasterize(kd,r,field = kd$var1.pred, fun = mean, update=TRUE, updatevalue = "NA"); #plot(krig)
  shape_Di <- raster::crop(krig, myShape[[1]], snap = "out") #get the polygon to crop in, crop the rectangle
  train <- mask(x=krig, mask=myShape[[1]]); #plot(train)
  Kriging <- data.frame(rasterToPoints(train)); data <-try
  
  #Plot the distinctiveness
  colpal<- c("midnightblue", "ghostwhite", "coral4") # Color palette for the map
  Di_map <- ggplot(data, aes(x= lon, y = lat)) + 
    geom_tile(data=Kriging,aes(x,y,fill=layer))+
    scale_fill_gradientn(name = NIS[i], colours=colpal, na.value = 'white')+
    borders(fill="grey55", colour = "grey5") +
    theme(panel.background = element_rect(fill='lightsteelblue1', colour = 'black')) +
    coord_quickmap(xlim=c(min(data$lon)-1, max(data$lon)+1),ylim= c(min(data$lat)-1,max(data$lat)+1))+
    labs(x = "Lon",y="Lat")+
    theme(legend.position = c(1.2,0.2))+
    theme(legend.title = element_text(size = 11),legend.text = element_text(size = 10))+
    guides(shape = guide_legend(override.aes = list(size = 1))); #Di_map
  
  
  #Plot the distribution
  if(NIS[i] == "Marenzelleria"){colpal<- "red"}
  if(NIS[i] == "Mya.arenaria"){colpal<- "darkorange1" }
  if(NIS[i] == "Potamopyrgus.antipodarum"){colpal<- "navy" }
  if(NIS[i] == "Streblospio.benedicti"){colpal<- "mediumpurple" }
  if(NIS[i] == "Polydora.cornuta"){colpal<- "limegreen" }
  
   dist <- ggplot(data, aes(x= lon, y = lat)) + 
    geom_point(color = colpal, size = 3)+
    borders(fill="grey55", colour = "grey5") +
    theme(panel.background = element_rect(fill='lightsteelblue1', colour = 'black')) +
    coord_quickmap(xlim=c(min(data$lon)-1, max(data$lon)+1),ylim= c(min(data$lat)-1,max(data$lat)+1))+
    labs(x = "Lon",y="Lat")+
    theme(legend.position = c(0.84,0.16))+
    theme(legend.title = element_text(size = 11),legend.text = element_text(size = 10))+
    guides(shape = guide_legend(override.aes = list(size = 1))); dist
  
  p <- list(list(Di_map,dist)); names(p) <- NIS[i]
  plots <- append(plots, p)
  #ggsave(file=paste(name, ".svg", sep = ""), plot = p[[1]], width=10, height=10)
  
}
#plots

crossPlot <- ggarrange(plots$Polydora.cornuta, plots$Mya.arenaria, plots$Marenzelleria, plots$Potamopyrgus.antipodarum, plots$Streblospio.benedicti,
                       ncol=3, nrow=2, common.legend = F, legend = "none"); crossPlot
#ggsave(file="description_metrics.svg", plot = crossPlot, width=20, height=20)




















