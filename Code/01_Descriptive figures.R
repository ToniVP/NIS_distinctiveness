###############################!
#1: Descriptive figures #######!
###############################!
remove(list=ls())
#Load libraries and functions
source("./Code/Functions.R")
source("./Code/Libraries.R")

# Fig. 2: Description maps ----
# NOTE: The stations map (panel A) was saved in vector format and then the names of the regions were added manually in Adobe Illustrator.
# The metrics shown in the other panels (B-E) were obtained later in the process, the full process is detailed in the script "03_Distinctiveness_local"

# Fig. 2A: Stations map ----
data<-read.csv("./Data/Raw_data/species_AFDW_2005_2020.txt",header=T,dec=".",sep="\t", check.names = FALSE)

data$year <- as.factor(data$year)

p <- ggplot(data, aes(x = lon,y = lat)) + 
  geom_point(size = 1.5, color = "darkgreen") +
  borders(fill="grey55", colour = "grey5") +
  theme(panel.background = element_rect(fill='lightsteelblue1', colour = 'black')) +
  coord_quickmap(xlim=c(min(data$lon), max(data$lon)),ylim= c(min(data$lat),max(data$lat)))+
  labs(x = "Lon",y="Lat")+
  theme(legend.position = "none")+
  theme(legend.title = element_text(size = 11),legend.text = element_text(size = 10))+
  guides(shape = guide_legend(override.aes = list(size = 1))); p

#ggsave(file="map_stations.svg", plot=p, width=8, height=10) 

# Fig. 2B-E: Richness, Shannon, relative biomass of NIS and % of NIS in local communities ----
data<-read.csv("./Data/Di_metrics_station.txt",header=T,dec=".",sep="\t", check.names = FALSE)
datat <- data

#Linear interpolation for all the variables
metrics <- c("richness", "Shannon", "ratio_NISvsNAT", "ratio_sp"); plots = list()

#NOTE: This process might not work some times, if it does not work, run the iterations manually
for (i in 1:length(metrics)){
  temp <- datat; colnames(temp)[which(colnames(temp) == metrics[i])] <- "var"
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

crossPlot <- ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], 
                       ncol=2, nrow=2, common.legend = F, legend = "right"); crossPlot
#ggsave(file="description_metrics.svg", plot = crossPlot, width=20, height=20)

# Fig. S1: Species biomass ----
data<-read.csv("./Data/Raw_data/species_AFDW_2005_2020.txt",header=T,dec=".",sep="\t", check.names = FALSE)
status<-read.csv("./Data/Raw_data/sp_status.txt",header=T,dec=".",sep="\t", check.names = FALSE)

# Biomass (AFDW) of all the species
biomass_cum <- as.data.frame(data %>% 
                               group_by(taxon) %>% 
                               summarise_at(.vars = "AFDW", sum, na.rm = TRUE))

biomass <- merge(biomass_cum, status, by="taxon"); biomass$rel_abun <- (biomass$AFDW/sum(biomass$AFDW))*100

# Donut plot for the relative abundances - Natives vs NIS ----
#Compute the cumulative percentages (top of each rectangle)
biom <- as.data.frame(biomass %>% 
                        group_by(status) %>% 
                        summarise_at(.vars = c("rel_abun","AFDW"), sum, na.rm = TRUE)); biom$ymax = cumsum(biom$rel_abun)

# Compute the rectangle bottom
biom$ymin = c(0, head(biom$ymax, n=-1))
#Comptue label position
biom$labelPosition <- (biom$ymax + biom$ymin) / 2
biom$label1 <- paste0(round(biom$rel_abun, digits = 2), "%")

# Make the plot-NIS VS NATIVE
natnis <- ggplot(biom, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=status)) +
  geom_rect(size = 1.5, aes(color = status, alpha = 0.5,
                            fill = after_scale(desaturate(lighten(color, .5), .5)))) +
  geom_text(x=2.2, aes(y = labelPosition,label=status, color=status, fontface = "bold"),size=11, check_overlap = TRUE) + # x here controls label position (inner / outer)
  geom_label(x=3.5, aes(y = labelPosition,label=label1, color = status, 
                        fill = after_scale(desaturate(lighten(color, .9), .9))),size=12) +
  scale_fill_aaas() +
  scale_color_aaas() +
  coord_polar(theta="y") +
  theme_void() +
  theme(panel.grid=element_blank()) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank(),
        panel.border = element_blank()) + theme(legend.position = "none"); natnis


#ggsave(file="NATNIS_AFDW.svg", plot=natnis, width=10, height=8)

# Donut plot for the relative abundances - NATIVES ----
nat <- biomass[which(biomass$status == "Native"),] #Ratio groups of native organisms
nat <- as.data.frame(biomass %>% 
                       group_by(phylum) %>% 
                       summarise_at(.vars = c("rel_abun","AFDW"), sum, na.rm = TRUE))

nat$rel_abun <- (nat$AFDW/sum(nat$AFDW))*100; 
nat_other <- nat[which(nat$rel_abun <=1),]; nat <- nat %>% add_row(phylum = "Other", AFDW = sum(nat_other$AFDW),
                                                                   rel_abun = sum(nat_other$rel_abun))
nat <- anti_join(nat, nat_other); nat$ymax = cumsum(nat$rel_abun)

# Compute the rectangle bottom
nat$ymin = c(0, head(nat$ymax, n=-1))
#Comptue label position
nat$labelPosition <- (nat$ymax + nat$ymin) / 2
nat$label1 <- paste0(round(nat$rel_abun, digits = 2), "%")

#Make the plot biomass NATIVES
NAT = ggplot(nat, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=phylum)) +
  geom_rect(size = 1.5, aes(color = phylum, alpha = 0.6, fill = after_scale(desaturate(lighten(color, .4), .4)))) +
  scale_fill_uchicago(palette = "light") +
  # x here controls label position (inner / outer)
  geom_text_repel(min.segment.length = 0,
                  aes(x = 3, y = labelPosition,label=phylum, colour = factor(phylum), fontface = "bold.italic",
                      segment.size = 0.8, segment.linetype = 3), size=14, max.overlaps = Inf)+
  geom_label_repel(x=3.7, aes(y = labelPosition, label=label1, color = phylum, size  = 12,
                              fill = after_scale(desaturate(lighten(color, .9), .9))),
                   size=5, fontface = "bold")+
  scale_color_uchicago(palette = "dark") +
  coord_polar(theta="y", start = 0) +
  theme_void() +
  theme(panel.grid=element_blank()) +
  theme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank(),
        panel.border = element_blank()) + theme(legend.position = "none"); NAT

#ggsave(file="NAT_AFDW.svg", plot=NAT, width=10, height=8)

# Donut plot for the relative abundances - NIS ----
nis <- biomass[which(biomass$status == "Non-indigenous"),] #Ratio groups of native organisms

nis$rel_abun <- (nis$AFDW/sum(nis$AFDW))*100; nis$ymax = cumsum(nis$rel_abun)

# Compute the rectangle bottom
nis$ymin = c(0, head(nis$ymax, n=-1))
#Comptue label position
nis$labelPosition <- (nis$ymax + nis$ymin) / 2
nis$label1 <- paste0(round(nis$rel_abun, digits = 2), "%")

#Make the plot biomass NIS
NIS <- ggplot(nis, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=taxon)) +
  geom_rect(size = 1.5, aes(color = taxon, alpha = 0.5, fill = after_scale(desaturate(lighten(color, .4), .4)))) +
  scale_fill_uchicago(palette = "light") +
  # x here controls label position (inner / outer)
  geom_text_repel(min.segment.length = 0,
                  aes(x = 3, y = labelPosition,label=taxon, colour = factor(taxon), fontface = "bold",
                      segment.size = 0.8, segment.linetype = 3),
                  size=7, max.overlaps = Inf)+
  geom_label_repel(x=3.7, aes(y = labelPosition, label=label1, color = taxon, 
                              fill = after_scale(desaturate(lighten(color, .9), .9))),
                   size=6, fontface = "bold")+
  scale_color_uchicago(palette = "dark") +
  coord_polar(theta="y") +
  theme_void() +- theme(panel.grid=element_blank()) +vtheme(axis.text=element_blank()) +
  theme(axis.ticks=element_blank(), panel.border = element_blank()) + theme(legend.position = "none"); NIS

#ggsave(file="NIS_AFDW.svg", plot=NIS, width=10, height=8)

# Fig. S4; Maps of environmental variables ----
data<-read.csv("./Data/Di_metrics_station.txt",header=T,dec=".",sep="\t", check.names = FALSE); coln <- c(colnames(data))

NIS <- c("Marenzelleria", "Mya.arenaria", "Potamopyrgus.antipodarum","Streblospio.benedicti", "Polydora.cornuta")#put the names of your NIS here
colnames(data)
Di_unified <- data.frame()

# Put the species as a factor rather than separate columns -- if necessary
for(i in 1:length(NIS)){
  
  temp <- data %>% dplyr::select(c(1:21), mean_NIS_Di, NIS[i], paste(NIS[i],"relab", sep = "_"), paste(NIS[i],"nw", sep = "_"));
  temp <- subset(temp, !is.na(temp[,NIS[i]]))
  colnames(temp)[c(23,24,25)] <- c("NIS_Di", "NIS_relab", "NIS_Di_nw"); temp$species <- NIS[i]
  Di_unified <- rbind(Di_unified, temp)
  
}

#Some data manipulation, join the Di in one column and add species as a factor
data <- Di_unified; data$species <- as.factor(data$species);  coln <- c(colnames(data))

#Join the environmental data with the metrics
env <- read.csv("./Data/Raw_data/env_data.txt",header=T,dec=".",sep="\t", check.names = FALSE)
head(env); env <- env %>% dplyr::select(station, year, lon, lat, month, bot_oxy, bot_sal, bot_T, depth)

data <- merge(data, env, by = c("station", "year", "lon", "lat")); #data$Sediment <- as.factor(data$Sediment)
data$station <- as.factor(data$station)

data <- data[which(!duplicated(data[,coln])),]

#Loop for the linear interpolation of each variable
metrics <- c("bot_sal", "bot_T", "depth", "bot_oxy"); plots = list(); datat = data

#NOTE: This process might not work some times, if it does not work, run the iterations manually
for (i in 1:length(metrics)){
temp <- datat; colnames(temp)[which(colnames(temp) == metrics[i])] <- "var"
temp <- temp[which(!is.na(temp$var)),]
coords <- cbind(lon = temp$lon, lat = temp$lat)

#Cut a polygon to crop the interpolated data
myShape <- getDynamicAlphaHull(coords, partCount = 1, buff=500, clipToCoast = FALSE) #create a polygon to define the training area
#plot(myShape[[1]])
try <- temp

#Krigging: linear interpolation
coordinates(temp) <- ~ lon+lat
kriging_result= autoKrige(var~1, temp)
Kriging=as.data.frame(kriging_result$krige_output); Kriging<-fortify(Kriging)

#Rasterize the linear interpolation to crop it for our interest area
coords <- cbind(Kriging$x1, Kriging$x2)
crs<- CRS("+proj=longlat")# proj4string of coords
kd<-SpatialPointsDataFrame(coords = coords,data = Kriging,proj4string = crs)
r<-raster(ncol=90, nrow=90, 
          xmn = min(Kriging$x1), ymn = min(Kriging$x2), xmx = max(Kriging$x1), ymx = max(Kriging$x2))
krig <- rasterize(kd,r,field = kd$var1.pred, fun = mean, update=TRUE, updatevalue = "NA"); #plot(krig)
shape_Di <- raster::crop(krig, myShape[[1]], snap = "out") #get the polygon to crop in, crop the rectangle
train <- mask(x=krig, mask=myShape[[1]]); #plot(train)
Kriging <- data.frame(rasterToPoints(train)); temp <-try

#Plot the results
if(metrics[i] == "bot_sal"){colpal<- c("ghostwhite","cornflowerblue", "midnightblue") }
if(metrics[i] == "bot_oxy"){colpal<- c("ghostwhite","cornflowerblue", "midnightblue") }
if(metrics[i] == "bot_T"){colpal<- c("cornflowerblue","ghostwhite", "red") }
if(metrics[i] == "depth"){colpal<- c("ghostwhite", "plum1", "darkorchid4") }

env_map <- ggplot(temp, aes(x= lon, y = lat)) + 
  geom_tile(data=Kriging,aes(x,y,fill=layer))+
  scale_fill_gradientn(name = metrics[i], colours=colpal, na.value = 'white')+
  borders(fill="grey55", colour = "grey5") +
  theme(panel.background = element_rect(fill='ghostwhite', colour = 'black')) +
  coord_quickmap(xlim=c(min(data$lon)-1, max(data$lon)+1),ylim= c(min(data$lat)-1,max(data$lat)+1))+
  labs(x = "Lon",y="Lat")+
  theme(legend.position = "right")+
  theme(legend.title = element_text(size = 11),legend.text = element_text(size = 10))+
  guides(shape = guide_legend(override.aes = list(size = 1))); env_map

p <- list(env_map); names(p) <- metrics[i]
name <- metrics[i]
plots <- append(plots, p) 

}

crossPlot <- ggarrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]], ncol=2, nrow=2, common.legend = F, legend = "right"); crossPlot

#ggsave(file="maps_env.svg", plot=crossPlot, width=10, height=11)

legend <- cowplot::get_legend(bs)
grid.draw(legend)
#ggsave(file="legendbotT.svg", plot=legend, width=10, height=8)

# Fig. S6; Depth distribution for all NIS ----

plots = list()
for(i in 1:length(NIS)){
  
  temp <- data[which(data$species == NIS[i]),]
  
  if(NIS[i] == "Marenzelleria"){fill = "rosybrown1"; col = "darkred"}
  if(NIS[i] == "Mya.arenaria"){fill = "lightgoldenrod" ; col = "darkgoldenrod"}
  if(NIS[i] == "Potamopyrgus.antipodarum"){fill = "skyblue1" ; col = "navy" }
  if(NIS[i] == "Streblospio.benedicti"){fill<- "orchid1"; col = "darkorchid" }
  if(NIS[i] == "Polydora.cornuta"){fill<- "lightgreen"; col = "darkgreen" }
  
  depth_NIS <- ggplot(temp, aes(x=depth)) + 
    geom_bar(stat="bin", bins = 15, fill = fill, color = col) +
    labs(x = "Depth (m)", y = "Occurrences") + ggtitle(NIS[i]) +
    theme_minimal() + theme(legend.position = "right", text = element_text(size = 15), axis.text.y = element_text(size = 15),
                            axis.text.x = element_text(size = 15)); depth_NIS
  
  p <- list(depth_NIS); names(p) <- NIS[i]
  name <- NIS[i]
  plots <- append(plots, p)
  
}

depth_all <- ggplot(data, aes(x=depth)) + 
  geom_bar(stat="bin", bins = 15, fill ="gray90", color = "black" ) +
  labs(x = "Depth (m)", y = "Occurrences") + ggtitle("All NIS") +
  theme_minimal() + theme(legend.position = "right", text = element_text(size = 15), axis.text.y = element_text(size = 15),
                          axis.text.x = element_text(size = 15)); depth_all

plots = append(list(depth_all), plots)

d <- grid.arrange(grobs = plots, nrow = 2); d
#ggsave(file="depth_NIS.svg", plot=d, width=10, height=8)

