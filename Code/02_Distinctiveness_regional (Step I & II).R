############################################################################!
# 2: Compute functional distinctiveness for the regional species pool ######!
############################################################################!
remove(list=ls())
#Load libraries and functions
source("./Code/Functions.R")
source("./Code/Libraries.R")

#Set the number of cores for parallel tasks
numcores <- detectCores()
registerDoParallel(numcores -1)

# 2.0: Prepare the data for the calculation of the distance matrix ----
# 2.0.1: Trait data cleaning ----
#Put all the traits and modalities in a list
TRAITS <- list ("Size" = c("0-10","11-20","21-50","51-100",">100"),
                "Reproductive type"= c("Semelparous-Monotelic","Iteroparous-Polytelic","Semi-continous"),
                "Adult life span" = c("<1yr","1-3yrs","3-6yrs","6-10yrs",">10yrs"),
                "Developmental mechanism" = c("Fragmentation/Fission","Direct","Lecitotrophic","Planktotrophic","Ovoviviparous"),
                "Environmental position" = c("Deep","Middle","Top","Interface","Epibenthic","Bentho-pelagic","Epilithic", "Epifaunal"),
                "Living Habits" =  c("Attached","Tube dweller","Burrow dweller","Crevic dweller","Case builder","Free","Parasite/commensal"),
                "Feeding method" = c("Suspension/filter feeder","Deposit feeder (inkl. Both)","Predator","Scavenger","Herbivore","Miner/Borer","Parasite", "Grazer"),
                "Mobility" = c("Sessile","Semi-motile","Motile"),
                "Movement method" = c("No movement","Swimmer","Crawler","Rafter/Drifter/Byssus","Tube-builder","Burrower","Temporary attachment"),
                "Bioturbation" = c("No transport","Diffusive mixing","Surface deposition","Conveyor belt transport","Reverse conveyor belt transport"))
#saveRDS(TRAITS, "Trait_modalities")

#Change the scores for categories on the traits that only have one modality
traits<-read.csv("./Data/Raw_data/sp_traits_raw.txt",header=T,dec=".",sep="\t", check.names = FALSE)
traits$taxon = rownames(traits)
TRAITS <- readRDS("./Data/Trait_modalities")

#Size - convert modalities to ordinal
a <- unlist(TRAITS$Size)
prova <- traits %>% gather(key = "trait", value = "dummy", all_of(a)) #Extract the value 0 or 1 per each species and trait modality
prova <- prova[which(prova$dummy == 1),] #select only those species which have at least one modality present
size <- prova %>% 
  group_by(taxon) %>% 
  summarise_at(.vars = "trait", funs(trimws(paste(., collapse = ';')))); colnames(size)[2] <- "Size"

traits <- merge(traits, size, by = "taxon"); traits$Size <- as.factor(traits$Size)
levels(traits$Size) <- c(5,1,2,3,4); traits$Size <- as.numeric(as.character(traits$Size))

#Adult life span
a <- unlist(TRAITS$`Adult life span`)
prova <- traits %>% gather(key = "trait", value = "dummy", all_of(a)) #Extract the value 0 or 1 per each species and trait modality
prova <- prova[which(prova$dummy == 1),] #select only those species which have at least one modality present
lifespan <- prova %>% 
  group_by(taxon) %>% 
  summarise_at(.vars = "trait", funs(trimws(paste(., collapse = ';')))); colnames(lifespan)[2] <- "Adult life span"

traits$`Adult life span` <- NULL; traits <- merge(traits, lifespan, by = "taxon"); traits$`Adult life span` <- as.factor(traits$`Adult life span`)
levels(traits$`Adult life span`) <- c(1,5,2,3,4); traits$`Adult life span` <- as.numeric(as.character(traits$`Adult life span`))

# Remove the previous categories from the trait database
size <- unlist(TRAITS$Size); lifespan <- unlist(TRAITS$`Adult life span`)

traits <- traits %>% dplyr::select(-size,-lifespan); rownames(traits) <- traits$taxon; traits$taxon = NULL

#Data with the species to include in as columns
data<-read.csv("./Data/Raw_data/species_site_AFDW_2005.txt",header=T,dec=".",sep="\t", check.names = FALSE) 

# We need to be sure that the species in the species x site matrix are in the trait database
trait <- get_common_species(data,traits, num = 5)

# 2.1: Functional distinctiveness at regional scale (regional species pool) ----

#Compute the overall pairwise functional distances matrix between all species
TRAITS[1] <- "Size"; TRAITS[3] <- "Adult life span"

dist_matrix <- funct_distances(TRAITS, trait); dist_matrix[1:5,1:5]
#write.table(dist_matrix,file="./Data/dist_matrix_ovrll.txt",sep="\t",row.names = TRUE)

# 2.1.1: Distinctiveness NOT WEIGHTED BY RELATIVE ABUNDANCE----------------
dist_matrix<-read.csv("./Data/dist_matrix_ovrll.txt",header=T,dec=".",sep="\t", check.names = FALSE)
standard.Di.gow<-colSums(dist_matrix)/(nrow(dist_matrix)-1) #Mean of the computed distances for each species, to obtain the

#average distance to the other species based on all trait combinations
Int_Di_nw<-as.data.frame(standard.Di.gow)

Int_Di_nw$taxon<-rownames(Int_Di_nw); colnames(Int_Di_nw)[1]<-c("int_Di")

# 2.1.2: Distinctiveness WEIGHTED BY RELATIVE ABUNDANCE (N individuals // biomass (AFDW)) ----------------
data<-read.csv("./Data/Raw_data/species_AFDW_2005_2020.txt",header=T,dec=".",sep="\t", check.names = FALSE) #Load the abundance database to compute the relative abundance of each species

# Abundance (N individuals), this was done only for Figure 3A ----
rel_abund <- as.data.frame(data %>% 
                             group_by(taxon) %>% 
                             summarise_at(.vars = c("abundance"), sum, na.rm = TRUE))

rel_abund$rel.abund <- (rel_abund$abundance/sum(rel_abund$abundance))*100 #Use the preferred abundance metric

rel_abund <- rel_abund[which(!is.na(match(rel_abund$taxon,rownames(trait)))),] #Make sure that the species match between the trait data and the relative abundance
rownames(rel_abund) <- rel_abund$taxon; 
rel_abund$taxon <- NULL; rel_abund$abundance <- NULL; rel_abund$AFDW <- NULL #Make sure to delete everything but the relative abundance
colnames(rel_abund) <- "rel.abund"

dist_mat_w <- matrix(0, nrow(dist_matrix), ncol(dist_matrix)) 
colnames(dist_mat_w) <- colnames(dist_matrix); rownames(dist_mat_w) <- rownames(dist_matrix)

#Here we multiply the functional pairwise distance by the relative abundance of each species
for (j in 1:ncol(dist_matrix)){ 
  for(i in 1:nrow(dist_matrix)){
    temp <- dist_matrix[i,j]*rel_abund[which(rownames(rel_abund) == rownames(dist_matrix[i,])),]
    dist_mat_w[i,j] <- temp
  }
}

weight.dist <- data.frame(dist = colSums(dist_mat_w))

Int_Di_ab <- data.frame(); 

#Distinctiveness weighted by species abundance
for(i in 1:ncol(dist_mat_w)){
  sumDi <- weight.dist[which(rownames(weight.dist) == colnames(dist_mat_w)[i]),] #Here we multiply the functional distance 
  sumrelAB <- sum(rel_abund[-which(rownames(rel_abund) == colnames(dist_mat_w)[i]),])
  
  Di<-sumDi/sumrelAB
  
  temp <- data.frame(taxon = colnames(dist_mat_w)[i], Int_Di = Di)
  Int_Di_ab <- rbind(Int_Di_ab, temp)
}

colnames(Int_Di_ab)<-c("taxon","int_Di")

# Biomass (AFDW) ----
rel_abund <- as.data.frame(data %>% 
                             group_by(taxon) %>% 
                             summarise_at(.vars = c("AFDW"), sum, na.rm = TRUE))

rel_abund$rel.abund <- (rel_abund$AFDW/sum(rel_abund$AFDW))*100 #Use the preferred abundance metric

rel_abund <- rel_abund[which(!is.na(match(rel_abund$taxon,rownames(trait)))),] #Make sure that the species match between the trait data and the relative abundance
rownames(rel_abund) <- rel_abund$taxon; 
rel_abund$taxon <- NULL; rel_abund$abundance <- NULL; rel_abund$AFDW <- NULL #Make sure to delete everything but the relative abundance
colnames(rel_abund) <- "rel.abund"

dist_mat_w <- matrix(0, nrow(dist_matrix), ncol(dist_matrix)) 
colnames(dist_mat_w) <- colnames(dist_matrix); rownames(dist_mat_w) <- rownames(dist_matrix)

#Here we multiply the functional pairwise distance by the relative abundance of each species
for (j in 1:ncol(dist_matrix)){ 
  for(i in 1:nrow(dist_matrix)){
    temp <- dist_matrix[i,j]*rel_abund[which(rownames(rel_abund) == rownames(dist_matrix[i,])),]
    dist_mat_w[i,j] <- temp
  }
}

weight.dist <- data.frame(dist = colSums(dist_mat_w))

Int_Di_AFDW <- data.frame(); 

#Distinctiveness weighted by species AFDW
for(i in 1:ncol(dist_mat_w)){
  sumDi <- weight.dist[which(rownames(weight.dist) == colnames(dist_mat_w)[i]),] #Here we multiply the functional distance 
  sumrelAB <- sum(rel_abund[-which(rownames(rel_abund) == colnames(dist_mat_w)[i]),])
  
  Di<-sumDi/sumrelAB
  
  temp <- data.frame(taxon = colnames(dist_mat_w)[i], Int_Di = Di)
  Int_Di_AFDW <- rbind(Int_Di_AFDW, temp)
}

colnames(Int_Di_AFDW)<-c("taxon","int_Di")

#Include the deciles defining the most/less distinct species in the regional pool
spe_index = mutate(Int_Di_AFDW, quartile = as.factor(ntile(Int_Di_AFDW$int_Di,4)),decile = as.factor(ntile(Int_Di_AFDW$int_Di,10)))
#write.table(spe_index,file="./Data/spe_index.txt",sep="\t")

# Fig. 3A; distribution of distinctiveness values based on different weightings ----
Int_Di_nw$weight <- "Non-weighted"; Int_Di_AFDW$weight <- "Biomass (AFDW)"; Int_Di_ab$weight <- "Abundance (N individuals)"
Int <- rbind(Int_Di_nw, Int_Di_AFDW, Int_Di_ab)

b <- Int %>%
  ggplot( aes(x=int_Di, group = weight, fill = weight, color  =weight)) + 
  xlim(0,0.5) + scale_color_aaas()+ scale_fill_aaas() +
  geom_density(alpha = 0.1, size = 1.5) + labs(x = "Functional distinctiveness (Di)", y = "Density") +
  theme_minimal() + 
  theme(legend.background = element_rect(fill = "transparent", color = NA), legend.position = "bottom", legend.title = element_blank(),
        text = element_text(size = 20), axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank()); b

#ggsave(file="density_distinctiveness_weight.svg", plot=b, width=10, height=8)

# Fig. 3B: Plot where the NIS are placed within D distribution ----
status<-read.csv("./Data/Raw_data/sp_status.txt",header=T,dec=".",sep="\t", check.names = FALSE)

p <- spe_index %>%
  ggplot( aes(x=int_Di)) +
  xlim(0,0.5) +
  geom_vline(aes(xintercept = median(int_Di)), size = 1.5, linetype = "solid", color = "black")+
  geom_vline(as.data.frame(spe_index[which(!is.na(match(spe_index$taxon,
                                                        status[which(status$status == "Non-indigenous"),"taxon"]))),]),
             mapping = aes(xintercept = int_Di, color = taxon), size = 1.2, linetype = "dashed") +
  scale_color_manual(values= c("firebrick1", "orange2", "chartreuse3", "royalblue3", "thistle3")) +
  geom_density(alpha = 0.1, size = 1, color = "black", fill = "gray60") + labs(x = "Functional distinctiveness (Di)", y = "Density") +
  theme_minimal() + 
  theme(legend.background = element_rect(fill = "transparent", color = NA), legend.position = "none", legend.title = element_blank(),
        text = element_text(size = 20), axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank()); p

#ggsave(file="density_distinctiveness_sp_weighting.svg", plot=crossPlot, width=10, height=16)

legend <- cowplot::get_legend(p)
grid.draw(legend)
#ggsave(file="legend_dens_sp.svg", plot=legend, width=10, height=8)

# 2.1.3: Figure S2; Sensitivity test about using different trait combinations ----

sensitivity_list <- list()
# Defining the total number of column : N; it MUST be 
N<-length(TRAITS)
listcol<-seq(1,N,1)
# Number of columns corresponding to traits
listk<-seq(1,N,1)
traits_N <- c(4,6,8); TRAITS[1] <- "Size"; TRAITS[3] <- "Adult life span"

# Integrated distinctiveness with different number of trait combinations (Fig. S2A-C) ----------
for (t in 1:length(traits_N)){ #start loop combination of traits
  
  k = traits_N[t]
  print (k)
  
  #Here are determined all possible 4 trait combinations for obtaining the average distance
  #between species, here are expressed only the position that are occupying the traits in the data frame
  dist.matr <- list()
  combination<-data.frame(t(combn(listk,k))) 
  
  if (is.null(combination$X5==T)){
    combination$X5<-NA
  }
  if (is.null(combination$X6==T)){
    combination$X6<-NA
  }
  if (is.null(combination$X7==T)){
    combination$X7<-NA
  }
  if (is.null(combination$X8==T)){
    combination$X8<-NA
  }
  if (is.null(combination$X9==T)){
    combination$X9<-NA
  }
  if (is.null(combination$X10==T)){
    combination$X10<-NA
  }
  if (is.null(combination$X11==T)){
    combination$X11<-NA
  }
  
  #This list of conditionals should be equal or have more elements than your trait number, NEVER less
  
  for (i in 1:nrow(combination)){ 
    data <- trait
    codex<-as.vector(t(combination[i,]))
    new_traits <- TRAITS[codex]
    codex <- unlist(new_traits)
    
    #Here, each row (in other words, each combination of 4 to 10 traits) is selected for a 
    #total calculation of distances based on all combinations possibles of the traits
    data<-dplyr::select(data,codex)
    
    #Calculation of distinctiveness using Gower's distance
    gow<-compute_dist_matrix(data, metric = "gower") # a distance matrix is calculated for each combination of traits
    m.gow<-as.matrix(gow)
    d.gow<-as.data.frame(m.gow)
    dist.matr <- c(dist.matr,list(d.gow))
  }
  
  dist_matrix <- as.data.frame(mean_matrix(dist.matr))
  
  standard.Di.gow<-colSums(dist_matrix)/(nrow(dist_matrix)-1) #Mean of the computed distances for each species, to obtain the
  #average distance to the other species based on all trait combinations
  Int_Di<-as.data.frame(standard.Di.gow)
  
  Int_Di$taxon<-rownames(Int_Di); colnames(Int_Di)[1]<-c("int_Di")

  spe_index = mutate(Int_Di, quartile = as.factor(ntile(Int_Di$int_Di,4)),
                     decile = as.factor(ntile(Int_Di$int_Di,10)), combination = as.factor(paste("Number_traits_combinations_", k)))
  spe_dist <- list(list(spe_index, dist_matrix)); 
  sensitivity_list[t] <- spe_dist; names(sensitivity_list)[t] <- paste("Number_traits_combinations_", k)
  
}

# Integrated distinctiveness with all the trait combinations together (Fig. S2D) ---------
dist_matrix <- funct_distances(TRAITS,trait);
standard.Di.gow<-colSums(dist_matrix)/(nrow(dist_matrix)-1) #Mean of the computed distances for each species, to obtain the
#average distance to the other species based on all trait combinations
Int_Di<-as.data.frame(standard.Di.gow)

Int_Di$taxon<-rownames(Int_Di); colnames(Int_Di)[1]<-c("int_Di")

spe_index = mutate(Int_Di, quartile = as.factor(ntile(Int_Di$int_Di,4)),decile = as.factor(ntile(Int_Di$int_Di,10)),
                   combination = as.factor("Integrated"))

spe_dist <- list(list(spe_index, dist_matrix))

sensitivity_list[4] <- spe_dist; names(sensitivity_list)[4] <- "integrated"

# Non-integrated distinctiveness (Fig. S2E) ----------
gaw<-compute_dist_matrix(trait, metric = "gower")
m.gaw<-as.matrix(gaw)
d.gaw<-as.data.frame(m.gaw); dist_matrix <- d.gaw

standard.Di.gow<-colSums(dist_matrix)/(nrow(dist_matrix)-1) #Mean of the computed distances for each species, to obtain the
#average distance to the other species based on all trait combinations
Int_Di<-as.data.frame(standard.Di.gow)

Int_Di$taxon<-rownames(Int_Di); colnames(Int_Di)[1]<-c("int_Di")

spe_index = mutate(Int_Di, quartile = as.factor(ntile(Int_Di$int_Di,4)),decile = as.factor(ntile(Int_Di$int_Di,10)),
                   combination = as.factor("Non-integrated"))

spe_dist <- list(list(spe_index, dist_matrix))

sensitivity_list[5] <- spe_dist; names(sensitivity_list)[5] <- "non-integrated"

dist_sens <- lapply(sensitivity_list, function(list){list[[1]]})

dist_sens <- do.call(rbind, dist_sens); rownames(dist_sens) <- NULL

# Density plots for all ways of computing distinctiveness ----

p <- dist_sens %>%
  ggplot(aes(x=int_Di, group = combination, fill = combination, color = combination)) + 
  xlim(0.1,0.4) + scale_color_aaas()+ scale_fill_aaas()+
  geom_density(alpha = 0.1) + labs(x = "Functional distinctiveness (Di)", y = "Density") +
  #ggtitle("Functional distinctiveness distribution") +
  theme_minimal() + 
  facet_wrap(~combination,ncol=2,scales = "free_x")+
  theme(legend.background = element_rect(fill = "transparent", color = NA), legend.position = "none", legend.title = element_blank(),
        text = element_text(size = 20), axis.title.y = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank()); p

#ggsave(file="sensitivity.svg", plot=p, width=10, height=8)

# 2.2: Calculate the community trait space; PCOA of functional distances ##########
dist_matrix<- read.csv("./Data/dist_matrix_ovrll.txt",header=T,dec=".",sep="\t", check.names = FALSE)
spe_index <- read.csv("./Data/spe_index.txt",header=T,dec=".",sep="\t", check.names = FALSE)
status<-read.csv("./Data/Raw_data/sp_status.txt",header=T,dec=".",sep="\t", check.names = FALSE)

gaw.pco<-cmdscale(dist_matrix, eig = T)

var<-round(gaw.pco$eig/sum(gaw.pco$eig)*100,1)

PCOA1<-data.frame(gaw.pco$points[,1])[,1]; PCOA2<-data.frame(gaw.pco$points[,2])[,1]; sum(gaw.pco$values$Cum_corr_eig)
PCOA<-data.frame(PCOA1,PCOA2)

PCOA_IDi<-cbind(spe_index,PCOA) #Bind each species with its position in the PCOA axis

row <- match(rownames(trait),spe_index$taxon) #The number of species needs to be the same in the distance dataframe and in the trait database
row <- which(!is.na(row));trait<-trait[row,]

gaw.env <- envfit(gaw.pco, trait, permutations = 999, na.rm = T) # With na.rm = TRUE it removes the rows containing NA values

# Plot to see the different loadings of modalities related with the position of species
funct<-spe_index
funct$group[funct$decile==1]<-"common"
funct$group[funct$decile==10]<-"distinct"
funct$group[funct$decile!=10 & funct$decile!=1]<-"intermediate"
grp=as.data.frame(funct$group)
colnames(grp)="group"

df1<-PCOA
df2<-data.frame(gaw.env$vectors$arrows)
df2<-rbind(df2,data.frame(gaw.env$factors$centroids[1:nrow(gaw.env$factors$centroids),]))

# Fig. 3C: Community trait space with the most functionally common/distinct species and NIS ----
tsp <- ggplot(aes(PCOA1,PCOA2,color = funct$group), data=df1)+
  geom_point(size=3)+
  scale_color_manual(values=c("#00AFBB","grey","#FA4848"), limits=c("common","intermediate","distinct"))+
  geom_text(data = PCOA_IDi[which(!is.na(match(PCOA_IDi$taxon, status[which(status$status == "Non-indigenous"),"taxon"]))),],
            aes(label=taxon),
            nudge_x = 0.03, nudge_y = 0.03,
            color = "darkgreen", size = 5)+
  geom_point(data = PCOA_IDi[which(!is.na(match(PCOA_IDi$taxon, status[which(status$status == "Non-indigenous"),"taxon"]))),], colour = "green", size = 3)+
  labs(x=paste("PCOA1 (",var[1],"%",")", sep = ""),y=paste("PCOA2 (",var[2],"%",")", sep = ""),color="Functional distinctiveness group")+
  geom_hline(yintercept=0,linetype=3,size=1)+ 
  geom_vline(xintercept=0,linetype=3,size=1)+
  theme_bw(15)+
  theme(panel.grid=element_blank(), legend.position="none", text = element_text(size = 20)); tsp

# Fig. 3D: Biplot from the community trait space ----
biplot <-ggplot() +
  geom_segment(data = df2,aes(x = 0, y = 0, xend = Dim1, yend = Dim2), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),type = "closed"),
               linetype=1, size=0.6,colour = "blue"
  )+
  geom_point(data = PCOA_IDi[which(!is.na(match(PCOA_IDi$taxon, 
                                                status[which(status$status == "Non-indigenous"),"taxon"]))),],
             mapping = aes(x = PCOA1, y = PCOA2), colour = "green", size = 5) +
  geom_text_repel(data = df2,aes(Dim1,Dim2,label=row.names(df2)),max.overlaps=27,size=4)+
  geom_text_repel(data = PCOA_IDi[which(!is.na(match(PCOA_IDi$taxon, status[which(status$status == "Non-indigenous"),"taxon"]))),],
                  aes(x = PCOA1, y = PCOA2, label=taxon),
                  color = "darkgreen", size = 5)+
  labs(x=paste("PCOA1 (",var[1],"%",")", sep = ""),y=paste("PCOA2 (",var[2],"%",")", sep = ""))+
  geom_hline(yintercept=0,linetype=3,size=1) + 
  geom_vline(xintercept=0,linetype=3,size=1)+
  theme_bw(15)+
  theme(panel.grid=element_blank(), text = element_text(size = 20)); biplot


#Merge both panels together
crossPlot <- ggarrange(tsp, biplot, ncol=2, nrow=1, common.legend = TRUE, legend = "none"); crossPlot
#ggsave(file="Trait_space_biplot.svg", plot=crossPlot, width=20, height=10)

# 2.3: Calculate the effect that each trait has on the distance between species ----
#Load the data and the species traits again
trait<-read.csv("./Data/sp_traits.txt",header=T,dec=".",sep="\t", check.names = FALSE) 

#Here the file rel_abund is the same used when weighting species Di by AFDW
data<-read.csv("./Data/Raw_data/species_AFDW_2005_2020.txt",header=T,dec=".",sep="\t", check.names = FALSE)
rel_abund <- as.data.frame(data %>% 
                             group_by(taxon) %>% 
                             summarise_at(.vars = c("AFDW"), sum, na.rm = TRUE))

rel_abund$rel.abund <- (rel_abund$AFDW/sum(rel_abund$AFDW))*100 #Use the preferred abundance metric

rel_abund <- rel_abund[which(!is.na(match(rel_abund$taxon,rownames(trait)))),]
rownames(rel_abund) <- rel_abund$taxon; 
rel_abund$taxon <- NULL; rel_abund$AFDW <- NULL #Make sure to delete everything but the relative abundance
colnames(rel_abund) <- "rel.abund"

# traits to be removed
traits.effects<- trait; traits.effects[,1:ncol(traits.effects)] <- NULL
trait.removed<-matrix(0,length(TRAITS)+1,length(TRAITS))
colnames<-c(paste("X",sep ="", seq(1:length(TRAITS))))

trait.removed<-as.data.frame(trait.removed)

for (a in (1: ncol(trait.removed))){
  trait.removed[a]<-a
  trait.removed[a,a]<-NA
}

#Here we set a cell as NA, it will be the column that will be removed when computing the distance repeatedly
#Then, the columns of the trait database will be selected based on these "rows" values; thus one trait will be removed each time

TRAITS[1] <- "Size"; TRAITS[3] <- "Adult life span"

#IMPORTANT: This is a loop that runs in parallel, within several cores in your computer, in order to reduce the computing time. 
# The number of cores of your computer destined to run this piece of code are set in the beginning; in case that you would not like to run this in parallel
# then "%dopar%" needs to be changed to "%do%". 
# For more information about parallel tasks in R, check this link: https://nceas.github.io/oss-lessons/parallel-computing-in-r/parallel-computing-in-r.html

# Calculating all the distances removing a trait in each case
traits.effects <- foreach (m= 1:nrow(trait.removed), .combine = cbind,.packages = c("funrar","dplyr"), .errorhandling = "stop") %dopar% {
  traits.selected<-as.vector(t(trait.removed[m,]))
  traits.selected<-na.omit(traits.selected)
  traits.selected<-as.vector(traits.selected)
  
  data<-trait #select the traits of interest
  codex<-traits.selected
  name_traits <- TRAITS[codex]
  names <- as.vector(unlist(name_traits))
  data<-dplyr::select(data,all_of(names))
  
  # calculating Di with gawdis
  dist_matrix <- funct_distances(name_traits,data);
  
  dist_mat_w <- matrix(0, nrow(dist_matrix), ncol(dist_matrix)) 
  colnames(dist_mat_w) <- colnames(dist_matrix); rownames(dist_mat_w) <- rownames(dist_matrix)
  
  for (j in 1:ncol(dist_matrix)){ #Here we multiply the functional pairwise distance by the relative abundance of each species
    for(i in 1:nrow(dist_matrix)){
      temp <- dist_matrix[i,j]*rel_abund[which(rownames(rel_abund) == rownames(dist_matrix[i,])),]
      dist_mat_w[i,j] <- temp
    }
  }
  
  weight.dist <- data.frame(dist = colSums(dist_mat_w))
  
  teff <- data.frame(); 
  
  #Di weighted by species abundances
  for(i in 1:ncol(dist_mat_w)){
    sumDi <- weight.dist[which(rownames(weight.dist) == colnames(dist_mat_w)[i]),] #Here we multiply the functional distance 
    sumrelAB <- sum(rel_abund[-which(rownames(rel_abund) == colnames(dist_mat_w)[i]),])
    
    Di<-sumDi/sumrelAB
    
    temp <- data.frame(taxon = colnames(dist_mat_w)[i], Int_Di = Di)
    teff <- rbind(teff, temp)
  }
  
  teff <- as.data.frame(teff[,2]); rownames(teff) <- rownames(data)
  teff
}

colnames(traits.effects)[1:(ncol(traits.effects)-1)]<-names(TRAITS)#Here we obtain the mean distance of one species to the others when deleting one trait
colnames(traits.effects)[ncol(traits.effects)] <- "None" #Each column represents the trait that was removed, and the distances without that trait
#write.table(traits.effects,file="./Data/traits.effects.txt",sep="\t", row.names = TRUE)

# Fig. 4; Traits effects on functional distinctiveness -----------
data<-read.csv("./Data/traits.effects.txt",header=T,dec=".",sep="\t", check.names = FALSE) 
traits <- c(colnames(data)[1:ncol(data)-1])

#Obtain the difference between the full distinctiveness and distinctiveness without each trait
delta <- lapply(traits, function(trait){
  
  eff <- as.data.frame(data$None - data[,trait])
  colnames(eff) <- trait; rownames(eff) <- rownames(data)
  eff
  
})

delta <- do.call(cbind, delta); delta$taxon <- rownames(delta)
delta <- melt(delta)

#Obtain the percentage of variation 
perc_delta <- lapply(traits, function(trait){
  
  eff <- as.data.frame(((data$None - data[,trait])/data$None)*100)
  
  colnames(eff) <- trait
  rownames(eff) <- rownames(data)
  eff
  
})

perc_delta <- do.call(cbind, perc_delta); perc_delta$taxon <- rownames(perc_delta); sort(colMeans(perc_delta[,1:10]))

perc_delta <- melt(perc_delta)

#See which traits have a bigger effect on NIS distinctiveness
NIS_delta <- perc_delta[perc_delta$taxon %in% c("Marenzelleria","Mya arenaria","Streblospio benedicti",
                                                "Polydora cornuta","Potamopyrgus antipodarum"),]; #colMeans(NIS_delta[,1:10])
NIS_median <- lapply(traits, function(trait){
  
  eff <- as.data.frame(median(NIS_delta[which(NIS_delta$variable == trait), 3]))
  
  colnames(eff) <- trait
  #rownames(eff) <- rownames(NIS_delta)
  eff
  
})
NIS_median <- do.call(cbind, NIS_median); NIS_median <- melt(NIS_median)

NIS_mean <- lapply(traits, function(trait){
  
  eff <- as.data.frame(mean(NIS_delta[which(NIS_delta$variable == trait), 3]))
  
  colnames(eff) <- trait
  #rownames(eff) <- rownames(NIS_delta)
  eff
  
})
NIS_mean <- do.call(cbind, NIS_mean); NIS_mean <- melt(NIS_mean)

# Figure 4A; overall effect of each trait on distinctiveness
# NOTE: The final scale color was adjusted in Illustrator
all <-  ggplot(perc_delta, aes(x=as.factor(variable), y=value, fill = as.factor(variable))) + 
  labs(x = "", y = "Variation in distinctiveness (%)") +
  theme_minimal() + 
  theme(legend.position = "none", text = element_text(size = 15), 
        axis.text.x = element_text(vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12)) + 
  geom_hline(aes(yintercept = 0), size = 1,
             color = "firebrick", linetype = "dashed") +
  geom_boxplot(color="black", alpha=0.4) + 
  geom_point(data = NIS_median, aes(x=as.factor(variable), y=value), color = "black", shape = 18, size = 7) +
  coord_flip()+
  scale_fill_rickandmorty(); all


# Figure 4B; effect of each trait on NIS distinctiveness
NIS_eff <-  ggplot(NIS_delta, aes(x=as.factor(variable), y=value, fill = as.factor(variable))) +
  labs(x = "", y = "Variation in distinctiveness (%)") +
  theme_minimal() + 
  theme(legend.position = "none", text = element_text(size = 15), 
        axis.text.x = element_text(vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12)) + 
  geom_hline(aes(yintercept = 0), size = 1,
             color = "firebrick", linetype = "dashed") +
  geom_point(size = 7, aes(color = taxon)) +  
  scale_color_manual(values= c("firebrick1", "orange2", "chartreuse3", "royalblue3", "thistle3")) + 
  geom_point(data = NIS_mean, aes(x=as.factor(variable), y=value), color = "black", shape = 18, size = 7) +
  coord_flip(); NIS_eff

crossPlot <- ggarrange(all,NIS_eff,ncol=2, nrow=1, common.legend = F, legend = "none"); crossPlot

#ggsave(file="Trait_effect.svg", plot=crossPlot, width=20, height=8)




