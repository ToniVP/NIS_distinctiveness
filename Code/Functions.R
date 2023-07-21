#########################!
## NEEDED FUNCTIONS #####!
#########################

# Function to select common species from the trait and the occurrences database---- 

get_common_species = function(x, y, num) { 
  
  if(class(x) == "data.frame"){colnames(x)<- gsub("\\.", " ", colnames(x)); mat_columns = colnames(x)[num:ncol(x)]}
  else {mat_columns = x} #the abundance matrix or a vector containing names
  
  trait_species = rownames(y) #the trait dataset
  
  # Extract matching species
  common_species = match(trait_species, mat_columns)
  rows = which(!is.na(common_species))
  
  return(y[rows,]) 
  
} 
# x as the abundance matrix fo the year/station that we want, as there are the selected species
# y as the raw trait matrix, to select species present in the abundance matrix
# num as the column number where the first species is

# Function to standardize (only if needed) ----

standr = function(x){(x-min(x))/(max(x)-min(x))} 

# Obtain the mean of each element from a list of matrices (needed to compute the overall functional distances matrix) ----

mean_matrix = function(x){
  y <- array(unlist(x) , c(dim(x[[1]]),length(x)))
  y <- apply(y, 1:2 , mean )
  colnames(y)<-colnames(x[[1]])
  rownames(y)<-rownames(x[[1]])
  return(y)
} 

# Function to calculate the overall pairwise matrix of functional distances for the regional species pool in an integrated way----

funct_distances = function(TRAITS, trait){
  
  dist.matr <- list()
  
  # N is the total number of traits (columns)
  N<-length(TRAITS)
  
  # Number of columns corresponding to traits
  listk<-seq(1,N,1)
  
  if (length(TRAITS) >= 4){#start condition for more than 4 trait categories
    k=4
    for (k in 4:N){ #start loop combination of traits
      
      #Here are determined all possible 4 trait combinations for obtaining the average distance
      #between species, here are expressed only the position that are occupying the traits in the data frame
      
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
        #gow<-compute_dist_matrix(data, metric = "euclidean") #euclidean if all your traits are numeric
        m.gow<-as.matrix(gow)
        d.gow<-as.data.frame(m.gow)
        #standard.gow<-standr(d.gow) # distances are standardized to fit between 1 and 0 (if other distance metric is used)
        dist.matr <- c(dist.matr,list(d.gow))
      }
      
      print (k)
    }
    
    dist_matrix <- as.data.frame(mean_matrix(dist.matr))
    
  }#end condition for more than 4 trait categories
  
  else{#start condition for less than 3 categories
    
    #Calculation of distinctiveness using Gower's distance
    gow<-compute_dist_matrix(trait, metric = "gower") # a distance matrix is calculated for each combination of traits
    #gow<-compute_dist_matrix(data, metric = "euclidean") #euclidean if all your traits are numeric
    m.gow<-as.matrix(gow)
    dist_matrix<-as.data.frame(m.gow)
    #dist_matrix<-standr(d.gow) # distances are standardized to fit between 1 and 0
  } #end condition for less than 3 categories
  
  return(dist_matrix)
} 
#"TRAITS" is a list of names for all your traits and modalities, ALL modalities should be inside this list.

#Here "trait" is a species x trait information matrix, previously cleaned with the traits of interest selected.

# Function to calculate distinctiveness between different species groups (NIS vs natives) at local communities----

mean_diss = function(dist_matrix, data, status, num){
  
  if(ncol(data) > 1){
    if (all(is.na(match(colnames(data), rownames(status[which(status$status == "Non-indigenous"),]))))){NIS <- c()} else {
      NIS <- c(colnames(data)[which(!is.na(match(colnames(data), rownames(status[which(status$status == "Non-indigenous"),]))))])}
    
    Allsp <- c(colnames(data)[num:ncol(data)]); Allsp <- gsub("\\.", " ", Allsp)
    NAT <- setdiff(Allsp,NIS)
    
    ratio_sp <- (length(NIS)/length(Allsp))*100
    
    #mean distance of NIS to other native species per each location/cell
    
    NIS_Di <- foreach(i = 1:length(NIS), .combine = cbind, .errorhandling = "pass") %do% {
      dis <- c()
      sp <- Allsp[!Allsp %in% NIS[i]]
      for(j in 1:length(sp)) {val <- dist_matrix[NIS[i],sp[j]]; names(val) <- NIS[i]; dis <- c(dis,val)}
      disW <- sum(dis); DiW <- disW/length(sp);  dist <- as.data.frame(DiW); colnames(dist) <- NIS[i]; dist}
    
    
    if(length(NAT) == 0){Di_NAT <- NA} else{ 
      Di_NAT <- c(); for (j in 1:length(NAT)){ #Here we multiply the functional pairwise distance by the relative abundance of each species
        dis <- c()
        sp <- Allsp[!Allsp %in% NAT[j]]
        for(i in 1:length(sp)){
          val <- dist_matrix[sp[i],NAT[j]]; names(val) <- NAT[j]; dis <- c(dis,val)
        }
        disW <- sum(dis); DiW <- disW/length(sp);  dist <- as.data.frame(DiW); colnames(dist) <- NAT[j] 
        Di_NAT <- c(Di_NAT,dist)
      }
      Di_NAT <- unlist(Di_NAT); Di_NAT <- mean(Di_NAT)
      #val <- Int_Di[which(!is.na(match(Int_Di$taxon,NAT))),]
      #Di_NAT <- mean(Di_NAT)
    }
    
    
    Di_ovrll<-c()
    for (j in 1:length(Allsp)){ #Here we multiply the functional pairwise distance by the relative abundance of each species
      dis <- c()
      sp <- Allsp[!Allsp %in% Allsp[j]]
      for(i in 1:length(sp)){
        val <- dist_matrix[sp[i],Allsp[j]]; names(val) <- Allsp[j]; dis <- c(dis,val)
      }
      disW <- sum(dis); DiW <- disW/length(sp);  dist <- as.data.frame(DiW); colnames(dist) <- Allsp[j] 
      Di_ovrll <- c(Di_ovrll,dist)
    }
    Di_ovrll <- unlist(Di_ovrll); Di_ovrll <- mean(Di_ovrll)
  } else {NIS_Di = NA; Di_NAT = NA; Di_ovrll = NA; ratio_sp = NA}
  
  return(list(NIS_Di, Di_NAT, Di_ovrll, ratio_sp))
}   
#"dist_matrix" accounts for the overall pairwise matrix of functional distances
#"data" is the species occurrences/abundance data for a certain location/cell where we have different sampling locations
#"status" is a list for all the species we have in TOTAL where their status (NIS or native) is defined
#"num" accounts for the column number that corresponds to the first species in your sites x abundances matrix

# Function to calculate distinctiveness between different species groups weighted by relative abundance at local communities----

mean_diss_rel_ab = function(dist_matrix, data, status, num){
  
  if(ncol(data) > 1){
    abund <- sum(data)
    
    rel_ab <- data.frame(t((sapply(data, function (x) (x/abund)*100)))); rownames(rel_ab)<-NULL
    colnames(rel_ab) <- colnames(data)
    
    if (all(is.na(match(colnames(data),rownames(status[which(status$status == "Non-indigenous"),]))))){
      NIS <- c()} else {NIS <- c(colnames(data)[which(!is.na(match(colnames(data), 
                                                                   rownames(status[which(status$status == "Non-indigenous"),]))))])}
    
    Allsp <- c(colnames(data)[num:ncol(data)]); Allsp <- gsub("\\.", " ", Allsp)
    NAT <- setdiff(Allsp,NIS)
    
    ratio_sp <- (length(NIS)/length(Allsp))*100
    
    #mean distance of NIS to other native species per each location/cell
    if(length(NIS) == 0){NIS_Di <- NA} else {NIS_Di <- foreach(i = 1:length(NIS), .combine = cbind, .errorhandling = "pass") %do% {
      dis <- c()
      abj <- c()
      sp <- Allsp[!Allsp %in% NIS[i]]
      for(j in 1:length(sp)) {val <- dist_matrix[sp[j],NIS[i]]; names(val) <- sp[i]; dis <- c(dis,val)
      ab <- rel_ab[,which(colnames(rel_ab) == sp[j])]; abj<-c(abj,ab)
      }
      #dis <- mean(dis)
      disW <- sum(dis*abj); DiW <- disW/sum(abj);  dist <- as.data.frame(DiW); colnames(dist) <- NIS[i]; NIS_relab <- rel_ab[NIS[i]]
      name = paste(NIS[i],"relab", sep="_"); colnames(NIS_relab) <- name
      dist <- cbind(dist, NIS_relab);
      dist}
    if(length(NIS) > 1){NIS_Di$mean_NIS_Di <- rowMeans(NIS_Di[,c(NIS)])}else{NIS_Di$mean_NIS_Di <- NIS_Di[,c(NIS)]}
    
    }
    
    if(length(NAT) == 0){Di_NAT <- NA} else{ Di_NAT <- c(); for (j in 1:length(NAT)){ #Here we multiply the functional pairwise distance by the relative abundance of each species
      dis <- c()
      abj <- c()
      sp <- Allsp[!Allsp %in% NAT[j]]
      for(i in 1:length(sp)){
        val <- dist_matrix[sp[i],NAT[j]]; names(val) <- NAT[j]; dis <- c(dis,val)
        ab <- rel_ab[,which(colnames(rel_ab) == sp[i])]; abj<-c(abj,ab)
      }
      disW <- sum(dis*abj); DiW <- disW/sum(abj);  dist <- as.data.frame(DiW); colnames(dist) <- NAT[j] 
      Di_NAT <- c(Di_NAT,dist)
    }
    Di_NAT <- unlist(Di_NAT); Di_NAT <- mean(Di_NAT)
    #val <- Int_Di[which(!is.na(match(Int_Di$taxon,NAT))),]
    #Di_NAT <- mean(Di_NAT)
    }
    
    Di_ovrll<-c()
    for (j in 1:length(Allsp)){ #Here we multiply the functional pairwise distance by the relative abundance of each species
      dis <- c()
      abj <- c()
      sp <- Allsp[!Allsp %in% Allsp[j]]
      for(i in 1:length(sp)){
        val <- dist_matrix[sp[i],Allsp[j]]; names(val) <- Allsp[j]; dis <- c(dis,val)
        ab <- rel_ab[,which(colnames(rel_ab) == sp[i])]; abj<-c(abj,ab)
      }
      disW <- sum(dis*abj); DiW <- disW/sum(abj);  dist <- as.data.frame(DiW); colnames(dist) <- Allsp[j] 
      Di_ovrll <- c(Di_ovrll,dist)
    }
    Di_ovrll <- unlist(Di_ovrll); Di_ovrll <- mean(Di_ovrll)
    # val <- Int_Di[which(!is.na(match(Int_Di$taxon,Allsp))),]
    # Di_ovrll <- mean(val$int_Di)
    
    
  } else {NIS_Di = NA; Di_NAT = NA; Di_ovrll = NA; ratio_sp = NA}
  
  
  return(list(NIS_Di, Di_NAT, Di_ovrll, ratio_sp))
}
#"dist_matrix" accounts for the overall pairwise matrix of functional distances
#"data" is the species occurrences/abundance data for a certain location/cell where we have different sampling locations
#"status" is a list for all the species we have in TOTAL where their status (NIS or native) is defined
#"num" accounts for the column number that corresponds to the first species in your sites x abundances matrix




