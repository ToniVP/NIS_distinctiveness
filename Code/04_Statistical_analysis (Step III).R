###########################################################################!
#4: Statistical analysis: Multimodel approach and single NIS models #######!
###########################################################################!
remove(list=ls())
#Load libraries and functions
source("./Code/Functions.R")
source("./Code/Libraries.R")

#Set the number of cores for parallel tasks
numcores <- detectCores()
registerDoParallel(numcores -1)

# 4.1: Data loading and exploration ----
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

data <- data[which(!duplicated(data[,coln])),] #Remove possible duplicates due to the merge function
test <- data

#test for correlation between variables
variables <- test %>% dplyr::select(Shannon, bot_sal, bot_oxy, bot_T, nbsp, depth)
library(psych); #install.packages("psych")
pairs.panels(variables, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)

#Variance inflation factors
library(car)
vif <- lm(NIS_Di ~ depth + Shannon + nbsp +bot_sal + bot_oxy + bot_T, test); summary (vif)

vif_values <- vif(vif)
barplot(vif_values, main = "VIF Values", horiz = TRUE, col = "steelblue")
abline(v = 5, lwd = 3, lty = 2)


# 4.2: Multi-model approach ----

#With this function we will extract an equal subset of observations for each species and fit a model each time
#IMPORTANT: the model that you want to fit (GLM, GAM, RF...) NEEDS to be specified in this function, as it will run such model for each iteration

resample_species <- function (NIS, num, data) {
  Di_unified_res <- data.frame()
  
  for (i in 1:length(NIS)){
    
    temp <- data[which(data$species == NIS[i]),]
    temp <- temp[sample(nrow(temp), num),]
    Di_unified_res <- rbind(Di_unified_res, temp)
    
  }
  
  test <- Di_unified_res
  
  system.time(fitTL<- mgcv::gam(NIS_Di ~ 
                            s(Shannon, k = 3) +
                            s(nbsp, k = 3) +
                            s(depth, k = 3) + 
                            s(bot_T, k = 3) + 
                            s(bot_sal, k = 3) + 
                            s(bot_oxy, k = 3)+
                            s(species, bs = 're') +
                            s(station, bs = 're'),
                          family=betar(link="logit"),
                          niterPQL=20,
                          data = test))
  return(fitTL)
} 
#"NIS" accounts for a vector with the species names, as they appear in your database
#"num" is the number of observations to randomly select

#Now we will resample 40 observations for each NIS, and fit a model each time with "species" as random factor, to detect potential 
#overall patterns on NIS distinctiveness, not driven by the most abundant species
num = 40; NIS = NIS; data = data
mod.list <- list();
q = 1

system.time (repeat { 
  tryCatch({
    print(q)
    model <- resample_species(NIS, num, data)
    
    mod.list[[q]] <- model; names(mod.list)[q] <- paste("model", q, sep = "_")
    
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n"); mod.list[[q]] = NA; summ.list[[q]] = NA}) 
  q = q+1
  if (q >400){break}
}) #Make as many models as you would like to, using the equalized database and the species and location as random factor, then put all of them into a list

mod.list <- mod.list[lengths(mod.list)!= 0] #If the model cannot be shown then remove it

#saveRDS(mod.list, file="./Data/Models/multi_model_list")
mod.list<- readRDS("./Data/Models/multi_model_list");

#Extract several parameters for each model in order to compare all of them and extract possible overall patterns
# R squared
models.r.sq <- sapply(mod.list, function(i)summary(i)$r.sq)

# Deviance explained
models.deviance <- sapply(mod.list, function(i) summary(i)$dev.expl)

# Pvalue for smooth terms and AIC
summ <- summary(mod.list[[1]])
variables <- c(names(summ$chi.sq)); variables <- gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", variables, perl=T)
models.pvalue <- data.frame()

for (i in 1:length(mod.list)){
  print(i)
  model <- mod.list[[i]]
  p.val <- t(data.frame(p.value = summary(model)$s.pv)); colnames(p.val) <- variables; rownames(p.val)<- names(mod.list[i])
  aic <- AIC(model); p.val <- cbind(p.val, AIC = aic)
  models.pvalue <- rbind(models.pvalue,p.val)
  
}

#Bind all the model descriptors in one database
df.model.eval <- cbind(r.sq=models.r.sq, dev.explain=models.deviance, models.pvalue)

#Compute how many models account for significant smooth terms per each factor
tot_pval <- data.frame()

for (i in 1:length(variables)){
  var = variables[i]
  
  sign <- sum(df.model.eval[,var] < 0.05 & df.model.eval[,var] >= 0.01); non_sign <- sum(df.model.eval[,var] >= 0.1)
  sign0.1 <- sum(df.model.eval[,var] >= 0.05 & df.model.eval[,var] < 0.1)
  sign0.01 <- sum(df.model.eval[,var] >= 0.001 & df.model.eval[,var] < 0.01)
  sign0.001 <- sum(df.model.eval[,var] < 0.001)
  
  temp <- data.frame(`<0.001`=sign0.001, `<0.01`=sign0.01, `<0.05` = sign, `<0.1` = sign0.1, 
                     `Non-significant` = non_sign, check.names= FALSE); temp$variable <- var
  tot_pval <- rbind(tot_pval, temp)
}
perc_pval <- tot_pval
tot_pval <- melt(tot_pval); colnames(tot_pval)[1:2] <- c("variable", "significance")
rownames(perc_pval) <- perc_pval$variable; perc_pval$variable <- NULL
perc_pval <- as.data.frame((perc_pval/400)*100) #Obtain the percentage of the total models where each variable was significant

# Fig. 5C-D; Multi-model results and overall trends ----

# Fig. 5C; Plot how many models show effects from the variables
#remove the random factors for the plot
tot_pval_noR <- tot_pval[which(tot_pval$variable != "species" & tot_pval$variable != "station"),]

pval <- ggplot(tot_pval_noR, aes(fill=significance, y=value, x=variable)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_manual(values = c("mistyrose1", "lightcoral", "darkred", "orange", "royalblue4"))+
  theme_minimal() + geom_hline(aes(yintercept = 0.5), size = 2,
                               color = "limegreen", linetype = "dashed") +
  #facet_wrap( ~ variable, nrow = 2, scales = "free") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()); pval

#ggsave(file="Multimodel_pval", plot=pval, width = 11, height = 7)

#Fig. 5D; Smooth lines from all models to see general trends
lines <- data.frame()

# Converting GAM objects to a gamViz objects, we end up with a big dataframe of all the points that define our smooth lines
for (i in 1:length(mod.list)){
  
  print(i); model <- mod.list[[i]]
  #with this function we extract all the points defining the smooth lines per each GAM, then save the generated plot as a list
  viz_fit <- mgcViz::getViz(model); trt_fit <- plot(viz_fit, allTerms = T) + l_fitLine() 
  
  for (j in 1:length(trt_fit[["plots"]])){ #We do it for every factor that we have in our models
    plot <- as.data.frame(trt_fit[["plots"]][[j]]$data$fit %>% mutate( fit = paste("Model", i, sep = "_")))
    plot$variable <- trt_fit[["plots"]][[j]][["ggObj"]][["labels"]][[1]]
    lines <- bind_rows(lines, plot)
  }
  
}

lines <- lines[which(lines$variable != "Gaussian quantiles"),]

#plot all the smooths from all the models and variables per each factor
lines_plot <- ggplot(data = lines, aes(x = x, y = y)) +
  geom_line(alpha = 1, color = "royalblue4") +
  geom_smooth(method = mgcv::gam, size = 2, color = "firebrick") +
  #labs(x = var, y = paste("s(", var, ")", sep = "")) + 
  theme_minimal(base_size = 12) + 
  facet_wrap( ~ variable, nrow = 2, scales = ("free_x")) + 
  theme(axis.title.x=element_blank(), axis.title.y = element_blank()); lines_plot

# Combine the pvalues plot with the smooth lines one
library(cowplot)
comb_plot <- plot_grid(pval, lines_plot, ncol = 1); comb_plot

#ggsave(file="Multimodel_results.svg", plot=comb_plot, width=11, height=14)

# 4.3: Single-species models (Appendix 2) ----
# Now single GAMs will be fit for each NIS, to see individual drivers of distinctiveness; "data" is the whole dataframe with environmental variables
NIS <- c("Marenzelleria", "Mya.arenaria", "Potamopyrgus.antipodarum","Streblospio.benedicti", "Polydora.cornuta")
NIS_mod <- list()

for (i in 1:length(NIS)){
  print(NIS[i])
  test <- data[which(data$species == NIS[i]),]
  system.time(fitTL<- gam(NIS_Di ~ s(Shannon, k = 3) +
                            s(richness, k = 3) +
                            s(depth, k = 3) + 
                            s(bot_T, k = 3) + 
                            s(bot_sal, k = 3) +
                            s(bot_oxy, k = 3) +
                            s(station, bs = 're', k = 3),
                          na.action = "na.omit",
                          family = betar(link = "logit"), 
                          data = test))
  
  NIS_mod[[i]] <- fitTL; names(NIS_mod)[i] <- NIS[i]
}

#saveRDS(NIS_mod, "./Data/Models/NIS_single_models")

NIS_mod<- readRDS("./Data/Models/NIS_single_models") #transformed one
fitTL = NIS_mod[[1]] #Check models individually

# Inspect model summary
summary(fitTL)

par(mfrow=c(4,2),mar=c(4,4,1,1))
plot(fitTL,shade=T,shade.col="grey",res=T,rug=F,pch=20)

# Check diagnostics
par(mfrow=c(2,2))
gam.check(fitTL, old.style = T);#gam.check(fitTL.1$gam)

par(mfrow=c(1,1),mar=c(5,4,2,2))

# Fig. S5; Single species distinctiveness trends ----

#Combine the smooth lines for the different models for each factor
mod.list <- readRDS("./Data/Models/NIS_single_models") #important that the order of the models in this list is the same as the NIS names in the vector
NIS <- c("Marenzelleria", "Mya.arenaria", "Potamopyrgus.antipodarum","Streblospio.benedicti", "Polydora.cornuta")#put the names of your NIS here

#Plot all the smooth lines from all models into just one plot to see general trends
lines <- data.frame(); resid <- data.frame()

# Converting GAM objects to a gamViz objects, we end up with a big dataframe of all the points that define our smooth lines
for (i in 1:length(mod.list)){
  
  print(i); model <- mod.list[[i]]
  #with this function we extract all the points defining the smooth lines per each GAM, then save the generated plot as a list
  viz_fit <- mgcViz::getViz(model); trt_fit <- plot(viz_fit, select = c(1:6), res = T) + l_fitLine() + l_points() + l_ciLine() + l_rug()
  
  for (j in 1:length(trt_fit[["plots"]])){ #We do it for every factor that we have in our models
    plot <- as.data.frame(trt_fit[["plots"]][[j]]$data$fit %>% mutate( fit = paste(NIS[i])))
    res <- as.data.frame(trt_fit[["plots"]][[j]]$data$res %>% mutate( fit = paste(NIS[i])))
    plot$variable <- trt_fit[["plots"]][[j]][["ggObj"]][["labels"]][[1]]
    res$variable <- trt_fit[["plots"]][[j]][["ggObj"]][["labels"]][[1]]
    resid <- bind_rows(resid, res)
    lines <- bind_rows(lines, plot)
  }
  
}

# Plotting all the species lines together for each variable
var <- lines; var.resid <- resid

sing_NIS <- ggplot() + 
  labs(x = "", y = "s(variable)") +
  ylim(-1,1)+
  geom_point(data = var.resid, aes(x,y, group = fit, colour = fit, shape = fit), alpha = 0.1, size = 0.5) +
  scale_color_viridis(discrete = T, name="Species", labels=levels(as.factor(var$fit))) +
  scale_fill_viridis(discrete = T, name="Species", labels=levels(as.factor(var$fit))) +
  geom_line(data = var, aes(x, y, group = fit, colour = fit),  size = 0.7) + 
  theme_minimal() + 
  geom_ribbon(data = var, aes(x, ymin = y-se, ymax = y+se, fill = fit, color = fit), size = 0.5, alpha = 0.3, linetype = 2, show.legend = F) +
  #guides(shape = "none") +
  scale_shape_manual(values= c(15, 16, 9, 17, 4), name="Species", labels=levels(as.factor(var.resid$fit))) +
  facet_wrap( ~ variable, nrow = 2, scales = ("free_x")) + 
  theme(legend.background = element_rect(fill = "transparent", color = NA), legend.key = element_rect(fill = "transparent", color = NA),
        text = element_text(size = 11), legend.position = "bottom") +
  labs(color  = var$fit, linetype = var$fit, shape = var$fit); sing_NIS

#ggsave(file="GAM_sing_NIS.svg", plot=sing_NIS, width=10, height=8)

legend <- cowplot::get_legend(p)
grid.draw(legend)
#ggsave(file="legendGAM.svg", plot=legend, width=10, height=8)
