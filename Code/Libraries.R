#####################!
# LIST OF LIBRARIES #!
#####################

listoflibrary<-c("data.table","vegan","utils","tidyr","stats","grid", "doParallel", "cowplot",
                 "ggsci","ggplot2","ggrepel","ggpubr","funrar","foreach","dplyr", "raster", "mgcViz", "gridExtra",
                 "sp", "rangeBuilder", "colorspace", "automap", "ape", "seegSDM", "memisc", "mgcv", "ggnewscale")

for (pkg in listoflibrary){
  if(!eval(bquote(require(.(pkg))))) {eval(bquote(install.packages(.(pkg))))}
  eval(bquote(library(.(pkg))))
}