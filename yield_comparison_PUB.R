# Comparison of GAEZ 2015 crop yields to the Global dataset of historical yields (GDHY) v1.3
# GDHY data downloaded from https://doi.pangaea.de/10.1594/PANGAEA.909132 

# Danielle S Grogan
# 2021

##################################################################################################################
### LIBRARIES AND SOURCE FILES ###
##################################################################################################################
library(raster)
library(rgdal)
library(rgeos)
library(ggplot2)
library(rasterVis)
library(RColorBrewer)
library(gridExtra)
library(stringr)
library(tidyr)
library(ggthemes)

##################################################################################################################
# function to converty gdhy to same latlong system as gaez (gdhy has x from 0 to 360; gaez is -180 to 180)
set_latlong = function(gdhy.data){
  gdhy.part1 = crop(gdhy.data, extent(0, 180, -90, 90))
  gdhy.part2 = crop(gdhy.data, extent(180, 360, -90, 90))
  extent(gdhy.part2)[1] = -180
  extent(gdhy.part2)[2] = 0
  
  gdhy.new = raster::merge(gdhy.part2, gdhy.part1)
}

##################################################################################################################
# make color breaks around 0 for use with levelplot()
color_scale_0<-function(max.val, min.val, zero.gap){
  
  col.scale = unique(c(seq(min.val-3,       min.val-1*(min.val/4),  length=100), 
                       seq(min.val-1*(min.val/4), min.val-2*(min.val/4),  length=100), 
                       seq(min.val-2*(min.val/4), min.val-3*(min.val/4),  length=100), 
                       seq(min.val-3*(min.val/4), 0 - zero.gap,   length=100), 
                       seq(0 - zero.gap,  0 + zero.gap,   length=100),
                       seq(0 + zero.gap,  max.val-3*(max.val/4),  length=100),
                       seq(max.val-3*(max.val/4), max.val-2*(max.val/4),  length=100),
                       seq(max.val-2*(max.val/4), max.val-1*(max.val/4),  length=100),
                       seq(max.val-1*(max.val/4), max.val+3,            length=100)))
  # the min-3 and max+3 ensures that the edges of the color labels are not cut off
}
################################################################################################################
# File paths

# gdhy.path = "path to folder containing GDHY data"
# gaez.harv.path = "path to folder containing GAEZ+_2015 crop harvest area data"  # download from https://doi.org/10.7910/DVN/XGGJAV
# gaez.prod.path = "path to folder containing GAEZ+_2015 crop production data"    # download from https://doi.org/10.7910/DVN/KJFUO1
# gaez.yield.path= "path to folder containing GAEZ+_2015 crop yield data"         # download from https://doi.org/10.7910/DVN/XGGJAV
 


################################################################################################################
# Grid cell yield comparison

# List of crops to evaluate
crops = c("maize", "rice", "soybean", "wheat")
# UNITS: GDHY and GAEZ yield data are both in t/ha

for(i in 1:length(crops)){
  # gdhy yield data
  gdhy = set_latlong(raster(file.path(gdhy.path, crops[i], "yield_2015.nc4")))
  
  # calculate GAEZ yields on 0.5 degree resolution harvested area and production (accounts for differences yield across in areas)
  gaez.harv = raster::aggregate(raster(paste(gaez.harv.path, "/GAEZAct2015_HarvArea_",   str_to_title(crops[i]), "_Total.tif", sep="")), fact=6, fun=sum)
  gaez.prod = raster::aggregate(raster(paste(gaez.prod.path, "/GAEZAct2015_Production_", str_to_title(crops[i]), "_Total.tif", sep="")), fact=6, fun=sum)
  gaez = gaez.prod/gaez.harv
  
  ### Global grid cell comparison
  gdhy.list = values(gdhy)
  gaez.list = values(gaez)
  both.list = as.data.frame(cbind(gdhy.list, gaez.list))
  
  colnames(both.list) = c("GDHY_Yield", "GAEZ_Yield")
  both.list = both.list[complete.cases(both.list[ , 1:2]),] # remove rows where both values are NA
  both.list$Crop = str_to_title(crops[i])
  
  if(i == 1){
    out.df = as.data.frame(both.list)
  }else{
    out.df = rbind(out.df, both.list)
  }
  
  print(i)
}
write.csv(out.df, "results/yield_validation/GAEZ_and_GDHY_yield_by_grid.csv")


##################################################################################################################
### Hexbin plots ###
##################################################################################################################

out.df = read.csv("results/yield_validation/GAEZ_and_GDHY_yield_by_grid.csv")

p=ggplot(data = out.df, aes(x = GDHY_Yield, y = GAEZ_Yield))+
  geom_hex(bins = 30, aes(fill=log(..count..)))+
  scale_fill_continuous(high = "#010B12", low = "#9FD3FA")+
  geom_smooth(method = lm, colour="black", size=0.5)+
  xlab("GDHY Crop Yield (t/ha)")+
  ylab("GAEZ+ 2015 Crop Yield (t/ha)")+
  theme_bw(base_size = 12)+ 
  facet_wrap(~Crop, scales = "free")

ggsave(paste("results/figures/yield_validation/GAEZ_vs_GDHY_hexbin_v2.png", sep=""),
       plot = p,
       width = 6, height = 4.5, units = "in",
       dpi = 300) 


# Linear regression:
maize = subset(out.df, out.df$Crop == "Maize") 
summary(lm(maize$GAEZ_Yield ~ maize$GDHY_Yield))

rice = subset(out.df, out.df$Crop == "Rice") 
summary(lm(rice$GAEZ_Yield ~ rice$GDHY_Yield))

soybean = subset(out.df, out.df$Crop == "Soybean") 
summary(lm(soybean$GAEZ_Yield ~ soybean$GDHY_Yield))

wheat = subset(out.df, out.df$Crop == "Wheat") 
summary(lm(wheat$GAEZ_Yield ~ wheat$GDHY_Yield))

out.table = rbind(maize.total, rice.total, soybean.total, wheat.total)
write.csv(out.table, "results/yield_validation/GDHY_and_GAEZ_yield_totals.csv")

##################################################################################################################
### Maps ###
##################################################################################################################
# for country boarders, use the GAUL_2012 national units layer provided here: https://data.jrc.ec.europa.eu/dataset/jrc-10112-10004
# country.borders = readOGR("path to GAUL_2012 shapefile", layer = "gaul0_asap") 

zmax = c(43, 24, 16, 19)

# list all gdhy file paths
gdhy.list = unlist(lapply(seq(1:4), FUN = function(i) file.path(gdhy.path, crops[i], "yield_2015.nc4")))

# load gdhy yield data
gdhy = do.call(stack,
               lapply(gdhy.list,
                      raster::brick)
)

# set to same lat/long system as gaez data
gdhy.yield = set_latlong(gdhy)


# list all gaez file paths
gaez.list = unlist(lapply(seq(1:4), 
                          FUN = function(i) paste(gaez.yield.path, "GAEZAct2015_Yield_", str_to_title(crops[i]), "_Mean.tif", sep="")))

# load gdhy yield data
gaez.yield = do.call(stack,
                     lapply(gaez.list,
                            raster::brick)
)


# Keep the high resolution of gaez yields, so dissaggregate the gdhy to same reslution (values will be the same)
gdhy.yield = raster::disaggregate(gdhy.yield, fact=6)

for(i in 1:length(crops)){
  
  # subset to a single crop
  gdhy.c = subset(gdhy.yield, i)
  gaez.c = subset(gaez.yield, i)  
  c.stack = stack(gdhy.c, gaez.c)  
  
  # saturate colors at 90% of max
  mx.90 = 0.90 * (max(as.matrix(c.stack), na.rm=T))
  c.stack[c.stack > mx.90] = mx.90
  
  # set names and map labels
  names(c.stack) = c("GDHY_Yield", "GAEZ_2015_Yield")
  labs = c(
    GDHY_Yield      = paste("GDHY", str_to_title(crops[i]), "Yield"),
    GAEZ_2015_Yield = paste("GAEZ+ 2015", str_to_title(crops[i]), "Yield")
  )
  
  # maps
  m.plot = 
    gplot(c.stack) + 
    geom_raster(aes(fill=value)) +
    borders(colour="grey50", size=0.3) + 
    coord_quickmap()+
    ylim(-60,90)+
    scale_fill_distiller(na.value="white", palette="YlGnBu", direction = 1, name="Yield (t/ha)") +
    facet_wrap(~variable, labeller = labeller(variable = labs)) +
    theme_map(base_size = 20)
  
  # save to file
  ggsave(paste("results/figures/yield_validation/GAEZ_vs_GDHY_", crops[i], "_maps.png", sep=""),
         plot = m.plot,
         width = 16, height = 4, units = "in",
         dpi = 300) 
}
