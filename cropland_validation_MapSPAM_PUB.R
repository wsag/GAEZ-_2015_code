# Cropland validation
# Compare the presence/absence of cropland in GAEZ+_2015 to MapSPAM data

# Danielle S Grogan
# 2021-01

##################################################################################################################
### LIBRARIES AND SOURCE FILES ###
library(raster)
library(rgdal)
library(rgeos)
library(ggplot2)
library(rasterVis)
library(rcartocolor) # color blind friendly color palettes
##################################################################################################################
# MapSPAM data citation:
# International Food Policy Research Institute. Global Spatially-Disaggregated Crop Production Statistics Data for 2010 Version 2.0. 
# Harvard Dataverse https://doi.org/10.7910/DVN/PRFF8V (2019). 

# MapSPAM file naming conventions:
# I = irrigated
# H = rainfed, high inputs (fertilizer)
# L = rainfed, low inputs
# S = subsistence
# R =  H + L + S (all rainfed)
# A = all (irrigated + rainfed)

# file paths and variable names

# for presence/absence, use GAEZ harvested area, as this is what can be downloaded from the FAO website for year 2010
# gaez.path.2010 = "path to GAEZ v4 year 2010 harvested area data" download data from here: https://gaez-data-portal-hqfao.hub.arcgis.com/pages/theme-details-theme-5
# mapspam.path = "path to MapSPAM data /harv_area/"

# crops to compare: maize, rice, soybean, wheat
crops.short = c("MAIZ", "RICE", "SOYB", "WHEA")
crops.long  = c("Maize", "Rice", "Soybean", "Wheat")

# make a stack of difference maps
for(i in 1:4){
  # load mapspam (all: irrigated + rainfed) and gaez crop data
  mapspam = raster(paste(mapspam.path, "spam2010V2r0_global_H_", crops.short[i], "_A.tif", sep=""))
  g.2010  = raster(paste(gaez.path.2010, "GAEZAct2010_HarvArea_", crops.long[i], "_Total.tif", sep=""))
  
  # make values 1 for presence and 0 for absence of the crop
  mapspam[mapspam > 0] = 1
  g.2010[g.2010 > 0] = 1
  c.stack = stack(mapspam, g.2010)  
  
  # calculate difference
  diff = mapspam - g.2010
  
  # agreement
  sum.c.stack = stackApply(c.stack, indices = c(1,1), fun = sum)
  agree = (sum.c.stack == 2)
  
  # mapspam only
  maps.only = (diff == 1)
  
  # gaez 2010 only
  gaez.only = (diff == -1)
  
  # a single layer that identifies agreememt and differences
  diff.flat = agree*3 + maps.only*2 + gaez.only
  
  r = diff.flat
  levels(r)=data.frame(ID=1:3, code=c("GAEZ 2010 only", "MapSPAM only", "Agree"))
  
  if(i == 1){
    diff.flat.stack = r
  }else{
    diff.flat.stack = stack(diff.flat.stack, r)
  }
}
names(diff.flat.stack) = paste(crops.long)

##################################################################################################################
### Map

# Administrative unit boundaries can be downloaded from: 
#     Urbano,  Ferdinando (2018):  Global  administrative  boundaries.  
#     European  Commission,  Joint  Research  Centre  (JRC)[Dataset] 
#     PID: http://data.europa.eu/89h/jrc-10112-10004

# country = readOGR("data/gaul0_asap/", "gaul0_asap")  # country outlines

### FIGURE 10
png("results/figures/presence_validation/GAEZ_2010_vs_MAPSPAM_presence_maps.png", res=300,
    units="in", width=8, height=5.5)
rasterVis::levelplot(diff.flat.stack, 
                     par.settings = list(axis.line = list(col = "transparent")),
                     scales=list(draw=FALSE),
                     # par.strip.text = p.strip,
                     col.regions = carto_pal(n = 3, name = "Safe"),
                     colorkey=list(space="bottom", height=0.6, width=0.8), 
                     ylim=c(-60,90), ylab='', xlab='') +
  layer(sp.polygons(country, col='darkgrey', fill='transparent', lwd=0.5))
dev.off()

