# compare FAOSTAT harvested area to GAEZ 2010 harvested area
# Danielle S Grogan
# 2021-01

##################################################################################################################
### Load Data ###
# NB: this data is not included in the GAEZ+_2015 data package
# FAOSTAT harvest area data can be downloaded from: http://www.fao.org/faostat/en/#data
# the GAEZ 2010 harvest area data citation is:  
#     Fischer G, Nachtergaele FO, Prieler S, Teixeira E, Toth G, van Velthuizen H, Verelst L, & Wiberg D (IIASA/FAO). 2012. 
#     Global Agro-ecological Zones (GAEZ v3.0). IIASA, Laxenburg, Austria and FAO, Rome, Italy.  
#     Available at http://pure.iiasa.ac.at/id/eprint/13290/.
# GAEZ v4 2010 data is available from the authors (Fischer et al 2012) upon request
##################################################################################################################

# load harvest area data
# harv.data = read.csv("compiled FAOSTAT and GAEZ 2010 data (not provided))

# subset to columns needed
harv.area = harv.data[,c("ADM0_NAME", "GAEZ_Name", "GAEZ_v4_2010_area_ha", "Avg2010")]

# rename columns for easier reading
colnames(harv.area)[2] = "GAEZ_crop"
colnames(harv.area)[4] = "FAO_2010av_area"

# subset to rows with no NAs
harv.area = na.omit(harv.area)

##################################################################################################################
### Table 3: Compare by Crop Category ###
##################################################################################################################

crops = unique(harv.area$GAEZ_crop)
out.table = data.frame(matrix(nr=length(crops), nc=9))
out.table[,1] = crops
colnames(out.table) = c("ADM0_NAME", "GAEZ_2010_total_ha", "FAO_2010_total_ha", 
                        "GAEZ_minus_FAO_ha", "GAEZ_diff_FAO_percent",
                        "R2", "slope", "signif", "sigma")

for(c in 1:length(crops)){
  
  harv.crop = subset(harv.area, harv.area$GAEZ_crop == crops[c])

  s = summary(lm(harv.crop$GAEZ_v4_2010_area_ha ~ harv.crop$FAO_2010av_area))
  out.table$GAEZ_2010_total_ha[c] = sum(harv.crop$GAEZ_v4_2010_area_ha)
  out.table$FAO_2010_total_ha[c]  = sum(harv.crop$FAO_2010av_area)

  out.table$R2[c]              = s$r.squared
  out.table$slope[c]           = s$coefficients[2,1]
  out.table$signif[c]          = s$coefficients[2,4]
  out.table$sigma[c]           = s$sigma
}

out.table$GAEZ_minus_FAO_ha = out.table$GAEZ_2010_total_ha - out.table$FAO_2010_total_ha
out.table$GAEZ_diff_FAO_percent = 100*(out.table$GAEZ_minus_FAO_ha/out.table$FAO_2010_total_ha)

write.csv(out.table, "results/harvest_area_validation/GAEZ_2010_vs_FAOav_2010_by_crop.csv", row.names=F)

##################################################################################################################
### Table 4: Compare by Country ###
##################################################################################################################

countries = unique(harv.area$ADM0_NAME)
out.table = data.frame(matrix(nr=length(countries), nc=9))
out.table[,1] = countries
colnames(out.table) = c("ADM0_NAME", "GAEZ_2010_total_ha", "FAO_2010_total_ha", 
                        "GAEZ_minus_FAO_ha", "GAEZ_diff_FAO_percent",
                        "R2", "slope", "signif", "sigma")

for(c in 1:length(countries)){
  
  harv.country = subset(harv.area, harv.area$ADM0_NAME == countries[c])
  
  if(nrow(harv.country) > 1 &                       # require > 1 data point
     sum(harv.country$GAEZ_v4_2010_area_ha > 0)){   # require non-zero data
    
    s = summary(lm(harv.country$GAEZ_v4_2010_area_ha ~ harv.country$FAO_2010av_area))
    out.table$GAEZ_2010_total_ha[c] = sum(harv.country$GAEZ_v4_2010_area_ha)
    out.table$FAO_2010_total_ha[c]  = sum(harv.country$FAO_2010av_area)
    
    out.table$R2[c]              = s$r.squared
    out.table$slope[c]           = s$coefficients[2,1]
    out.table$signif[c]          = s$coefficients[2,4]
    out.table$sigma[c]           = s$sigma
  }
}

out.table$GAEZ_minus_FAO_ha = out.table$GAEZ_2010_total_ha - out.table$FAO_2010_total_ha
out.table$GAEZ_diff_FAO_percent = 100*(out.table$GAEZ_minus_FAO_ha/out.table$FAO_2010_total_ha)

write.csv(out.table, "results/harvest_area_validation/GAEZ_2010_vs_FAOav_2010_by_country.csv", row.names=F)

# EOF
