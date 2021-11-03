# Compare FAOSTAT 2015 harvested area to GAEZ 2015 harvested area
# Danielle S Grogan
# 2021

##################################################################################################################
### LIBRARIES AND SOURCE FILES ###
##################################################################################################################
library(raster)
library(rgdal)
library(rgeos)

##################################################################################################################
# function to aggregate values by country using GAUL country brick
country_agg = function(country.brick,
                       country.names,
                       brk,
                       out.nm){
  
  for(c in 1:nlayers(country.brick)){
    
    # make a single country layer
    country.raster = subset(country.brick, c) 
    country.raster[country.raster == 0] = NA
    grid_mult = brk * country.raster  # multiply crop layers by a single country layer
    
    out.data = as.data.frame(t(cellStats(grid_mult, sum)))
    out.data = cbind(as.character(country.names$ADM0_NAME[c]), out.data)
    
    if(c == 1){
      colnames(out.data)[1] = "ADM0_NAME"
      colnames(out.data)[2:ncol(out.data)] = names(brk) # crops
      write.table(out.data, out.nm, quote=FALSE, sep=",")
    }else{
      write.table(out.data, out.nm, quote=FALSE, sep=",", append = TRUE, row.names = F, col.names = F)
    }
    print(c)
  }
}


##################################################################################################################
### Load Country Data ###
##################################################################################################################
# country brick: use the FAO Global Administrative Unit (GAUL) 2012 data
# country.brick = brick("path to GAUL 2012 gridded data")

# country data - link country.brick level names (which are ADM0 codes) to "human readable" country names
cell.fraction = read.delim("data/ExtractedCellTableRAWv3.csv", sep="\t") # this file is provided in the data/ folder of the git repo
country.names = subset(cell.fraction, select=c("ADM0_CODE", "ADM0_NAME", "DISP_AREA"))
d = duplicated(country.names)
country.names = country.names[d==F,]

# slight modification to names and extent of country brick
names(country.brick) = country.names$ADM0_CODE

##################################################################################################################
### Spatial aggregation ###
##################################################################################################################

### GAEZ 2015 harvested darea: download data from https://doi.org/10.7910/DVN/KAGRFI 
gaez.2015.path = "path to GAEZ+2015 harvested area data folder"
gaez.files = list.files(gaez.2015.path, pattern='.tif', full.names = T)
g.stack = do.call(stack,
                  lapply(gaez.files, 
                         raster::brick))

# rename layers with crop names
for(i in 1:length(gaez.files)){
  names(g.stack)[i] = sub("GAEZAct2015_HarvArea_", "", sub(".tif", "", strsplit(gaez.files[i], "/")[[1]][11]))
}

country_agg(country.brick,
            country.names,
            brk = g.stack,
            out.nm = "results/harvested_area/Harvested_area_by_country_GAEZ_2015.csv")


##################################################################################################################
### Read & summarize country/crop aggregates ###
##################################################################################################################

# read gaez aggregates by crop and country
gaez.2015.harv = read.csv("results/harvested_area/Harvested_area_by_country_GAEZ_2015.csv")

# summarize by crop (for Table 3)
gaez.2015.harv.crop = as.data.frame(colSums(gaez.2015.harv[,2:ncol(gaez.2015.harv)]))
gaez.2015.harv.crop.total = subset(gaez.2015.harv.crop, grepl("Total", rownames(gaez.2015.harv.crop)))
rownames(gaez.2015.harv.crop.total) = sub("_Total", "", rownames(gaez.2015.harv.crop.total))
colnames(gaez.2015.harv.crop.total) = c("Harvested Area 1000 ha")
write.csv(gaez.2015.harv.crop.total, "results/harvested_area/Harvested_area_by_crop_global_total_GAEZ_2015.csv")

# summarize by country (for Table 4)
gaez.2015.harv.total = subset(gaez.2015.harv, select = c(grepl("Total", colnames(gaez.2015.harv))))
rownames(gaez.2015.harv.total) = gaez.2015.harv$ADM0_NAME

# remove foddercrops because FAO does not report them
rm.col = which(colnames(gaez.2015.harv.total) == "Foddercrops_Total")
gaez.2015.harv.total = subset(gaez.2015.harv.total, select = -c(rm.col))  

gaez.2015.harv.country = as.data.frame(rowSums(gaez.2015.harv.total))
colnames(gaez.2015.harv.country) = c("Harvested Area 1000 ha")
write.csv(gaez.2015.harv.country, "results/harvested_area/Harvested_area_by_country_total_GAEZ_2015.csv")

##################################################################################################################
### FAO data ###
##################################################################################################################

# load crop list that matches FAO crop names to GAEZ crop names
# this data is provided in the data/ folder of the git repo, and was downloaded from FAOSTAT: fao.org/faostat/en
fao.harv.data = read.csv("data/FAOSTAT/CompareGAEZFAOSTAT_kha_2.csv") 
fao.harv.2015 = fao.harv.data[,c(1,2,3,10)]

# summarize by crop (for Table 3)
crops = unique(fao.harv.2015$GAEZNameX)
fao.harv.crop = data.frame(matrix(nr = length(crops), nc=1))
rownames(fao.harv.crop) = crops
for(c in 1:length(crops)){
  fao.harv.sub = subset(fao.harv.2015, fao.harv.2015$GAEZNameX == crops[c])
  fao.harv.crop[c,1] = sum(fao.harv.sub$Avg2015, na.rm=T)/1000 # divide by 1000 to convert from ha to 1000 ha
}
rownames(fao.harv.crop) = sub("_Total", "", rownames(fao.harv.crop))
colnames(fao.harv.crop) = c("Harvested Area 1000 ha")
write.csv(fao.harv.crop, "results/harvested_area/Harvested_area_by_crop_global_total_FAO_2015.csv")

# summarize by country (for Table 4)
countries = unique(fao.harv.2015$ADM0_NAME)
fao.harv.country = data.frame(matrix(nr = length(countries), nc=1))
rownames(fao.harv.country) = countries
for(c in 1:length(countries)){
  fao.harv.sub = subset(fao.harv.2015, fao.harv.2015$ADM0_NAME == countries[c])
  fao.harv.country[c,1] = sum(fao.harv.sub$Avg2015, na.rm=T)/1000 # divide by 1000 to convert from ha to 1000 ha
  print(c)
}
colnames(fao.harv.country) = c("Harvested Area 1000 ha")
write.csv(fao.harv.country, "results/harvested_area/Harvested_area_by_country_total_FAO_2015.csv")


##################################################################################################################
### # Compare harvested area data ###
##################################################################################################################

# By crop (Table 3)

# load gaez and fao 2015 harvested area by crop
gaez.harv = read.csv("results/harvested_area/Harvested_area_by_crop_global_total_GAEZ_2015.csv")
fao.harv  = read.csv("results/harvested_area/Harvested_area_by_crop_global_total_FAO_2015.csv")

# Crop data harmonization: in FAO, Fruits_&_Nuts are reported separately, but in GAEZ they are combined with CropsNES
fao.fruitNut = fao.harv[which(fao.harv$X == "Fruits_&_Nuts"),2]           # find harvested area for Fruits_&_Nuts
fao.CropsNES = fao.harv[which(fao.harv$X == "CropsNES"),2]                # find harvested area for CropsNES
fao.harv[which(fao.harv$X == "CropsNES"),2] = fao.CropsNES + fao.fruitNut # sum the two together, replace CropsNES with the summed value
fao.harv = fao.harv[-c(which(fao.harv$X == "Fruits_&_Nuts")),]  # remove the Fruits_&_Nuts category

# Fix crop list names
gaez.harv$X = sub("Yamsandotherroots", "Yams_and_other_roots", gaez.harv$X)
gaez.harv$X = sub("PotatoAndSweetpotato", "Potato_&_Sweet_potato", gaez.harv$X)
gaez.harv$X = sub("Othercereals", "Other_cereals", gaez.harv$X)
gaez.harv$X = sub("Oilpalmfruit", "Oil_palm_fruit", gaez.harv$X)

# merge by crop names
harv.merge = merge(gaez.harv, fao.harv, by = "X")

# calculate differences
colnames(harv.merge) = c("Crop", "GAEZ_2015_1000ha", "FAO_2015_1000ha")
harv.merge$GAEZ_diff_1000ha = harv.merge$GAEZ_2015_1000ha - harv.merge$FAO_2015_1000ha
harv.merge$GAEZ_diff_FAO_percent = 100*(harv.merge$GAEZ_diff_1000ha/harv.merge$FAO_2015_1000ha)

# sort by largest area
harv.merge = harv.merge[order(-harv.merge$GAEZ_2015_1000ha),]

# save to file
write.csv(harv.merge, "results/harvested_area/Table_3_GAEZ_vs_FAO_2015_crop_global_total.csv", row.names = F)



#################################################################################################
# By country (Table 4)

# load gaez and fao 2015 harvested area by crop
gaez.harv = read.csv("results/harvested_area/Harvested_area_by_country_total_GAEZ_2015.csv")
fao.harv  = read.csv("results/harvested_area/Harvested_area_by_country_total_FAO_2015.csv")

# fix country name mismatches
fao.harv$X = as.character(fao.harv$X)
gaez.harv$X = as.character(gaez.harv$X)
fao.harv$X[64]   = "Côte d'Ivoire"
fao.harv$X[201]  = "Réunion"
fao.harv$X[160]  = "Moldova Republic of"
gaez.harv$X[193] = "Moldova Republic of"

# merge by country names
harv.merge = merge(gaez.harv, fao.harv, by = "X", all=T)

# calculate differences
colnames(harv.merge) = c("Country", "GAEZ_2015_1000ha", "FAO_2015_1000ha")
harv.merge$GAEZ_diff_1000ha = harv.merge$GAEZ_2015_1000ha - harv.merge$FAO_2015_1000ha
harv.merge$GAEZ_diff_FAO_percent = 100*(harv.merge$GAEZ_diff_1000ha/harv.merge$FAO_2015_1000ha)

# sort by largest area
harv.merge = harv.merge[order(-harv.merge$GAEZ_2015_1000ha),]

# save to file
write.csv(harv.merge, "results/harvested_area/Table_4_GAEZ_vs_FAO_2015_country_total.csv", row.names = F)
