#### FULL CIRCUM-ATLANTIC SUBSET ####
library(magrittr)
library(tidyverse)
library(geoChronR)
library(lubridate)
library(rgdal)


#### Extracting Atlantic Series ####

# Import the country and continent borders from the maptools package
data("wrld_simpl", package = "maptools")

# Load and extract Iso2k
# Replace the working directory to match the folder where you've locally
#   stored the Iso2k database
if (Sys.info()['sysname'] == "Darwin"){
  setwd("~/Box Sync/Konecky Lab/User Storage/Andrew Flaim/Rcode/Iso2k") # MacOS
} else if (Sys.info()['login'] == "andre"){
  setwd("C:/Users/andre/Box Sync/Konecky Lab/User Storage/Andrew Flaim/Rcode/Iso2k") # Laptop
} else{
  setwd("D:/Box/Box_Sync/Konecky Lab/User Storage/Andrew Flaim/Rcode/Iso2k") # Desktop
}

load("iso2k1_0_0.RData")
# Remove extraneous objects
rm(D, TS)
# Apply geoChronR function to extract primary TS
all_iso_ts = sTS[which(pullTsVariable(sTS, 
                                      variable = "paleoData_iso2kPrimaryTimeseries") == "TRUE")]

# Extract Atlantic series
all_iso_ts_green = all_iso_ts[which(pullTsVariable(all_iso_ts, 
                                                   variable = "geo_longitude") > -105)]
all_iso_ts_green = all_iso_ts_green[which(pullTsVariable(all_iso_ts_green, 
                                                         variable = "geo_longitude") < 50)]
all_iso_ts_green = all_iso_ts_green[which(pullTsVariable(all_iso_ts_green, 
                                                         variable = "geo_latitude") > 0)]
all_iso_ts_green = all_iso_ts_green[which(pullTsVariable(all_iso_ts_green, 
                                                         variable = "geo_latitude") < 85)]

# Define intervals of interest
mca.start = 950
mca.end = 1250

lia.start = 1650
lia.end = 1850

ant.start = 536
ant.end = 660

# Restrict to records that have data extending to at least 536 CE.
#   This threshold could be adjusted based on the desired calibration window.
endYear = 536 # Late antique little ice age start

allRecs = matrix(NA, length(all_iso_ts_green), 11) %>%
  set_colnames(c("record", "resolution", "duration", 
                 "archive", "lat", "lon", "infMat", "var", "interp", 
                 "season", "start_year"))
for(i in 1:(length(all_iso_ts_green))) {
  thisYearVec = na.omit(all_iso_ts_green[[i]]$year)
  thisYearVec = subset(thisYearVec, thisYearVec >= mca.start)
  if(max(thisYearVec) < endYear) {
    next
  } else {
    allRecs[i, 1] = all_iso_ts_green[[i]]$paleoData_TSid
    allRecs[i, 2] = (max(thisYearVec)-min(thisYearVec))/length(thisYearVec)
    allRecs[i, 3] = max(thisYearVec)-min(thisYearVec)
    allRecs[i, 4] = all_iso_ts_green[[i]]$archiveType
    allRecs[i, 5] = all_iso_ts_green[[i]]$geo_latitude
    allRecs[i, 6] = all_iso_ts_green[[i]]$geo_longitude
    allRecs[i, 7] = all_iso_ts_green[[i]]$paleoData_inferredMaterial
    allRecs[i, 8] = all_iso_ts_green[[i]]$paleoData_variableName
    allRecs[i, 9] = all_iso_ts_green[[i]]$isotopeInterpretation1_variableGroup
    allRecs[i, 10] = if(!is.null(all_iso_ts_green[[i]]$isotopeInterpretation1_seasonality)){
      all_iso_ts_green[[i]]$isotopeInterpretation1_seasonality} else {NA}
    allRecs[i, 11] = min(thisYearVec)
  }
}

allRecs = as.data.frame(allRecs) %>%
  mutate_at(c("resolution", "lat", "lon", "duration"), as.numeric)

#
#### Map All Recs ####
# Define archive point types
shapes = c("Coral" = 15, "GlacierIce" = 20, "Sclerosponge" = 15, "Wood" = 18,
           "LakeSediment" = 17, "Speleothem" = 15, "MolluskShells" = 15,
           "MarineSediment" = 17)

# Map the records by archive type and resolution
ggplot() +
  geom_polygon(data = wrld_simpl, aes(x = long, y = lat, group = group), 
               fill = "grey", colour = "black", alpha = 0.2) +
  geom_point(data = allRecs, mapping = aes(x = lon, y = lat, 
                                           shape = archive), color = "black", size = 7) +
  geom_point(data = allRecs, mapping = aes(x = lon, y = lat, 
                                               shape = archive, color = resolution), size = 6) +
  scale_shape_manual(values = shapes) + 
  scale_color_gradient2(low="blue", mid = "white", high="red") +
  # Removes Axes and labels
  scale_x_continuous(breaks = NULL) +
  xlab("") + 
  ylab("") +
  # Change theme to remove axes and ticks
  theme(panel.background = element_blank(),
        axis.ticks=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank())
#
###############################################################################
#### Average value for the MCA in each record ####
mca_include = vector()
mean_mca_value = vector()
for (i in 1:length(all_iso_ts_green)){
  if(is.null(all_iso_ts_green[[i]]$year)){
    next
  }
  thisrec = na.omit(data.frame(year = all_iso_ts_green[[i]]$year, 
                               val = all_iso_ts_green[[i]]$paleoData_values,
                               zscore = (all_iso_ts_green[[i]]$paleoData_values - mean(all_iso_ts_green[[i]]$paleoData_values))/sd(all_iso_ts_green[[i]]$paleoData_values)))
  mca_years = which(thisrec$year >= mca.start & thisrec$year <= mca.end)
  if (length(mca_years) > 0){
    mca_include = append(mca_include,i)
    mean_mca_value = append(mean_mca_value, mean(thisrec$zscore[mca_years]))
  }
}

mca_subset = all_iso_ts_green[mca_include]

mcaRecs = matrix(NA, length(mca_subset), 12) %>%
  set_colnames(c("record", "resolution", "duration", 
                 "archive", "lat", "lon", "infMat", "var", "interp", 
                 "season", "start_year", "mean_zscore"))
for(i in 1:(length(mca_subset))) {
  thisYearVec = na.omit(mca_subset[[i]]$year)
  thisYearVec = subset(thisYearVec, thisYearVec >= mca.start)
  if(max(thisYearVec) < endYear) {
    next
  } else {
    mcaRecs[i, 1] = mca_subset[[i]]$paleoData_TSid
    mcaRecs[i, 2] = (max(thisYearVec)-min(thisYearVec))/length(thisYearVec)
    mcaRecs[i, 3] = max(thisYearVec)-min(thisYearVec)
    mcaRecs[i, 4] = mca_subset[[i]]$archiveType
    mcaRecs[i, 5] = mca_subset[[i]]$geo_latitude
    mcaRecs[i, 6] = mca_subset[[i]]$geo_longitude
    mcaRecs[i, 7] = mca_subset[[i]]$paleoData_inferredMaterial
    mcaRecs[i, 8] = mca_subset[[i]]$paleoData_variableName
    mcaRecs[i, 9] = mca_subset[[i]]$isotopeInterpretation1_variableGroup
    mcaRecs[i, 10] = if(!is.null(mca_subset[[i]]$isotopeInterpretation1_seasonality)){
      mca_subset[[i]]$isotopeInterpretation1_seasonality} else {NA}
    mcaRecs[i, 11] = min(thisYearVec)
    mcaRecs[i,12] = mean_mca_value[i]
  }
}

mcaRecs = as.data.frame(mcaRecs) %>%
  mutate_at(c("resolution", "lat", "lon", "duration", "mean_zscore"), as.numeric)

#### Plot MCA mean values ####
shapes = c("Coral" = 15, "GlacierIce" = 20, "Sclerosponge" = 15, "Wood" = 18,
           "LakeSediment" = 17, "Speleothem" = 15, "MolluskShells" = 15,
           "MarineSediment" = 17)

ggplot() +
  geom_polygon(data = wrld_simpl, aes(x = long, y = lat, group = group), 
               fill = "grey", colour = "black", alpha = 0.2) +
  geom_point(data = mcaRecs, mapping = aes(x = lon, y = lat, 
                                           shape = archive), color = "black", size = 7) +
  geom_point(data = mcaRecs, mapping = aes(x = lon, y = lat, 
                                           shape = archive, color = mean_zscore), size = 6) +
  scale_shape_manual(values = shapes) + 
  #scale_fill_gradient(limits = c(-3,3)) +
  scale_color_gradient2(limits = c(-3, 3), low="blue", mid = "white", high="red") +
  # Removes Axes and labels
  scale_x_continuous(breaks = NULL) +
  xlab("") + 
  ylab("") +
  # Change theme to remove axes and ticks
  theme(panel.background = element_blank(),
        axis.ticks=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank())
#

###############################################################################
#### Average value for the LIA in each record ####
lia_include = vector()
mean_lia_value = vector()
for (i in 1:length(all_iso_ts_green)){
  if(is.null(all_iso_ts_green[[i]]$year)){
    next
  }
  thisrec = na.omit(data.frame(year = all_iso_ts_green[[i]]$year, 
                               val = all_iso_ts_green[[i]]$paleoData_values,
                               zscore = (all_iso_ts_green[[i]]$paleoData_values - mean(all_iso_ts_green[[i]]$paleoData_values))/sd(all_iso_ts_green[[i]]$paleoData_values)))
  lia_years = which(thisrec$year >= lia.start & thisrec$year <= lia.end)
  if (length(lia_years) > 0){
    lia_include = append(lia_include,i)
    mean_lia_value = append(mean_lia_value, mean(thisrec$zscore[lia_years]))
  }
}

lia_subset = all_iso_ts_green[lia_include]

liaRecs = matrix(NA, length(lia_subset), 12) %>%
  set_colnames(c("record", "resolution", "duration", 
                 "archive", "lat", "lon", "infMat", "var", "interp", 
                 "season", "start_year", "mean_zscore"))
for(i in 1:(length(lia_subset))) {
  thisYearVec = na.omit(lia_subset[[i]]$year)
  thisYearVec = subset(thisYearVec, thisYearVec >= lia.start)
  if(max(thisYearVec) < endYear) {
    next
  } else {
    liaRecs[i, 1] = lia_subset[[i]]$paleoData_TSid
    liaRecs[i, 2] = (max(thisYearVec)-min(thisYearVec))/length(thisYearVec)
    liaRecs[i, 3] = max(thisYearVec)-min(thisYearVec)
    liaRecs[i, 4] = lia_subset[[i]]$archiveType
    liaRecs[i, 5] = lia_subset[[i]]$geo_latitude
    liaRecs[i, 6] = lia_subset[[i]]$geo_longitude
    liaRecs[i, 7] = lia_subset[[i]]$paleoData_inferredMaterial
    liaRecs[i, 8] = lia_subset[[i]]$paleoData_variableName
    liaRecs[i, 9] = lia_subset[[i]]$isotopeInterpretation1_variableGroup
    liaRecs[i, 10] = if(!is.null(lia_subset[[i]]$isotopeInterpretation1_seasonality)){
      lia_subset[[i]]$isotopeInterpretation1_seasonality} else {NA}
    liaRecs[i, 11] = min(thisYearVec)
    liaRecs[i,12] = mean_lia_value[i]
  }
}

liaRecs = as.data.frame(liaRecs) %>%
  mutate_at(c("resolution", "lat", "lon", "duration", "mean_zscore"), as.numeric)

#### Plot lia mean values ####
shapes = c("Coral" = 15, "GlacierIce" = 20, "Sclerosponge" = 15, "Wood" = 18,
           "LakeSediment" = 17, "Speleothem" = 15, "MolluskShells" = 15,
           "MarineSediment" = 17)

ggplot() +
  geom_polygon(data = wrld_simpl, aes(x = long, y = lat, group = group), 
               fill = "grey", colour = "black", alpha = 0.2) +
  geom_point(data = liaRecs, mapping = aes(x = lon, y = lat, 
                                           shape = archive), color = "black", size = 7) +
  geom_point(data = liaRecs, mapping = aes(x = lon, y = lat, 
                                           shape = archive, color = mean_zscore), size = 6) +
  scale_shape_manual(values = shapes) + 
  scale_color_gradient2(limits = c(-3,3), low="blue", mid = "white", high="red") +
  # Removes Axes and labels
  scale_x_continuous(breaks = NULL) +
  xlab("") + 
  ylab("") +
  # Change theme to remove axes and ticks
  theme(panel.background = element_blank(),
        axis.ticks=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank())
#

###############################################################################
#### Average value for the LAlalia in each record ####
lalia_include = vector()
mean_lalia_value = vector()
for (i in 1:length(all_iso_ts_green)){
  if(is.null(all_iso_ts_green[[i]]$year)){
    next
  }
  thisrec = na.omit(data.frame(year = all_iso_ts_green[[i]]$year, 
                               val = all_iso_ts_green[[i]]$paleoData_values,
                               zscore = (all_iso_ts_green[[i]]$paleoData_values - mean(all_iso_ts_green[[i]]$paleoData_values))/sd(all_iso_ts_green[[i]]$paleoData_values)))
  lalia_years = which(thisrec$year >= ant.start & thisrec$year <= ant.end)
  if (length(lalia_years) > 0){
    lalia_include = append(lalia_include,i)
    mean_lalia_value = append(mean_lalia_value, mean(thisrec$zscore[lalia_years]))
  }
}

lalia_subset = all_iso_ts_green[lalia_include]

laliaRecs = matrix(NA, length(lalia_subset), 12) %>%
  set_colnames(c("record", "resolution", "duration", 
                 "archive", "lat", "lon", "infMat", "var", "interp", 
                 "season", "start_year", "mean_zscore"))
for(i in 1:(length(lalia_subset))) {
  thisYearVec = na.omit(lalia_subset[[i]]$year)
  thisYearVec = subset(thisYearVec, thisYearVec >= lalia.start)
  if(max(thisYearVec) < endYear) {
    next
  } else {
    laliaRecs[i, 1] = lalia_subset[[i]]$paleoData_TSid
    laliaRecs[i, 2] = (max(thisYearVec)-min(thisYearVec))/length(thisYearVec)
    laliaRecs[i, 3] = max(thisYearVec)-min(thisYearVec)
    laliaRecs[i, 4] = lalia_subset[[i]]$archiveType
    laliaRecs[i, 5] = lalia_subset[[i]]$geo_latitude
    laliaRecs[i, 6] = lalia_subset[[i]]$geo_longitude
    laliaRecs[i, 7] = lalia_subset[[i]]$paleoData_inferredMaterial
    laliaRecs[i, 8] = lalia_subset[[i]]$paleoData_variableName
    laliaRecs[i, 9] = lalia_subset[[i]]$isotopeInterpretation1_variableGroup
    laliaRecs[i, 10] = if(!is.null(lalia_subset[[i]]$isotopeInterpretation1_seasonality)){
      lalia_subset[[i]]$isotopeInterpretation1_seasonality} else {NA}
    laliaRecs[i, 11] = min(thisYearVec)
    laliaRecs[i,12] = mean_lalia_value[i]
  }
}

laliaRecs = as.data.frame(laliaRecs) %>%
  mutate_at(c("resolution", "lat", "lon", "duration", "mean_zscore"), as.numeric)

#### Plot lalalia mean values ####
shapes = c("Coral" = 15, "GlacierIce" = 20, "Sclerosponge" = 15, "Wood" = 18,
           "LakeSediment" = 17, "Speleothem" = 15, "MolluskShells" = 15,
           "MarineSediment" = 17)

ggplot() +
  geom_polygon(data = wrld_simpl, aes(x = long, y = lat, group = group), 
               fill = "grey", colour = "black", alpha = 0.2) +
  geom_point(data = laliaRecs, mapping = aes(x = lon, y = lat, 
                                           shape = archive), color = "black", size = 7) +
  geom_point(data = laliaRecs, mapping = aes(x = lon, y = lat, 
                                           shape = archive, color = mean_zscore), size = 6) +
  scale_shape_manual(values = shapes) + 
  scale_color_gradient2(limits = c(-3,3), low="blue", mid = "white", high="red") +
  # Removes Axes and labels
  scale_x_continuous(breaks = NULL) +
  xlab("") + 
  ylab("") +
  # Change theme to remove axes and ticks
  theme(panel.background = element_blank(),
        axis.ticks=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank())
#


# Regional boxplots of regional mean values
