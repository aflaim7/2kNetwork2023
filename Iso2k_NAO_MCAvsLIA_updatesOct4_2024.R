#### October 4, 2024 ####
# Andrew Flaim
# aflaim@wustl.edu
# 
# Iso2k MCA vs. LIA comparisons continued from PAGES Potsdam NAO meeting.
# Adjusting previous analysis to include shorter records in the LIA subset,
# as well as changing regional subsets based on modeled NAO-d18O fingerprints.

#
#### Libraries ####
library(magrittr)
library(tidyverse)
library(geoChronR)
library(lubridate)
library(rgdal)
library(progress)

#
#### Extracting Atlantic Series ####
# Import the country and continent borders from the maptools package
world <- map_data("world")
#Load Iso2k
setwd("/RAID/work/a.a.flaim/Iso2k/")
load("iso2k1_0_0.RData")
# Remove extraneous objects
rm(D, TS)
# Use geoChronR to extract primary TS
all_iso_ts = sTS[which(pullTsVariable(sTS, 
                                      variable = "paleoData_iso2kPrimaryTimeseries") == "TRUE")]

# Extract Atlantic series
all_iso_ts_green = all_iso_ts[which(pullTsVariable(all_iso_ts, 
                                                   variable = "geo_longitude") > -105)]
all_iso_ts_green = all_iso_ts_green[which(pullTsVariable(all_iso_ts_green, 
                                                         variable = "geo_longitude") < 50)]
all_iso_ts_green = all_iso_ts_green[which(pullTsVariable(all_iso_ts_green, 
                                                         variable = "geo_latitude") > 30)]
all_iso_ts_green = all_iso_ts_green[which(pullTsVariable(all_iso_ts_green, 
                                                         variable = "geo_latitude") < 85)]

# Define intervals of interest
mca_start = 1100
mca_end = 1250

lia_start = 1650
lia_end = 1850

# lalia.start = 536
# lalia.end = 660
# 
# hist.start = 1800
# hist.end = 2000

#
#### New subset for MCA duration ####
mcaRecs = matrix(NA, length(all_iso_ts_green), 12) %>%
  set_colnames(c("record", "resolution", "duration", 
                 "archive", "lat", "lon", "infMat", "var", "interp", 
                 "season", "start_year", "end_year"))
for(i in 1:(length(all_iso_ts_green))) {
  thisYearVec = na.omit(all_iso_ts_green[[i]]$year)
  # Filter for records with data covering the mca
  testrec = na.omit(data.frame(year = all_iso_ts_green[[i]]$year, 
                               val = all_iso_ts_green[[i]]$paleoData_values))
  # Test if start extends to MCA start
  if(min(testrec$year) > mca_end){next}
  if(length(which(testrec$year < mca_end & testrec$year > mca_start)) < 3){next}
  # Fill out dataframe
  mcaRecs[i, 1] = all_iso_ts_green[[i]]$paleoData_TSid
  mcaRecs[i, 2] = (max(thisYearVec)-min(thisYearVec))/length(thisYearVec)
  mcaRecs[i, 3] = max(thisYearVec)-min(thisYearVec)
  mcaRecs[i, 4] = all_iso_ts_green[[i]]$archiveType
  mcaRecs[i, 5] = all_iso_ts_green[[i]]$geo_latitude
  mcaRecs[i, 6] = all_iso_ts_green[[i]]$geo_longitude
  mcaRecs[i, 7] = all_iso_ts_green[[i]]$paleoData_inferredMaterial
  mcaRecs[i, 8] = all_iso_ts_green[[i]]$paleoData_variableName
  mcaRecs[i, 9] = all_iso_ts_green[[i]]$isotopeInterpretation1_variableGroup
  mcaRecs[i, 10] = if(!is.null(all_iso_ts_green[[i]]$isotopeInterpretation1_seasonality)){
    all_iso_ts_green[[i]]$isotopeInterpretation1_seasonality} else {NA}
  mcaRecs[i, 11] = min(thisYearVec)
  mcaRecs[i, 12] = max(thisYearVec)
  }
# Format columns
mcaRecs = as.data.frame((mcaRecs)) %>%
  mutate_at(c("resolution", "lat", "lon", "duration", "start_year"), as.numeric)
# Extract TS
mcaTS = all_iso_ts_green[which(!is.na(mcaRecs$record))]
# Remove NA
mcaRecs = mcaRecs[-which(is.na(mcaRecs$record)),]
# Remove marine sediments and d-excess records
rm_index = vector()
marine = which(mcaRecs$archive == "MarineSediment")
#excess = which(mcaRecs$var == "deuteriumExcess")
rm_index = append(rm_index, marine)
#rm_index = append(rm_index, excess)
# Clean up
mcaRecs = mcaRecs[-rm_index,]
mcaTS = mcaTS[-rm_index]
rm(rm_index, marine, excess)

#
#### New subset for lia duration ####
liaRecs = matrix(NA, length(all_iso_ts_green), 12) %>%
  set_colnames(c("record", "resolution", "duration", 
                 "archive", "lat", "lon", "infMat", "var", "interp", 
                 "season", "start_year", "end_year"))
for(i in 1:(length(all_iso_ts_green))) {
  thisYearVec = na.omit(all_iso_ts_green[[i]]$year)
  # Filter for records with data covering the lia
  testrec = na.omit(data.frame(year = all_iso_ts_green[[i]]$year, 
                               val = all_iso_ts_green[[i]]$paleoData_values))
  # Test if start extends to lia start
  if(min(testrec$year) > lia_end){next}
  if(length(which(testrec$year < lia_end & testrec$year > lia_start)) < 3){next}
  # Fill out dataframe
  liaRecs[i, 1] = all_iso_ts_green[[i]]$paleoData_TSid
  liaRecs[i, 2] = (max(thisYearVec)-min(thisYearVec))/length(thisYearVec)
  liaRecs[i, 3] = max(thisYearVec)-min(thisYearVec)
  liaRecs[i, 4] = all_iso_ts_green[[i]]$archiveType
  liaRecs[i, 5] = all_iso_ts_green[[i]]$geo_latitude
  liaRecs[i, 6] = all_iso_ts_green[[i]]$geo_longitude
  liaRecs[i, 7] = all_iso_ts_green[[i]]$paleoData_inferredMaterial
  liaRecs[i, 8] = all_iso_ts_green[[i]]$paleoData_variableName
  liaRecs[i, 9] = all_iso_ts_green[[i]]$isotopeInterpretation1_variableGroup
  liaRecs[i, 10] = if(!is.null(all_iso_ts_green[[i]]$isotopeInterpretation1_seasonality)){
    all_iso_ts_green[[i]]$isotopeInterpretation1_seasonality} else {NA}
  liaRecs[i, 11] = min(thisYearVec)
  liaRecs[i, 12] = max(thisYearVec)
}
# Format columns
liaRecs = as.data.frame((liaRecs)) %>%
  mutate_at(c("resolution", "lat", "lon", "duration", "start_year"), as.numeric)
# Extract TS
liaTS = all_iso_ts_green[which(!is.na(liaRecs$record))]
# Remove NA
liaRecs = liaRecs[-which(is.na(liaRecs$record)),]
# Remove marine sediments and d-excess records
rm_index = vector()
marine = which(liaRecs$archive == "MarineSediment")
#excess = which(liaRecs$var == "deuteriumExcess")
coral = which(liaRecs$archive == "Coral")
rm_index = append(rm_index, marine)
#rm_index = append(rm_index, excess)
rm_index = append(rm_index, coral)
# Clean up
liaRecs = liaRecs[-rm_index,]
liaTS = liaTS[-rm_index]
rm(rm_index, marine, excess)

#
#### Plot record coverage for MCA subset ####
colors = vector()
colors[which(mcaRecs$archive == "GlacierIce")] = "lightblue"
colors[which(mcaRecs$archive == "Wood")] = "brown"
colors[which(mcaRecs$archive == "Speleothem")] = "purple"
colors[which(mcaRecs$archive == "Coral")] = "coral"
colors[which(mcaRecs$archive == "Sclerosponge")] = "coral"
colors[which(mcaRecs$archive == "MolluskShells")] = "coral"
colors[which(mcaRecs$archive == "LakeSediment")] = "green"
colors[which(mcaRecs$archive == "MarineSediment")] = "blue"

# initialize
length = vector()
start_year = 1000
end_year = 2000
coverage_matrix = matrix(ncol = length(mcaTS), nrow = 2000)
coverage_cont = coverage_matrix

# order by minimum year
for(i in 1:length(mcaTS)){
  length = append(length, (min(na.omit(mcaTS[[i]]$year))))
}
ordered_mcaTS = mcaTS[(order(length))]
colors = colors[(order(length))]

# Initialize the plot and draw time period boundaries
numVar = length(mcaTS)
plot(0, 0, xlim = c(900, 2000), ylim = c(0, numVar),
     axes = F, xlab = "", ylab = "",  cex = 0.5)
par(new = T)
# Plot the years for each kept dataset
for (i in 1:numVar){
  y = rep(i, length(ordered_mcaTS[[i]]$paleoData_values))
  if (exists("year", where = ordered_mcaTS[[i]])){
    x = ordered_mcaTS[[i]]$year
  } else {x = ordered_mcaTS[[i]]$age}
  for (j in 1:length(ordered_mcaTS[[i]]$paleoData_values)){
    if (is.na(length(ordered_mcaTS[[i]]$paleoData_values[j]))){
      y[j] = -50
    }
  }
  plot(x, y, xlim = c(900, 2000), ylim = c(0, numVar),
       axes = F, xlab = "", ylab = "",  cex = 0.5, col = colors[i])#color_assign[i])#extract$colour[i])
  par(new = T)
  #print(i)
}
rect(950, -20, 1250, 275, border = "black", lwd = 3, lty = 2)
axis(1, at = seq(from = 900, to = 2000, by = 100))
axis(2)
#legend("topright", legend = rev(archives), fill = rev(color_codes))
#abline(v = 1823, lty = 2)

#
#### Plot record coverage for LIA subset ####
colors = vector()
colors[which(liaRecs$archive == "GlacierIce")] = "lightblue"
colors[which(liaRecs$archive == "Wood")] = "brown"
colors[which(liaRecs$archive == "Speleothem")] = "purple"
colors[which(liaRecs$archive == "Coral")] = "coral"
colors[which(liaRecs$archive == "Sclerosponge")] = "coral"
colors[which(liaRecs$archive == "MolluskShells")] = "coral"
colors[which(liaRecs$archive == "LakeSediment")] = "green"
colors[which(liaRecs$archive == "MarineSediment")] = "blue"

# initialize
length = vector()
start_year = 1000
end_year = 2000
coverage_matrix = matrix(ncol = length(liaTS), nrow = 2000)
coverage_cont = coverage_matrix

# order by minimum year
for(i in 1:length(liaTS)){
  length = append(length, (min(na.omit(liaTS[[i]]$year))))
}
ordered_liaTS = liaTS[(order(length))]
colors = colors[(order(length))]

# Initialize the plot and draw time period boundaries
numVar = length(liaTS)
plot(0, 0, xlim = c(900, 2000), ylim = c(0, numVar),
     axes = F, xlab = "", ylab = "",  cex = 0.5)
par(new = T)
# Plot the years for each kept dataset
for (i in 1:numVar){
  y = rep(i, length(ordered_liaTS[[i]]$paleoData_values))
  if (exists("year", where = ordered_liaTS[[i]])){
    x = ordered_liaTS[[i]]$year
  } else {x = ordered_liaTS[[i]]$age}
  for (j in 1:length(ordered_liaTS[[i]]$paleoData_values)){
    if (is.na(length(ordered_liaTS[[i]]$paleoData_values[j]))){
      y[j] = -50
    }
  }
  plot(x, y, xlim = c(900, 2000), ylim = c(0, numVar),
       axes = F, xlab = "", ylab = "",  cex = 0.5, col = colors[i])#color_assign[i])#extract$colour[i])
  par(new = T)
  print(i)
}
rect(1650, -20, 1850, 275, border = "black", lwd = 3, lty = 3)
axis(1, at = seq(from = 900, to = 2000, by = 100))
axis(2)
#legend("topright", legend = rev(archives), fill = rev(color_codes))
#abline(v = 1823, lty = 2)

#
#### Standardize and calculate MCA/LIA values ####
# MCA
mca_means = vector()
for(i in 1:length(mcaTS)){
  testrec = na.omit(data.frame(year = mcaTS[[i]]$year, 
                               val = mcaTS[[i]]$paleoData_values))
  # Trim to after 900CE
  include = which(testrec$year > 900)
  testrec = testrec[include,]
  # z-score
  testrec$val = scale(testrec$val)
  # mca subset
  testrec = testrec[which(testrec$year < mca_end & testrec$year > mca_start),]
  print(length(testrec$year))
  mca_means = append(mca_means, mean(testrec$val))
}
mcaRecs = cbind(mcaRecs, mca_means)

# LIA
lia_means = vector()
for(i in 1:length(liaTS)){
  testrec = na.omit(data.frame(year = liaTS[[i]]$year, 
                               val = liaTS[[i]]$paleoData_values))
  # Trim to after 900CE
  include = which(testrec$year > 900)
  testrec = testrec[include,]
  # z-score
  testrec$val = scale(testrec$val)
  # lia subset
  testrec = testrec[which(testrec$year < lia_end & testrec$year > lia_start),]
  lia_means = append(lia_means, mean(testrec$val))
}
liaRecs = cbind(liaRecs, lia_means)

#
#### Bin and test correlations with NAO indices ####
# Read in NAO
hernandez_NAO = read.csv("/RAID/work/a.a.flaim/Data/hernandez_nao.csv")

binvec =  seq(from = 900, to = 1999, by = 10)
binvec = as.numeric(hernandez_NAO$year)
binYears = rowMeans(cbind(binvec[-1], binvec[-length(binvec)]))

# Bin mca recs
binned_dec_mca = matrix(NA, length(binYears), length(mcaTS))
pb = progress_bar$new(total = ncol(binned_dec_mca))
for(i in 1:(ncol(binned_dec_mca))) {
  thisrec = na.omit(data.frame(year = mcaTS[[i]]$year, val = mcaTS[[i]]$paleoData_values))
  binned_dec_mca[,i] = geoChronR::bin(thisrec$year, thisrec$val, binvec)$y
  pb$tick()
  Sys.sleep(1/ncol(binned_ann_recs))
}

# Bin lia recs
binned_dec_lia = matrix(NA, length(binYears), length(liaTS))
pb = progress_bar$new(total = ncol(binned_dec_lia))
for(i in 1:(ncol(binned_dec_lia))) {
  thisrec = na.omit(data.frame(year = liaTS[[i]]$year, val = liaTS[[i]]$paleoData_values))
  binned_dec_lia[,i] = geoChronR::bin(thisrec$year, thisrec$val, binvec)$y
  pb$tick()
  Sys.sleep(1/ncol(binned_ann_recs))
}

# Bin NAO index
hernandez_bin_djf = geoChronR::bin(hernandez_NAO$year, hernandez_NAO$median, binvec)

# Test mca correlations
mca_corr = vector()
for(i in 1:dim(binned_dec_mca)[2]){
  testcor = cbind(hernandez_bin_djf$y, binned_dec_mca[,i])
  testcor = testcor[76:92,]
  testcor = na.omit(testcor)
  if (nrow(testcor) < 5){
    mca_corr = append(mca_corr, NA)
    next
  }
  corr = cor.test(testcor[,1], testcor[,2])
  mca_corr = append(mca_corr, corr$estimate)
}
mcaRecs = cbind(mcaRecs, mca_corr)

# Test lia correlations
lia_corr = vector()
for(i in 1:dim(binned_dec_lia)[2]){
  testcor = cbind(hernandez_bin_djf$y, binned_dec_lia[,i])
  testcor = testcor[76:92,]
  testcor = na.omit(testcor)
  if (nrow(testcor) < 5){
    lia_corr = append(lia_corr, NA)
    next
  }
  corr = cor.test(testcor[,1], testcor[,2])
  lia_corr = append(lia_corr, corr$estimate)
}
liaRecs = cbind(liaRecs, lia_corr)

#
#### Map MCA mean values ####
shapes = c("GlacierIce" = 20, "LakeSediment" = 17, "Speleothem" = 15, "MolluskShells" = 15,
           "MarineSediment" = 17, "Wood" = 18, "TerrestrialSediement" = 17)

ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "grey", colour = "black", alpha = 0.2) +
  geom_point(data = mcaRecs, mapping = aes(x = lon, y = lat,
                                                 shape = archive), color = "black", size = 8) +
  geom_point(data = mcaRecs, mapping = aes(x = lon, y = lat,
                                                 shape = archive, color = mca_means), size = 6) +
  scale_shape_manual(values = shapes) +
  scale_color_gradient2(limits = c(-1.5, 1.5), low="blue", mid = "white", high="red") +
  # Removes Axes and labels
  scale_x_continuous(breaks = NULL) +
  coord_equal(xlim = c(-105, 50), ylim = c(30, 90)) +
  xlab("") +
  ylab("") +
  # Change theme to remove axes and ticks
  theme(panel.background = element_blank(),
        axis.ticks=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank()) +
  # N.EU
  geom_rect(aes(xmin = -15, xmax = 40, ymin = 50, ymax = 85), 
            color = "black", fill = NA, lwd = 1.25) +
  # S.EU
  geom_rect(aes(xmin = -15, xmax = 40, ymin = 30, ymax = 50), 
            color = "black", fill = NA, lwd = 1.25) +
  # Greenland
  geom_rect(aes(xmin = -90, xmax = -15, ymin = 50, ymax = 85), 
            color = "black", fill = NA, lwd = 1.25)

#
#### Map lia mean values ####
shapes = c("GlacierIce" = 20, "LakeSediment" = 17, "Speleothem" = 15, "MolluskShells" = 15,
           "MarineSediment" = 17, "Wood" = 18, "TerrestrialSediement" = 17)

ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "grey", colour = "black", alpha = 0.2) +
  geom_point(data = liaRecs, mapping = aes(x = lon, y = lat,
                                           shape = archive), color = "black", size = 8) +
  geom_point(data = liaRecs, mapping = aes(x = lon, y = lat,
                                           shape = archive, color = lia_means), size = 6) +
  scale_shape_manual(values = shapes) +
  scale_color_gradient2(limits = c(-1.5, 1.5), low="blue", mid = "white", high="red") +
  # Removes Axes and labels
  scale_x_continuous(breaks = NULL) +
  coord_equal(xlim = c(-105, 50), ylim = c(30, 90)) +
  xlab("") +
  ylab("") +
  # Change theme to remove axes and ticks
  theme(panel.background = element_blank(),
        axis.ticks=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank()) +
  # N.EU
  geom_rect(aes(xmin = -15, xmax = 40, ymin = 50, ymax = 85), 
            color = "black", fill = NA, lwd = 1.25) +
  # S.EU
  geom_rect(aes(xmin = -15, xmax = 40, ymin = 30, ymax = 50), 
            color = "black", fill = NA, lwd = 1.25) +
  # Greenland
  geom_rect(aes(xmin = -90, xmax = -15, ymin = 50, ymax = 85), 
            color = "black", fill = NA, lwd = 1.25)
#
#### Map MCA NAO corr ####
shapes = c("GlacierIce" = 20, "LakeSediment" = 17, "Speleothem" = 15, "MolluskShells" = 15,
           "MarineSediment" = 17, "Wood" = 18, "TerrestrialSediement" = 17)

ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "grey", colour = "black", alpha = 0.2) +
  geom_point(data = mcaRecs, mapping = aes(x = lon, y = lat,
                                           shape = archive), color = "black", size = 8) +
  geom_point(data = mcaRecs, mapping = aes(x = lon, y = lat,
                                           shape = archive, color = mca_corr), size = 6) +
  scale_shape_manual(values = shapes) +
  scale_color_gradient2(limits = c(-1.5, 1.5), low="blue", mid = "white", high="red") +
  # Removes Axes and labels
  scale_x_continuous(breaks = NULL) +
  coord_equal(xlim = c(-105, 50), ylim = c(30, 90)) +
  xlab("") +
  ylab("") +
  # Change theme to remove axes and ticks
  theme(panel.background = element_blank(),
        axis.ticks=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank()) +
  # N.EU
  geom_rect(aes(xmin = -15, xmax = 40, ymin = 50, ymax = 85), 
            color = "black", fill = NA, lwd = 1.25) +
  # S.EU
  geom_rect(aes(xmin = -15, xmax = 40, ymin = 30, ymax = 50), 
            color = "black", fill = NA, lwd = 1.25) +
  # Greenland
  geom_rect(aes(xmin = -90, xmax = -15, ymin = 50, ymax = 85), 
            color = "black", fill = NA, lwd = 1.25)

#
#### Map lia NAO corr ####
shapes = c("GlacierIce" = 20, "LakeSediment" = 17, "Speleothem" = 15, "MolluskShells" = 15,
           "MarineSediment" = 17, "Wood" = 18, "TerrestrialSediement" = 17)

ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "grey", colour = "black", alpha = 0.2) +
  geom_point(data = liaRecs, mapping = aes(x = lon, y = lat,
                                           shape = archive), color = "black", size = 8) +
  geom_point(data = liaRecs, mapping = aes(x = lon, y = lat,
                                           shape = archive, color = lia_corr), size = 6) +
  scale_shape_manual(values = shapes) +
  scale_color_gradient2(limits = c(-1.5, 1.5), low="blue", mid = "white", high="red") +
  # Removes Axes and labels
  scale_x_continuous(breaks = NULL) +
  coord_equal(xlim = c(-105, 50), ylim = c(30, 90)) +
  xlab("") +
  ylab("") +
  # Change theme to remove axes and ticks
  theme(panel.background = element_blank(),
        axis.ticks=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank()) +
  # N.EU
  geom_rect(aes(xmin = -15, xmax = 40, ymin = 50, ymax = 85), 
            color = "black", fill = NA, lwd = 1.25) +
  # S.EU
  geom_rect(aes(xmin = -15, xmax = 40, ymin = 30, ymax = 50), 
            color = "black", fill = NA, lwd = 1.25) +
  # Greenland
  geom_rect(aes(xmin = -90, xmax = -15, ymin = 50, ymax = 85), 
            color = "black", fill = NA, lwd = 1.25)
#

#### Regional boxplots of regional mean values ####
# Define individual regions:
# Greenland
lon_min = -90
lon_max = -15
lat_min = 50
lat_max = 85
# Filter Greenland records
mcaRecs_green = mcaRecs[which((mcaRecs$lat > lat_min) & 
                                (mcaRecs$lat < lat_max) & 
                                (mcaRecs$lon > lon_min) & 
                                (mcaRecs$lon < lon_max)),]
liaRecs_green = liaRecs[which((liaRecs$lat > lat_min) & 
                                (liaRecs$lat < lat_max) & 
                                (liaRecs$lon > lon_min) & 
                                (liaRecs$lon < lon_max)),]

# N.EU
lon_min = -15
lon_max = 40
lat_min = 50
lat_max = 85
# Filter N.EU records
mcaRecs_NEU = mcaRecs[which((mcaRecs$lat > lat_min) & 
                                (mcaRecs$lat < lat_max) & 
                                (mcaRecs$lon > lon_min) & 
                                (mcaRecs$lon < lon_max)),]
liaRecs_NEU = liaRecs[which((liaRecs$lat > lat_min) & 
                                (liaRecs$lat < lat_max) & 
                                (liaRecs$lon > lon_min) & 
                                (liaRecs$lon < lon_max)),]

# S. EU
lon_min = -15
lon_max = 40
lat_min = 30
lat_max = 50

# Filter N.EU records
mcaRecs_SEU = mcaRecs[which((mcaRecs$lat > lat_min) & 
                                (mcaRecs$lat < lat_max) & 
                                (mcaRecs$lon > lon_min) & 
                                (mcaRecs$lon < lon_max)),]
liaRecs_SEU = liaRecs[which((liaRecs$lat > lat_min) & 
                                (liaRecs$lat < lat_max) & 
                                (liaRecs$lon > lon_min) & 
                                (liaRecs$lon < lon_max)),]

#### Create boxplots for MCA and LIA and Hist ###
box.list = list("Grn MCA" = mcaRecs_green$mca_means, "Grn LIA" = liaRecs_green$lia_means,
                "N.EU MCA" = mcaRecs_NEU$mca_means, "S.EU MCA" = mcaRecs_SEU$mca_means,
                "N.EU LIA" = liaRecs_NEU$lia_means, "S.EU LIA" = liaRecs_SEU$lia_means)
#"N.EU Hist" = histRecs_NEU$mean_zscore, "S.EU Hist" = histRecs_SEU$mean_zscore)

boxplot(box.list, outline = F, ylim = c(-1.5, 1.5))

stripchart(box.list,
           method = "jitter", 
           col = c("black"), 
           pch = 16,
           vertical = T,
           add = T)

#
