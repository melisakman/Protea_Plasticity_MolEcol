
# Calculate mean plasticity values for all maternal family lines
# 11/19/17
# Andrew Latimer

setwd("~/Google Drive/ProteaPaper/Andrew_analyses_makeover")

library(dplyr)
library(lme4)
library(lmerTest)
library(R.utils)
library(fields)

# Load data

# load trait data with Rachel's MDS values included
d <- read.csv("repens_drydown_alldata_mds.csv")
# get rid of junk rows
d <- d[,-grep("X.", names(d))]
# get rid of rows with no population value (not measured)
d <- d[!is.na(d$pop),]

# Separate data for sampling time 1 (day 6) and sampling time 2 (day 12)
d1 <- filter(d, round==1)
d2 <- filter(d, round==2)

# check combinations of maternal line and treatment
xtabs(~MaternalLine+drought_water, d1)
xtabs(~MaternalLine+drought_water, d2)

# Traits to calculate plasticity for 
d1traits <- c("height_cm", "conductance", "stomatal_size", "stomatal_density", "SLA", "LWR", "area_cm2", "H4", "dataME0", "dataME1","dataME2", "dataME3", "dataME4", "dataME5", "dataME6", "dataME7", "dataME8", "dataME9", "dataME10", "dataME11", "dataME12", "dataME13")

d2traits <- c("height_cm", "conductance", "stomatal_size", "stomatal_density", "SLA", "LWR", "area_cm2", "H4", "growth_no_neg", "root_length", "soluble_shoot_all_g", "soluble_root_all_g", "starch_shoot_all_g", "starch_root_all_g", "dataME0", "dataME1","dataME2", "dataME3", "dataME4", "dataME5", "dataME6", "dataME7", "dataME8", "dataME9", "dataME10", "dataME11", "dataME12", "dataME13")

# Function to calculate plasticity
# (family mean for wet treatment - family mean for dry treatment) / (family mean across both treatments)
calc.plast <- function(trt1, trt2) {
  return((trt1-trt2) / mean(c(trt1, trt2, na.rm=T)))
}

# Calculate family means for wet treatment and dry treatments
d1_mean <- aggregate(d1, by=list(d1$MaternalLine, d1$drought_water), FUN=mean, na.rm=T)
d1_mean$MaternalLine <- d1_mean$Group.1
d1_mean$drought_water <- d1_mean$Group.2
d2_mean <- aggregate(d2, by=list(d2$MaternalLine, d2$drought_water), FUN=mean, na.rm=T)
d2_mean$MaternalLine <- d2_mean$Group.1
d2_mean$drought_water <- d2_mean$Group.2


# Probably there's a better way, but for now, calculate family-level plasticities by looping over family

# Function to generate plasticity values from family mean data sets
get_plast_trait_vals <- function(famdata, traits) {
          # famdata is a data frame with family means for 2 treatments. Has to contain a column "MaternalLine", a column "drought_water" with treatment names "drought" and "water", and columns with family mean trait values
  # traits contains the names of the columns with the trait values
  fams <- unique(famdata$MaternalLine)
  n.fams <- length(fams)
  plast <- data.frame(MaternalLine=fams, pop=sapply(fams, substr, start=1, stop=1)) # storage data frame
  for (i in 1:n.fams) {
    trt1.row <- which(famdata$MaternalLine == fams[i] & famdata$drought_water == "water")
    trt2.row <- which(famdata$MaternalLine == fams[i] & famdata$drought_water == "drought")
    for (j in 1:length(traits)) {
      if (length(trt1.row) & length(trt2.row)) plast[i,traits[j]] <- calc.plast(famdata[trt1.row,traits[j]], famdata[trt2.row,traits[j]]) else plast[i,traits[j]] <- NA
    }
  }
return(plast)
}

# Do the calculations for both days 

d1_plast <- get_plast_trait_vals(famdata=d1_mean, traits=d1traits)
head(d1_plast)

d2_plast <- get_plast_trait_vals(famdata=d2_mean, traits=d2traits)
head(d2_plast)

# Checking how many measurements per pop
d2_plast$conductance_1 <- !is.na(d2_plast$conductance)
xtabs(~pop+conductance_1, d2_plast)
# No values for Alicedale on Day 2 for conductance
# Otherwise, between 2 and 5 family-level plasticity values per population. 

# Write results

write.csv(d1_plast, "family_plasticity_d1.csv")
write.csv(d2_plast, "family_plasticity_d2.csv")
