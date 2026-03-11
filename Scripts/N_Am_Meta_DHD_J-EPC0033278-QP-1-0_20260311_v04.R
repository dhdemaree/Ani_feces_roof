###### Script for full meta analysis Data Tables for Paper that includes North America #######
# By David Demaree
# ORISE Participant, US EPA, LASMB, Athens, GA.
# 2/1/2024

# The purpose of this R script is to conduct meta-analysis of prevalence data
# and generate a table of all the distributions.  Specifically this scrpt
# is for all the figures that I will be using for the journal article that I am
# going to be writing.

# This top script is the current full form the change log of past changes is
# below the script

# This script is large so make sure to run it locally on your computer and
# do not run it if files are on a network drive or OneDrive.  For EPA windows 
# computers make sure to right click on the work directory file that you are 
# using for this script and make sure that "Always keep on this device" is
# checked. Otherwise some parts of the script will take a long time to run 15+ min.

# Also while I could have made a script that will do everything I am not because
# of runtimes.

# This script is differet from the other scripts because it seperates out data
# from north America

# For this script I made some wording changes so I am going to re do a bunch of
# figures I am going to save up changes until I get to the end of this draft.

#-----------------------------
# Changes

# Change the word pathogen in all the figures to microbe.
# Change the prevalence scale from 0-1 to percentage.

#########################################################

# Directory

# This script sets the directory to the folder this script is in.
# This script will only work in rstudio and won't work in base R
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()

# Directory that was used for the project
# getwd()
# setwd("C:/Users/DDEMAREE/OneDrive - Environmental Protection Agency (EPA)/FI to Pathogen RainRWH/R Stuff/Metaanalysis (Complete)")






# Packages
library(fitdistrplus)
library(forestplot)
library(dplyr) 
library(tidyverse)
library(meta)
library(mc2d)
library(ggplot2)

# Sets up parallel processing to speed up aggregation
library(foreach)
library(doParallel)

# Remove data
rm(list = ls() )

# Parallel Processing set up.
cl <- detectCores()-5                                                             # Detects the number of cores and uses half of them I am doing this because I want to leave a little more for other stuff.
registerDoParallel(cl)                                                            # Registers the cores so that they can be used in do Parallel
getDoParWorkers()                                                                 # Detects number of cores in cluster




### Run once
# This is for stuff that should only be run one time to collect all the data.
FinalEnumList <-list()                                                             # Creates a list that will be used to create the table that has all the log normal means and standard deviations
FinalPrevalencelist <- list()                                                      # Creates the final prevalence list that everything will be saved to.   

#####Functions to calculate log normal mean and standard deviation#####

EnumFunAMeanASD <- function(Amean, ASD) {                                         # Arithmetic mean and Arithmetic Standard Deviation
  lnlnSD <- sqrt(log(1+((ASD^2)/(Amean^2))))
  lnlnmean <- log(Amean)-((lnlnSD^2)/2)
  lnmedian <- exp(lnlnmean) 
  return(c(lnmedian, lnlnmean, lnlnSD))
}

EnumFunGmeanASD <- function(geomedian, ASD) {                                     # Geometric mean/median and Standard Deviation
  lnlnSD <- sqrt(log(1+((ASD^2)/(geomedian^2)) - ((ASD^4)/(geomedian^4))))
  lnlnmean <- log(geomedian)
  lnmedian <- exp(lnlnmean) 
  return(c(lnmedian, lnlnmean, lnlnSD))
}

EnumFunModeASD <- function(Amode, ASD) {                                          # Mode and SD 
  lnlnSD <- sqrt(log(1+((ASD^2)/(Amode^2))- ((3*(ASD^4))/(Amode^4))))
  lnlnmean <- log(Amode)+(lnlnSD^2)
  lnmedian <- exp(lnlnmean) 
  return(c(lnmedian, lnlnmean, lnlnSD))
}

EnumFunMinMaxN <- function(Amin, Amax, Total) {                                   # Minimum, Maximum, and sample size 
  lnlnmean <- (log(Amin)+log(Amax))/2
  lnlnSD <- (log(Amin)+log(Amax))/(qnorm(Total/(Total+1))-qnorm(1/(Total+1)))
  lnmedian <- exp(lnlnmean) 
  return(c(lnmedian, lnlnmean, lnlnSD))
}

EnumFunMinMax <- function(Amin, Amax) {                                           # Minimum and Maximum
  lnlnmean <- (log(Amin)+log(Amax))/2
  lnlnSD <- (log(Amax)-log(Amin))/sqrt(12)
  lnmedian <- exp(lnlnmean) 
  return(c(lnmedian, lnlnmean, lnlnSD))
}

EnumFun95CL <- function(AConf95H, AConf95L) {                                     # 95% Confidence intervals
  lnlnmean <- (log(AConf95H)+log(AConf95L))/2                                     # This assumes that the 95% confidence intervals are similar to the min and max
  lnlnSD <- (log(AConf95H)-log(AConf95L))/3.92
  lnmedian <- exp(lnlnmean) 
  return(c(lnmedian, lnlnmean, lnlnSD))
}

EnumFunGmeanLODBelow <- function(geomedian, LOD, LODBelow) {                # Median/geometic mean, Limit of detection, percent below the limit of detection
  lnlnmean <- log(geomedian)
  lnlnSD <- (log(LOD)-log(geomedian))/qnorm(LODBelow)
  lnmedian <- exp(lnlnmean) 
  return(c(lnmedian, lnlnmean, lnlnSD))
}

EnumFunMeanLODBelow <- function(Amean, LOD, LODBelow) {                      # Mean, Limit of detection, percent below the limit of detection
  lnlnSD <- qnorm(LODBelow)+sqrt(((qnorm(LODBelow))^2)+(2*(log(Amean)-log(LOD))))
  lnlnmean <- log(Amean)-((lnlnSD^2)/2)
  lnmedian <- exp(lnlnmean) 
  return(c(lnmedian, lnlnmean, lnlnSD))
}

### GRADE function
#GRADE vector
GRADEvector <- c("Very Low","Low","Moderate","High")                              # Keeping this seperate from the function so that it doesn't keep recreating the vector
GRADEfunct <- function(x) {                                                       # Converts the score to a grade scores are from 0-3 that represent from "Very Low" to "High"
  a <- x+1                                                                        # Adding +1 because the table goes from 1-4 and GRADE scores from 0-3 
  GRADE <- GRADEvector[a]                                            
  return(c(GRADE))
}

###################
# At start of repeat

# Import different dataset each time you repeat

# Import different dataset each time you repeat
# I did it this way instead of using a forloop because the dataset was so big
# so I wanted to have some stops so that I could check the data and fix errors
# before continuing.

# Import data
# Coliform
RawData <- read.csv("Rawdata/Coliform_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv")

# E. Coli 
RawData <- read.csv("Rawdata/EColi_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv")

# Enterococci
RawData <- read.csv("Rawdata/Enteriococci_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv")

# Salmonella 
RawData <- read.csv("Rawdata/Salmonella_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv")

# Cryptosporidium 
RawData <- read.csv("Rawdata/Cryptosporidium_STUDYN_DHD_J-EPC0033278-QP-1-0_220230627_v03.csv")

# Giardia 
RawData <- read.csv("Rawdata/Giardia_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv")

# Campylobacter 
RawData <- read.csv("Rawdata/Campylobacter_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv")


#### Everything from here on is repeated verbatimum

# Change the raw data from characters to numeric
RawData["Presence.Absence"][RawData["Presence.Absence"] == "ND"] <- "0" 
RawDatacols.num <- c("MDL.ALOD.LOD", 
                     "Presence.Absence",
                     "Total",
                     "Individual.Sample",
                     "Mean", 
                     "Geometric.mean", 
                     "Standard.deviation", 
                     "Median", 
                     "Min", 
                     "Max", 
                     "X95..low",
                     "X95..high",
                     "CenMedianLow",
                     "CenMedianHigh",
                     "CenMinLow",
                     "CenMinHigh",
                     "CenMaxLow",
                     "CenMaxHigh",
                     "GRADES_Score")
RawData [RawDatacols.num] <- sapply(RawData[RawDatacols.num],as.numeric)

# This removed values that do not have an animal classification at the taxonomic
# levels that we are analyzing for.  For this script we are removing the rows that
# do not have an Animal common name which is the species level.  Uses tidy

RawData <- subset(RawData, grepl('[a-zA-Z]', class))

# This assigns number to everything so that a forloop can run through all of it.

RawData$Studies <- paste(RawData$Reference, RawData$Year, sep = " " )             # Study Citation
RawData$Studies.num <- as.numeric(as.factor(RawData$Studies))                     # Study Number
RawData$Animal.num <- as.numeric(as.factor(RawData$class))                        # Animal Number
RawData$Path.num <- as.numeric(as.factor(RawData$Pathogen))                       # Pathogen Number
RawData$Method.num <- as.numeric(as.factor(RawData$Conventional..Molecular..or.Immunologic))  # Method Type Number
RawData <- mutate(RawData,                                                        # Enum Type Number
                  Enum.num = case_when(
                    Enum.Stats == "StatsTotal"       ~ "1",
                    Enum.Stats == "StatsnoND"        ~ "2",
                    Enum.Stats == "Stats"            ~ "3",
                    Enum.Stats == "StatsIndiv"       ~ "4",
                    Enum.Stats == "Semiquantitative" ~ "5",
                    TRUE                             ~ "NA",
                  )
)
RawData$Enum.num <- as.numeric(RawData$Enum.num)
RawData <- mutate(RawData,                                                        # Censor Number
                  Cen.Num = case_when(
                    Censored.Data == "yes"       ~ "1",
                    Censored.Data == "no"        ~ "2",
                    TRUE                         ~ "NA",
                  )
)
RawData$Cen.Num <- as.numeric(RawData$Cen.Num)


# Numbers to set the bounds of the forloops.  I am doing it this way because I
# encountered issues when I entered the function directly into the forloop
# R can only display 50 columns at one time you have to press the arrows to see 
# the other columns

StudyNum <- max(RawData$Studies.num)
AniNum  <- max(RawData$Animal.num)
MethNum <- max(RawData$Method.num)
PathNum <- max(RawData$Path.num)
EnumNum <- max(RawData$Enum.num, na.rm = TRUE)                                    # Added na.rm to remove NA values.  NA was removed because this column only has values for enumerated data not presence data
CenNum <- max(RawData$Cen.Num, na.rm = TRUE)

Prevalencelist <- list()

# Prevalence based on method-pathogen-animal
for (i in 1:MethNum) {
  for (j in 1:PathNum) {
    for (k in 1:AniNum) {
      AniSubset <- RawData[(RawData$Animal.num==k & 
                              RawData$Presence=='yes' & 
                              RawData$Path.num==j & 
                              RawData$Method.num ==i),]
      if (nrow(AniSubset) > 0){
        AniSubset$Studies <- paste(AniSubset$Reference, AniSubset$Year, sep = " " )
        AniSubset$Studies.Num <- as.numeric(as.factor(AniSubset$Studies))
        AniMeta <- metaprop(event = Presence.Absence,
                            n = Total,
                            studlab = Studies,
                            data = AniSubset,
                            method = "Inverse",
                            sm = "PLOGIT",
                            common = TRUE,
                            random = TRUE,
                            hakn = TRUE,
                            title = "Pathogen Meta analysis",
                            control=list(stepadj=0.5, maxiter=1000))
        AniMetaEff <- c(AniMeta$TE.common,
                        AniMeta$lower.common,
                        AniMeta$upper.common,
                        AniMeta$TE.random,
                        AniMeta$lower.random,
                        AniMeta$upper.random)
        AniMetaEFF1 <- meta:::backtransf(AniMetaEff, sm = "PLOGIT")
        GRADES <- min(AniSubset$GRADES_Score)                                     #Takes the lowest GRADE
        GRADES <- GRADEfunct(GRADES)                                              #Gives a GRADE rating to the score
        AniStudiesUn <-unique(AniSubset$Studies)
        AniStudiesSub <- paste(unlist(AniStudiesUn), collapse = ", ")      # Collects the study names
        AniMetaDat <- c(AniSubset$class[[1]], 
                        AniSubset$Pathogen[[1]], 
                        AniSubset$Conventional..Molecular..or.Immunologic [[1]],
                        GRADES, 
                        AniMeta$k, 
                        max(AniSubset$Studies.Num),
                        sum(AniSubset$Presence.Absence), 
                        sum(AniSubset$Total), 
                        AniMetaEFF1,
                        AniSubset$Animal.num[[1]],
                        AniSubset$Method.num[[1]],
                        AniSubset$Path.num[[1]],
                        AniStudiesSub)
        Prevalencelist[[length(Prevalencelist)+1]] <- AniMetaDat
        print(paste(i,MethNum,"/",j,PathNum,"/",k,AniNum,"/"))
      }else{next}
    }
  }
}


# Presence for conventional methods, total animals, and by pathogen #####

for (i in 1:MethNum) {
  for (j in 1:PathNum) {
    AniSubset <- RawData[(RawData$Presence=='yes' & 
                            RawData$Path.num==j & 
                            RawData$Method.num ==i),]
    if (nrow(AniSubset) > 0){
      AniSubset$Studies <- paste(AniSubset$Reference, AniSubset$Year, sep = " " )
      AniSubset$Studies.Num <- as.numeric(as.factor(AniSubset$Studies))
      AniMeta <- metaprop(event = Presence.Absence,
                          n = Total,
                          studlab = Studies,
                          data = AniSubset,
                          method = "Inverse",
                          sm = "PLOGIT",
                          common = TRUE,
                          random = TRUE,
                          hakn = TRUE,
                          title = "Pathogen Meta analysis",
                          control=list(stepadj=0.5, maxiter=1000))
      AniMetaEff <- c(AniMeta$TE.common,
                      AniMeta$lower.common,
                      AniMeta$upper.common,
                      AniMeta$TE.random,
                      AniMeta$lower.random,
                      AniMeta$upper.random)
      AniMetaEFF1 <- meta:::backtransf(AniMetaEff, sm = "PLOGIT")
      GRADES <- min(AniSubset$GRADES_Score)                                     #Takes the lowest GRADE
      GRADES <- GRADEfunct(GRADES)                                              #Gives a GRADE rating to the score
      AniStudiesUn <-unique(AniSubset$Studies)
      AniStudiesSub <- paste(unlist(AniStudiesUn), collapse = ", ")      # Collects the study names
      AniMetaDat <- c("Total", 
                      AniSubset$Pathogen[[1]], 
                      AniSubset$Conventional..Molecular..or.Immunologic [[1]],
                      GRADES, 
                      AniMeta$k, 
                      max(AniSubset$Studies.Num),
                      sum(AniSubset$Presence.Absence), 
                      sum(AniSubset$Total), 
                      AniMetaEFF1,
                      AniSubset$Animal.num[[1]],
                      AniSubset$Method.num[[1]],
                      AniSubset$Path.num[[1]],
                      AniStudiesSub)
      Prevalencelist[[length(Prevalencelist)+1]] <- AniMetaDat
      print(paste(i,MethNum,"/",j,PathNum,"/",k,AniNum,"/"))
    }else{next}
  }
}

##### Presence in North America #####

for (i in 1:MethNum) {
  for (j in 1:PathNum) {
    AniSubset <- RawData[(RawData$Presence=='yes' & 
                            RawData$Path.num==j & 
                            RawData$Method.num ==i &
                            RawData$Global.Region == 'North America'),]
    if (nrow(AniSubset) > 0){
      AniSubset$Studies <- paste(AniSubset$Reference, AniSubset$Year, sep = " " )
      AniSubset$Studies.Num <- as.numeric(as.factor(AniSubset$Studies))
      AniMeta <- metaprop(event = Presence.Absence,
                          n = Total,
                          studlab = Studies,
                          data = AniSubset,
                          method = "Inverse",
                          sm = "PLOGIT",
                          common = TRUE,
                          random = TRUE,
                          hakn = TRUE,
                          title = "Pathogen Meta analysis",
                          control=list(stepadj=0.5, maxiter=1000))
      AniMetaEff <- c(AniMeta$TE.common,
                      AniMeta$lower.common,
                      AniMeta$upper.common,
                      AniMeta$TE.random,
                      AniMeta$lower.random,
                      AniMeta$upper.random)
      AniMetaEFF1 <- meta:::backtransf(AniMetaEff, sm = "PLOGIT")
      GRADES <- min(AniSubset$GRADES_Score)                                     #Takes the lowest GRADE
      GRADES <- GRADEfunct(GRADES)                                              #Gives a GRADE rating to the score
      AniStudiesUn <-unique(AniSubset$Studies)
      AniStudiesSub <- paste(unlist(AniStudiesUn), collapse = ", ")      # Collects the study names
      AniMetaDat <- c("North America", 
                      AniSubset$Pathogen[[1]], 
                      AniSubset$Conventional..Molecular..or.Immunologic [[1]],
                      GRADES, 
                      AniMeta$k, 
                      max(AniSubset$Studies.Num),
                      sum(AniSubset$Presence.Absence), 
                      sum(AniSubset$Total), 
                      AniMetaEFF1,
                      AniSubset$Animal.num[[1]],
                      AniSubset$Method.num[[1]],
                      AniSubset$Path.num[[1]],
                      AniStudiesSub)
      Prevalencelist[[length(Prevalencelist)+1]] <- AniMetaDat
      print(paste(i,MethNum,"/",j,PathNum,"/",k,AniNum,"/"))
    }else{next}
  }
}





FinalPrevalencelist <- c(FinalPrevalencelist, Prevalencelist)
PrevalenceData <- as.data.frame(do.call(rbind, Prevalencelist))
colnames(PrevalenceData) <- c("Animal", 
                              "Pathogen", 
                              "Method",
                              "GRADE", 
                              "Datasets",
                              "Studies",
                              "Positive",
                              "Sample(n)", 
                              "Fixed.effects.mean",
                              "Fixed.effects.95low",
                              "Fixed.effects.95high",
                              "Random.effects.mean",
                              "Random.effects.95low",
                              "Random.effects.95high",
                              "AniNum",
                              "MethNum",
                              "PathNum",
                              "Reference") 

# These lines change the columns in the new Prevalence Data to numeric since
# everything was set to character after compiling the table.

PrevalenceData.num <- c("Datasets",
                        "Studies",
                        "Positive",
                        "Sample(n)", 
                        "Fixed.effects.mean",
                        "Fixed.effects.95low",
                        "Fixed.effects.95high",
                        "Random.effects.mean",
                        "Random.effects.95low",
                        "Random.effects.95high",
                        "AniNum",
                        "MethNum",
                        "PathNum")
PrevalenceData[PrevalenceData.num] <- sapply(PrevalenceData[PrevalenceData.num],as.numeric)


##### Enumerated Data #####

### The script separated into steps ###

# Changing Rawdata to the subset that only has Enumerated datat
RawdataE <- RawData[(RawData$Enumerated =='yes'),]
# Creates the list for the enumerated data
Enumerationlist <- list()                                                         

### Stats and StasTotal 
EnumSubset <- RawdataE[(RawdataE$Enum.Stats=="StatsTotal" |                       #subset based on StatsTotal and Stats 
                          RawdataE$Enum.Stats=="Stats"),]
if (nrow(EnumSubset) > 0){ 
  Subrow <- nrow(EnumSubset)
  for (o in 1:Subrow) {
    Amean <- EnumSubset$Mean[o]                                                     # "arithmetic mean"
    ASD <- EnumSubset$Standard.deviation[o]                                         # arithmetic standard deviation
    Amode <- NA                                                                     # mode no data in lit review
    Total <- EnumSubset$Total[o]                                                    # sample size
    AConf95H <- EnumSubset$X95..high[o]                                             # upper 95% confidence interval
    AConf95L <- EnumSubset$X95..low[o]                                              # lower 95% confidence interval
    LOD  <- EnumSubset$MDL.ALOD.LOD[o]                                              # Limit of detection
    LoDPrev <- EnumSubset$Presence.Absence[o]                                       # Samples above Limit of detection
    LODBelow <- 1 - (LoDPrev/Total)
    Amax <- EnumSubset$Max[o]
    Log10Mean <- EnumSubset$log10mean[o]
    Log10SD <- EnumSubset$log10SD[o]
    if(!is.na(EnumSubset$Median[o])) {                                              # median or geometric mean
      geomedian <- EnumSubset$Median[o] 
    } else if (!is.na(EnumSubset$Geometric.mean[o])){
      geomedian <- EnumSubset$Geometric.mean[o]
    } else if (EnumSubset$CenMedianLow[o] == 0 &                                   # Median by substitution methods
               !is.na(EnumSubset$CenMedianHigh[o])){
      geomedian <- EnumSubset$CenMedianHigh[o]
      geomedian <- geomedian/sqrt(2)
    } else {
      geomedian <- NA
    }
    if(!is.na(EnumSubset$Min[o])){                                                 # Minimum
      Amin <- EnumSubset$Min[o]
    } else if (!is.na(EnumSubset$CenMinHigh[o])){                                     # Minimum by substitution methods
      Amin <- EnumSubset$CenMinHigh[o]
      Amin <- Amin/sqrt(2)
    } else {
      Amin <- NA
    }
    
    if(!is.na(Log10Mean) & 
       !is.na(Log10SD)){
      Log10Data <- rnorm(n = 10000, mean = Log10Mean, sd = Log10SD)
      Log10Data <- 10^Log10Data
      Log10EnumFit <- fitdist(Log10Data, "lnorm")
      EnumSubLNmeanSD <- c(exp(Log10EnumFit$estimate[[1]]), Log10EnumFit$estimate[[1]], Log10EnumFit$estimate[[2]], Total, "No")
    } else if(!is.na(Amean) &&                                                #Entering the created variables in to the Equations
              !is.na(ASD) ) { 
      EnumSubLNmeanSD <- EnumFunAMeanASD (Amean, ASD)               # Mean and ASD
      EnumSubLNmeanSD <- c(EnumSubLNmeanSD, Total, "No")            # Adds total sample size to EnumSubLNmeanSD
    } else if (!is.na(geomedian) && 
               !is.na(ASD) && 
               (geomedian >= ASD)) {
      EnumSubLNmeanSD <- EnumFunGmeanASD(geomedian, ASD)            # Geometric mean/median and Standard Deviation
      EnumSubLNmeanSD <- c(EnumSubLNmeanSD, Total, "No")            # Adds total sample size to EnumSubLNmeanSD
    } else if (!is.na(Amode) && 
               !is.na(ASD)) {
      EnumSubLNmeanSD <- EnumFunModeASD(Amode, ASD)                 # Mode and SD
      EnumSubLNmeanSD <- c(EnumSubLNmeanSD, Total, "No")            # Adds total sample size to EnumSubLNmeanSD
    } else if (!is.na(Amin) && 
               !is.na(Amax) && !is.na(Total)) {
      EnumSubLNmeanSD <- EnumFunMinMaxN(Amin, Amax, Total)          # Minimum, Maximum, and sample size
      EnumSubLNmeanSD <- c(EnumSubLNmeanSD, Total, "No")            # Adds total sample size to EnumSubLNmeanSD
    } else if (!is.na(Amin) && 
               !is.na(Amax)) {
      EnumSubLNmeanSD <- EnumFunMinMax(Amin, Amax)                  # Minimum and Maximum
      EnumSubLNmeanSD <- c(EnumSubLNmeanSD, Total, "No")            # Adds total sample size to EnumSubLNmeanSD
    } else if (!is.na(AConf95L) && 
               !is.na(AConf95H)) {
      EnumSubLNmeanSD <- EnumFun95CL(AConf95H, AConf95L)            # 95% Confidence intervals
      EnumSubLNmeanSD <- c(EnumSubLNmeanSD, Total, "No")            # Adds total sample size to EnumSubLNmeanSD
    } else if (!is.na(geomedian) && 
               !is.na(LOD) && 
               !is.na(LODBelow)) {
      if (LODBelow < 0.5){
        EnumSubLNmeanSD <- EnumFunGmeanLODBelow(geomedian, LOD, LODBelow)              # Median/geometic mean, Limit of detection, percent below the limit of detection
        EnumSubLNmeanSD <- c(EnumSubLNmeanSD, Total, "No")                            # Adds total sample size to EnumSubLNmeanSD
      } else {
        EnumSubLNmeanSD <- c("High Presence", NA, NA, Total, NA)
      }
    } else if (!is.na(Amean) && 
               !is.na(LOD) && 
               !is.na(LODBelow)) {
      EnumSubLNmeanSD <- EnumFunMeanLODBelow(Amean, LOD, LODBelow)      # Mean, Limit of detection, percent below the limit of detection
      EnumSubLNmeanSD <- c(EnumSubLNmeanSD, Total, "No")             # Adds total sample size to EnumSubLNmeanSD
    } else if (!is.na(geomedian) &&                                                    # Imputed (geomedian, LOD, %below LOD)
               !is.na(LOD)) {
      AniImput <- EnumSubset$Animal.num[o]                                             # This is a little funky where I am creating variables from the EnumSubset data frame
      PathImput <- EnumSubset$Path.num[o]
      MethImput <- EnumSubset$Method.num[o]
      LoDSub <- PrevalenceData[(PrevalenceData$AniNum == AniImput &                     # Then using the variables to create a subset fo the PrevaleanceData data frame.
                                  PrevalenceData$PathNum == PathImput &
                                  PrevalenceData$MethNum == MethImput ),]              # Added method type since it was not in the original, I might still remove it if I am going to combine all the prevalence data.
      if (nrow(LoDSub) > 0){
        LODBelow <- 1-LoDSub$Fixed.effects.mean[1]                                     # Then I am pulling from that data frame to find imputed values for the function.
        if (LODBelow < 0.5){
          EnumSubLNmeanSD <- EnumFunGmeanLODBelow(geomedian, LOD, LODBelow)              # Median/geometic mean, Limit of detection, percent below the limit of detection
          EnumSubLNmeanSD <- c(EnumSubLNmeanSD, Total, "Yes")                            # Adds total sample size to EnumSubLNmeanSD
        } else {
          EnumSubLNmeanSD <- c("High Presence", NA, NA, Total, NA)}
      } else {
        EnumSubLNmeanSD <- c("No Presence", NA, NA, Total, NA)                         # This is to give a warning that there was no presence metadata that could be used for the analysis
      }
    } else if (!is.na(Amean) && 
               !is.na(LOD)) {
      AniImput <- EnumSubset$Animal.num[o]                                             # This is a little funky where I am creating variables from the EnumSubset data frame
      PathImput <- EnumSubset$Path.num[o]
      MethImput <- EnumSubset$Method.num[o]                                            
      LoDSub <- PrevalenceData[(PrevalenceData$AniNum == AniImput &                     # Then using the variables to create a subset fo the PrevaleanceData data frame.
                                  PrevalenceData$PathNum == PathImput &
                                  PrevalenceData$MethNum == MethImput ),]
      if (nrow(LoDSub) > 0){
        LODBelow <- 1-LoDSub$Fixed.effects.mean[1]  
        EnumSubLNmeanSD <- EnumFunMeanLODBelow(Amean, LOD, LODBelow)                   # Mean, Limit of detection, percent below the limit of detection
        EnumSubLNmeanSD <- c(EnumSubLNmeanSD, Total, "Yes")                            # Adds total sample size to EnumSubLNmeanSD
      } else {
        EnumSubLNmeanSD <- c("No Presence", NA, NA, Total, NA)                         # This is to give a warning that there was no presence metadata that could be used for the analysis
      }                                         
    } else {EnumSubLNmeanSD <- c("No Function Stats Total", NA, NA, Total, NA)}
    EnumOutData <- c(EnumSubset$Studies[o],                             #This adds data to the output data sheet
                     EnumSubset$class[o],
                     EnumSubset$Pathogen[o],
                     EnumSubset$Conventional..Molecular..or.Immunologic[o],
                     EnumSubset$Enum.Stats[o],
                     EnumSubLNmeanSD,   
                     EnumSubset$Enum.num[o],
                     EnumSubset$Studies.num[o],
                     EnumSubset$Animal.num[o],
                     EnumSubset$Method.num[o],
                     EnumSubset$Path.num[o],
                     EnumSubset$GRADES_Score[o],
                     EnumSubset$USA[o],
                     EnumSubset$Global.Region[o])
    Enumerationlist[[length(Enumerationlist)+1]] <- EnumOutData
  }
}

# StatsIndiv

RawdataEI <- RawdataE[(RawdataE$Enum.Stats=="StatsIndiv"), ] 

if (nrow(RawdataEI) > 0){
  RawdataEI$Studies.numsub <- as.numeric(as.factor(RawdataEI$Studies))                     # Study Number
  RawdataEI$Animal.numsub <- as.numeric(as.factor(RawdataEI$class))           # Animal Number
  RawdataEI$Path.numsub <- as.numeric(as.factor(RawdataEI$Pathogen))                       # Pathogen Number
  
  StudyNumEI <- max(RawdataEI$Studies.numsub)
  AniNumEI  <- max(RawdataEI$Animal.numsub)
  PathNumEI <- max(RawdataEI$Path.numsub)
  
  for (n in 1:CenNum) {
    for (i in 1:MethNum) {
      for (j in 1:PathNumEI) {
        for (k in 1:AniNumEI) {
          for (l in 1:StudyNumEI) {
            EnumSubset <- RawdataEI[(RawdataEI$Cen.Num==n &
                                       RawdataEI$Path.numsub==j & 
                                       RawdataEI$Method.num ==i & 
                                       RawdataEI$Studies.numsub==l & 
                                       RawdataEI$Animal.numsub==k),]
            if (nrow(EnumSubset) > 0){                                            # If there is no data in a subset it moves on to the next condition
              if (nrow(EnumSubset) > 1){                                           # Remove data if there is only one data point
                if (EnumSubset$Censored.Data[1] == "no" ) {                        # no Censored Data in whole data sets
                  FitEnum <- EnumSubset$Individual.Sample
                  FitEnum1 <- fitdist(FitEnum, "lnorm")
                } else if (EnumSubset$Censored.Data[1]  == "yes" ) {               # Censored Data if there is any censored data point in the dataset
                  EnumIndivDataList <- list()
                  IndivRow <- nrow(EnumSubset)
                  for (p in 1:IndivRow) {
                    if (!is.na(EnumSubset$Individual.Sample[p])){
                      EnumIndiv <- c(EnumSubset$Individual.Sample[p],
                                     EnumSubset$Individual.Sample[p])
                    } else if (is.na(EnumSubset$Individual.Sample[p])){
                      EnumIndiv <- c(EnumSubset$Min[p],
                                     EnumSubset$Max[p])
                    } else {EnumIndiv <- c("No Individual data", NA)}
                    EnumIndivDataList[[length(EnumIndivDataList)+1]] <- EnumIndiv   #Creates List
                    EnumIndivData <- as.data.frame(do.call(rbind, EnumIndivDataList)) # Creates data frame from list
                    colnames(EnumIndivData) <- c("left",
                                                 "right")
                  }
                  FitEnum1 <- try(fitdistcens(EnumIndivData, "lnorm"),TRUE)
                }
                if (class(FitEnum1) == "try-error"){
                  next
                } else {
                  EnumSubLNmeanSD <- as.numeric(c(exp(FitEnum1$estimate[1]), 
                                                  FitEnum1$estimate[1], 
                                                  FitEnum1$estimate[2]))
                }
              } else {EnumSubLNmeanSD <- c("One Sample", NA, NA)}
              EnumSubLNmeanSD <- c(EnumSubLNmeanSD, nrow(EnumSubset), "No")           #Adds total to individual data
              EnumOutData <- c(EnumSubset$Studies[1],                             #This adds data to the output data sheet
                               EnumSubset$class[1],
                               EnumSubset$Pathogen[1],
                               EnumSubset$Conventional..Molecular..or.Immunologic[1],
                               EnumSubset$Enum.Stats[1],
                               EnumSubLNmeanSD,   
                               EnumSubset$Enum.num[1],
                               EnumSubset$Studies.num[1],
                               EnumSubset$Animal.num[1],
                               EnumSubset$Method.num[1],
                               EnumSubset$Path.num[1],
                               EnumSubset$GRADES_Score[1],
                               EnumSubset$USA[1],
                               EnumSubset$Global.Region[1])
              Enumerationlist[[length(Enumerationlist)+1]] <- EnumOutData
            }else{next}
          }
        }
      }
    }
  }
}

### Semiquantitative

RawdataES <- RawdataE[RawdataE$Enumerated =='yes' &
                        RawdataE$Enum.Stats=="Semiquantitative",]

if (nrow(RawdataES) > 0){
  RawdataES$Studies.numsub <- as.numeric(as.factor(RawdataES$Studies))                     # Study Number
  RawdataES$Animal.numsub <- as.numeric(as.factor(RawdataES$class))           # Animal Number
  RawdataES$Path.numsub <- as.numeric(as.factor(RawdataES$Pathogen))                       # Pathogen Number
  
  StudyNumES <- max(RawdataES$Studies.numsub)
  AniNumES  <- max(RawdataES$Animal.numsub)
  PathNumES <- max(RawdataES$Path.numsub)
  
  for (n in 1:CenNum) {
    for (i in 1:MethNum) {
      for (j in 1:PathNumES) {
        for (k in 1:AniNumES) {
          for (l in 1:StudyNumES) {
            EnumSubset <- RawdataES[(RawdataES$Cen.Num==n &
                                       RawdataES$Path.numsub==j & 
                                       RawdataES$Method.num ==i & 
                                       RawdataES$Studies.numsub==l & 
                                       RawdataES$Animal.numsub==k),]
            if (nrow(EnumSubset) > 0){                                            # If there is no data in a subset it moves on to the next condition
              EnumSemiQDataList <-list()
              SemiQRow <- nrow(EnumSubset)
              for (p in 1:SemiQRow) {                                             # Adds semiquantitative ranges to the table.
                EnumIndiv <- c(EnumSubset$Min[p],
                               EnumSubset$Max[p])
                EnumSemiQDataList[[length(EnumSemiQDataList)+1]] <- EnumIndiv   #Creates List
                EnumSemiQData <- as.data.frame(do.call(rbind, EnumSemiQDataList)) # Creates data frame from list
                colnames(EnumSemiQData) <- c("left",
                                             "right")
              }
              duptimes <- EnumSubset$Presence.Absence                               # the number of time s a row is duplicated
              idx <- rep(1:nrow(EnumSemiQData), duptimes)                           # Creates a vector of the replicated rows
              EnumSemiQData <- EnumSemiQData[idx,]                                  # Weird filtering method that replicates the columns 
              FitEnum1 <- try(fitdistcens(EnumSemiQData, "lnorm"),TRUE)           # Generating the Stats
              if (class(FitEnum1) == "try-error"){
                next
              } else {
                EnumSubLNmeanSD <- c(as.numeric(c(exp(FitEnum1$estimate[1]), 
                                                  FitEnum1$estimate[1], 
                                                  FitEnum1$estimate[2],
                                                  EnumSubset$Total[1])), 
                                     "No")
              }
              EnumOutData <- c(EnumSubset$Studies[1],                             #This adds data to the output data sheet
                               EnumSubset$class[1],
                               EnumSubset$Pathogen[1],
                               EnumSubset$Conventional..Molecular..or.Immunologic[1],
                               EnumSubset$Enum.Stats[1],
                               EnumSubLNmeanSD,   
                               EnumSubset$Enum.num[1],
                               EnumSubset$Studies.num[1],
                               EnumSubset$Animal.num[1],
                               EnumSubset$Method.num[1],
                               EnumSubset$Path.num[1],
                               EnumSubset$GRADES_Score[1],
                               EnumSubset$USA[1],
                               EnumSubset$Global.Region[1])
              Enumerationlist[[length(Enumerationlist)+1]] <- EnumOutData
            }else{next}
          }
        }
      }
    }
  }
}

##### Summary Statistics that do not include non detects #####
# Special case: Where there is a mean or individual value and only one positive
# sample

EnumSubset <- RawdataE[RawdataE$Enum.Stats=="StatsnoND" &
                         RawdataE$Presence.Absence == 1 &                           # If there is 1 presence data point and a mean then it is assumed that that mean represents the value of the sample
                         !is.na(RawdataE$Individual.Sample),]                      # I am making the assumption that the mean is always avalible if the is only one datapoint.

if (nrow(EnumSubset) > 0){                                                          # Since this is a rare case I added this so that it will skip over this step without doing it.
  Subrow <- nrow(EnumSubset)
  for (o in 1:Subrow) {                                             # This for loop was added just incase there was a situation where there were muiltiple rows of data which is an artifact of data entry where I can have multiple data sets on the same thing from the same study.
    ImputefitcensList <-list()
    RowData <- c(EnumSubset$Individual.Sample[o],                        # If the value is over the LLOD it adds a row with the values that set a left and right range
                 EnumSubset$Individual.Sample[o]) 
    ImputefitcensList[[length(ImputefitcensList)+1]] <- RowData   #Creates List
    ImpNum <- EnumSubset$Total[o] - 1                                                     # Creates the max number of runs that the for loop is going to run. The number of values under the limit of detection 
    for (i in 1:ImpNum) {                                                           # This loop adds a bunch of rows with the censored values.
      RowData <- c(NA,                                                              # This creates a vector for censored data with NA and the log of the LLOD.  NA is used since log 10 goes to negative infinity.  
                   EnumSubset$MDL.ALOD.LOD[o])
      ImputefitcensList[[length(ImputefitcensList)+1]] <- RowData
    }
    # Distributions Perameters
    Imputefitcens <- as.data.frame(do.call(rbind, ImputefitcensList))                 # Creates data frame from list
    colnames(Imputefitcens) <- c("left",                                              # Adds header to dataframe
                                 "right")
    ImpFitCensPeram <- fitdistcens(Imputefitcens, "lnorm")                             # This fits a distribution to the censored data in the table created above.
    Impmean <- as.numeric(ImpFitCensPeram$estimate[1])                                # Mean from the distfitcens
    ImpSD <- as.numeric(ImpFitCensPeram$estimate[2])                                  # Standard deviation from distfitcens                                                                       # 
    LogLLODImp <- EnumSubset$MDL.ALOD.LOD[1]                                          # The Limit of Detection 
    ImpTotDataList <- list()
    for (i in 1:10000) {                                                              # For loop that generates 10000 dataset
      c <- rtrunc("rlnorm", n=ImpNum, mean = Impmean, sd = ImpSD, lsup = LogLLODImp) # Generating the random numbers for the censored values using rtrunc from mc2d. n is the number of values generated, "norm" for normal distribution, mean for mean, sd for standard devation, lsup truncates the data set to below a certain value.
      d <- na.omit(Imputefitcens$left)                                                   # The non censored values
      e <- c(c,d)                                                                     # Creating a complete data set
      f <- fitdist(e, "lnorm")
      ImpTotRow <- c(median(e), 
                     f$estimate[1], 
                     f$estimate[2])                                                   # Creates a vector of the mean and standard deviation of the complete data set
      ImpTotDataList[[length(ImpTotDataList)+1]] <- ImpTotRow                         # Adds Row to a list
      ImpTotData <- as.data.frame(do.call(rbind, ImpTotDataList))                     # Creates dataframe from the list
      colnames(ImpTotData) <- c("Median",
                                "Mean",
                                "SD")
    }
    EnumSubLNmeanSD <- c(mean(ImpTotData$Median),                  # Output for the meta analysis, Median
                         mean(ImpTotData$Mean),                    # Mean
                         mean(ImpTotData$SD),                      # Standard Deviation
                         EnumSubset$Total[o],                      # Total sample Size
                         "Yes")                                    # If imputation was used to create results
    EnumOutData <- c(EnumSubset$Studies[o],                             #This adds data to the output data sheet
                     EnumSubset$class[o],
                     EnumSubset$Pathogen[o],
                     EnumSubset$Conventional..Molecular..or.Immunologic[o],
                     EnumSubset$Enum.Stats[o],
                     EnumSubLNmeanSD,   
                     EnumSubset$Enum.num[o],
                     EnumSubset$Studies.num[o],
                     EnumSubset$Animal.num[o],
                     EnumSubset$Method.num[o],
                     EnumSubset$Path.num[o],
                     EnumSubset$GRADES_Score[o],
                     EnumSubset$USA[o],
                     EnumSubset$Global.Region[o])
    Enumerationlist[[length(Enumerationlist)+1]] <- EnumOutData
  }
}

# Special case where there are there are two values over LOD and am min and a
# max which have to be the two values.

EnumSubset <- RawdataE[RawdataE$Enum.Stats=="StatsnoND" &
                         RawdataE$Presence.Absence == 2 &                           # If there are two data points this assumes that the min and the max are the two data points.
                         !is.na(RawdataE$Min) &
                         !is.na(RawdataE$Max),]                      

if (nrow(EnumSubset) > 0){                                            # If there is no data in a subset it moves on to the next condition
  Subrow <- nrow(EnumSubset)                                          # length of subset
  for (o in 1:Subrow) {                                               # This for loop was added just incase there was a situation where there were muiltiple rows of data which is an artifact of data entry where I can have multiple data sets on the same thing from the same study.
    ImputefitcensList <- list()                              # Creates a new data.frame for fitdistcens.
    
    RowData <- c(EnumSubset$Min[o],                        # If the value is over the LLOD it adds a row with the values that set a left and right range
                 EnumSubset$Min[o]) 
    ImputefitcensList[[length(ImputefitcensList)+1]] <- RowData
    RowData <- c(EnumSubset$Max[o],                        # If the value is over the LLOD it adds a row with the values that set a left and right range
                 EnumSubset$Max[o]) 
    ImputefitcensList[[length(ImputefitcensList)+1]] <- RowData            
    ImpNum <- EnumSubset$Total[o] - 2                                                     # Creates the max number of runs that the for loop is going to run
    for (i in 1:ImpNum) {                                                           # This loop adds a bunch of rows with the censored values.
      RowData <- c(NA,                                                              # This creates a vector for censored data with NA and the log of the LLOD.  NA is used since log 10 goes to negative infinity.  
                   EnumSubset$MDL.ALOD.LOD[o])
      ImputefitcensList[[length(ImputefitcensList)+1]] <- RowData
    }
    Imputefitcens <- as.data.frame(do.call(rbind, ImputefitcensList))                 # Creates data frame from list
    colnames(Imputefitcens) <- c("left",
                                 "right")
    # Distributions Perameters
    ImpFitCensPeram <- fitdistcens(Imputefitcens, "lnorm")                             # This fits a distribution to the censored data in the table created above.
    Impmean <- as.numeric(ImpFitCensPeram$estimate[1])                                # Mean from the distfitcens
    ImpSD <- as.numeric(ImpFitCensPeram$estimate[2])                                  # Standard deviation from distfitcens
    LogLLODImp <- EnumSubset$MDL.ALOD.LOD[1]                                          # The Limit of Detection      
    ImpTotDataList <- list()                                                          # Creates a new data.frame for fitdistcens.
    for (i in 1:10000) {                                                              # For loop that generates 10000 dataset
      c <- rtrunc("rlnorm", n=ImpNum, mean = Impmean, sd = ImpSD, lsup = LogLLODImp) # Generating the random numbers for the censored values using rtrunc from mc2d. n is the number of values generated, "norm" for normal distribution, mean for mean, sd for standard devation, lsup truncates the data set to below a certain value.
      d <- na.omit(Imputefitcens$left)                                                   # The non censored values
      e <- c(c,d)                                                                     # Creating a complete data set
      f <- fitdist(e, "lnorm")
      ImpTotRow <- c(median(e), 
                     f$estimate[1], 
                     f$estimate[2])                                                   # Creates a vector of the mean and standard deviation of the complete data set
      ImpTotDataList[[length(ImpTotDataList)+1]] <- ImpTotRow                         # Adds Row to a list
    }
    ImpTotData <- as.data.frame(do.call(rbind, ImpTotDataList))                     # Creates dataframe from the list
    colnames(ImpTotData) <- c("Median",
                              "Mean",
                              "SD")
    EnumSubLNmeanSD <- c(mean(ImpTotData$Median),                  # Output for the meta analysis, Median
                         mean(ImpTotData$Mean),                    # Mean
                         mean(ImpTotData$SD),                      # Standard Deviation
                         EnumSubset$Total[o],                      # Total sample Size
                         "Yes")                                    # If imputation was used to create results
    EnumOutData <- c(EnumSubset$Studies[o],                                           #This adds data to the output data sheet
                     EnumSubset$class[o],
                     EnumSubset$Pathogen[o],
                     EnumSubset$Conventional..Molecular..or.Immunologic[o],
                     EnumSubset$Enum.Stats[o],
                     EnumSubLNmeanSD,   
                     EnumSubset$Enum.num[o],
                     EnumSubset$Studies.num[o],
                     EnumSubset$Animal.num[o],
                     EnumSubset$Method.num[o],
                     EnumSubset$Path.num[o],
                     EnumSubset$GRADES_Score[o],
                     EnumSubset$USA[o],
                     EnumSubset$Global.Region[o])
    Enumerationlist[[length(Enumerationlist)+1]] <- EnumOutData
  }
}

# General case where there are >2 presence values that do not fall under
# any of the special cases above

EnumSubset <- RawdataE[RawdataE$Enum.Stats=="StatsnoND" &
                         RawdataE$Presence.Absence >= 2 &                           
                         RawdataE$Censored.Data == "no" &
                         !is.na(RawdataE$MDL.ALOD.LOD),] 

if (nrow(EnumSubset) > 0){                                            # If there is no data in a subset it moves on to the next condition
  Subrow <- nrow(EnumSubset)                                          # length of subset
  for (o in 1:Subrow) {                                               # Applying the functions on each row not by subset.
    Amean <- EnumSubset$Mean[o]                         # "arithmetic mean"
    ASD <- EnumSubset$Standard.deviation[o]                           # arithmetic standard deviation
    if(!is.na(EnumSubset$Median[o])) {                                # median or geometric mean
      geomedian <- EnumSubset$Median[o] 
    } else if (!is.na(EnumSubset$Geometric.mean[o])){
      geomedian <- EnumSubset$Geometric.mean[o]
    } else {geomedian <- EnumSubset$Median[o]}
    Amode <- NA                                                       # mode no data in lit review
    Amin <- EnumSubset$Min[o]                                         # minimum
    Amax <- EnumSubset$Max[o]                                         # maximum
    Total <- EnumSubset$Total[o]                                      # sample size
    AConf95H <- EnumSubset$X95..high[o]                               # upper 95% confidence interval
    AConf95L <- EnumSubset$X95..low[o]                                # lower 95% confidence interval
    LOD  <- EnumSubset$MDL.ALOD.LOD[o]                                # Limit of detection
    LoDPrev <- EnumSubset$Presence.Absence[o]                         # Samples above Limit of detection
    LoDPrevBelow <- Total - EnumSubset$Presence.Absence[1]
    LODBelow <- 1 - (LoDPrev/Total)
    if(LoDPrev == 2 &&                                                # I am keeping this because censored data should use a lot of the same variables and go row by row.
       !is.na(Amin) && 
       !is.na(Amax)) {
    } else if(!is.na(Amean) &&                                                # I am keeping this because censored data should use a lot of the same variables and go row by row.
              !is.na(ASD) ) { 
      EnumSubLNmeanSD <- EnumFunAMeanASD (Amean, ASD)               # Mean and ASD
    } else if (!is.na(geomedian) && 
               !is.na(ASD) && 
               (geomedian >= ASD)) {
      EnumSubLNmeanSD <- EnumFunGmeanASD(geomedian, ASD)            # Geometric mean/median and Standard Deviation
    } else if (!is.na(Amode) && 
               !is.na(ASD)) {
      EnumSubLNmeanSD <- EnumFunModeASD(Amode, ASD)                 # Mode and SD
    } else if (!is.na(Amin) && 
               !is.na(Amax) && !is.na(Total)) {
      EnumSubLNmeanSD <- EnumFunMinMaxN(Amin, Amax, Total)          # Minimum, Maximum, and sample size
    } else if (!is.na(Amin) && 
               !is.na(Amax)) {
      EnumSubLNmeanSD <- EnumFunMinMax(Amin, Amax)                  # Minimum and Maximum
    } else if (!is.na(AConf95L) && 
               !is.na(AConf95H)) {
      EnumSubLNmeanSD <- EnumFun95CL(AConf95H, AConf95L)            # 95% Confidence intervals
    } else if (!is.na(geomedian) && 
               !is.na(LOD) && 
               !is.na(LODBelow)) {
      EnumSubLNmeanSD <- EnumFunGmeanLODBelow(geomedian, LOD, LODBelow)        # Median/geometic mean, Limit of detection, percent below the limit of detection
    } else if (!is.na(Amean) && 
               !is.na(LOD) && 
               !is.na(LODBelow)) {
      if (LODBelow < 0.5){
        EnumSubLNmeanSD <- EnumFunGmeanLODBelow(geomedian, LOD, LODBelow)              # Median/geometic mean, Limit of detection, percent below the limit of detection
        EnumSubLNmeanSD <- c(EnumSubLNmeanSD, Total, "No")                            # Adds total sample size to EnumSubLNmeanSD
      } else {
        EnumSubLNmeanSD <- c("High Presence", NA, NA, Total, NA)
      }
    } else if (!is.na(geomedian) &&                                                    # Imputed (geomedian, LOD, %below LOD)
               !is.na(LOD)) {
      AniImput <- EnumSubset$Animal.num[o]                                             # This is a little funky where I am creating variables from the EnumSubset data frame
      PathImput <- EnumSubset$Path.num[o]                                              
      LoDSub <- PrevalenceData[(PrevalenceData$AniNum == AniImput &                     # Then using the variables to create a subset fo the PrevaleanceData data frame.
                                  PrevalenceData$PathNum == PathImput),]
      if (nrow(LoDSub) > 0){
        LODBelow <- 1-LoDSub$Fixed.effects.mean[1]                                     # Then I am pulling from that data frame to find imputed values for the function.
        if (LODBelow < 0.5){
          EnumSubLNmeanSD <- EnumFunGmeanLODBelow(geomedian, LOD, LODBelow)              # Median/geometic mean, Limit of detection, percent below the limit of detection
          EnumSubLNmeanSD <- c(EnumSubLNmeanSD, Total, "Yes")                            # Adds total sample size to EnumSubLNmeanSD
        } else {
          EnumSubLNmeanSD <- c("High Presence", NA, NA, Total, NA)
        }
      } else {next}
    } else if (!is.na(Amean) && 
               !is.na(LOD)) {
      AniImput <- EnumSubset$Animal.num[o]                                             # This is a little funky where I am creating variables from the EnumSubset data frame
      PathImput <- EnumSubset$Path.num[o]                                              
      LoDSub <- PrevalenceData[(PrevalenceData$AniNum == AniImput &                     # Then using the variables to create a subset fo the PrevaleanceData data frame.
                                  PrevalenceData$PathNum == PathImput),]
      if (nrow(LoDSub) > 0){
        LODBelow <- 1-LoDSub$Fixed.effects.mean[1]  
        EnumSubLNmeanSD <- EnumFunMeanLODBelow(Amean, LOD, LODBelow)                   # Mean, Limit of detection, percent below the limit of detection
      } else { next }                                         
    } else { next }
    Impmean <- EnumSubLNmeanSD[2]
    ImpSD <- EnumSubLNmeanSD[3]
    ImpTotDataList <- list()
    for (i in 1:10000) {                                                              # For loop that generates 10000 dataset
      c <- rtrunc("rlnorm", n=LoDPrevBelow, mean = Impmean, sd = ImpSD, lsup = LOD) # Generating the random numbers for the censored values using rtrunc from mc2d. n is the number of values generated, "norm" for normal distribution, mean for mean, sd for standard devation, lsup truncates the data set to below a certain value.
      d <- rtrunc("rlnorm", n=LoDPrev, mean = Impmean, sd = ImpSD, linf = LOD)                                                  # The non censored values
      e <- c(c,d)                                                                     # Creating a complete data set
      f <- fitdist(e, "lnorm")
      ImpTotRow <- c(median(e), 
                     f$estimate[1], 
                     f$estimate[2])                                                   # Creates a vector of the mean and standard deviation of the complete data set
      
      ImpTotDataList[[length(ImpTotDataList)+1]] <- ImpTotRow                         # Adds Row to a list
    }
    ImpTotData <- as.data.frame(do.call(rbind, ImpTotDataList))                     # Creates dataframe from the list
    colnames(ImpTotData) <- c("Median",
                              "Mean",
                              "SD")
    
    EnumSubLNmeanSD <- c(mean(ImpTotData$Median),                  # Output for the meta analysis, Median
                         mean(ImpTotData$Mean),                    # Mean
                         mean(ImpTotData$SD),                      # Standard Deviation
                         EnumSubset$Total[o],                      # Total sample Size
                         "Yes")                                    # If imputation was used to create results 
    EnumOutData <- c(EnumSubset$Studies[o],                                           #This adds data to the output data sheet
                     EnumSubset$class[o],
                     EnumSubset$Pathogen[o],
                     EnumSubset$Conventional..Molecular..or.Immunologic[o],
                     EnumSubset$Enum.Stats[o],
                     EnumSubLNmeanSD,   
                     EnumSubset$Enum.num[o],
                     EnumSubset$Studies.num[o],
                     EnumSubset$Animal.num[o],
                     EnumSubset$Method.num[o],
                     EnumSubset$Path.num[o],
                     EnumSubset$GRADES_Score[o],
                     EnumSubset$USA[o],
                     EnumSubset$Global.Region[o])
    Enumerationlist[[length(Enumerationlist)+1]] <- EnumOutData
  }
}

###Creates dataframe from the list data and adds headers to the datatable.
EnumOut <- as.data.frame(do.call(rbind, Enumerationlist))
colnames(EnumOut) <- c("Study",
                       "Animal",
                       "Pathogen",
                       "Method",
                       "Enumtype",
                       "Median",
                       "lnLognormalMean",
                       "lnLognormalSD",
                       "Sample(n)",
                       "Imputed",
                       "EnumNum",
                       "StudyNum",
                       "AniNum",
                       "MethNum",
                       "PathNum",
                       "GRADE_Score",
                       "USA",
                       "Global.Region")

##### Aggregates the data into a list
# This was not done in the most efficient way because I was trying to retain
# as much data as possible since I was not sure what I would need.

# I removed NA because errors will be generated if there are NA values for
# aggregation also the errors will be caught in the previous tables so I don't
# need to catch them in this one.

# I removed log normal means below 0 because negative values are not suppose
# to be possible for a log normal distrubition.  And represent a convergence
# error in the fitting.

# Run these steps individually Rstudio gets a little wierd when they are all
# run together.

EnumAniName <- c()
EnumMethName <- c()
EnumPathName <- c()
EnumStudyName <- c()
EnumDatasetName <- c()
EnumSampleName <- c()
EnumMeanName <- c()
EnumSDName <- c()
EnumGRADE <- c()
EnumSampleUN <- c()
EnumRubSE <-c()
EnumStudyAgg <- c()


EnumOutcols.num <- c("Median",                                                   # Making all the numbers in the table numeric
                     "lnLognormalMean",
                     "lnLognormalSD",
                     "Sample(n)",
                     "EnumNum", 
                     "StudyNum", 
                     "AniNum", 
                     "PathNum",
                     "MethNum",
                     "GRADE_Score")
EnumOut[EnumOutcols.num] <- sapply(EnumOut[EnumOutcols.num],as.numeric)
EnumOut1 <- EnumOut[!is.na(EnumOut$lnLognormalMean),]                             # Removing all NA data
EnumOut1$Weight <- 10^EnumOut1$GRADE_Score                                        # Adds weight which is based on the GRADES 1000 = high, 100= moderate, 10 = low, 1 = very low


if (nrow(EnumOut1) > 0){
  EnumOut1$Meth.numsub <- as.numeric(as.factor(EnumOut1$Method))                     # Study Number
  EnumOut1$Animal.numsub <- as.numeric(as.factor(EnumOut1$Animal))                   # Animal Number
  EnumOut1$Path.numsub <- as.numeric(as.factor(EnumOut1$Pathogen))                   # Pathogen Number
  
  MethNumEO1 <- max(EnumOut1$Meth.numsub)
  AniNumEO1  <- max(EnumOut1$Animal.numsub)
  PathNumEO1 <- max(EnumOut1$Path.numsub)
  for (i in 1:MethNumEO1) {
    for (j in 1:AniNumEO1) {
      for (k in 1:PathNumEO1) {                                                             # filter by animal
        EnumAniSubset <- EnumOut1[(EnumOut1$Meth.numsub==i &
                                     EnumOut1$Animal.numsub==j &
                                     EnumOut1$Path.numsub==k),]                                # Creates subsets
        if (nrow(EnumAniSubset) > 0){                                                   # Skips subsets that have no data
          EnumSum <- sum(EnumAniSubset$Weight)                                          # The total of the subset weights
          EnumAniSubset$WValues <- EnumAniSubset$Weight/EnumSum*10000                   # The number of weighted values that are going to be randomly generated
          EnumAniSubset$WValues <- round(EnumAniSubset$WValues, digits=0)               # This is rounded since only a whole number value can be generated.  There might be an error here since studies with larger than 10000 difference will have no representation in the final plot.
          EnumAniSubsetSub <- EnumAniSubset[(EnumAniSubset$WValues >0),]                  # Further subsetting to increase speed, since not all of the rows have data that is going to be used since many are knocked out by havng high
          
          EnumAniRow <- nrow(EnumAniSubset)                                             # Number of rows in the subset
          EnumAniRowSub <- nrow(EnumAniSubsetSub) 
          EnumAniSubsetMeans <- c()                                                     # Creates a blank vector to add stuff to add means to that also gets reset after each run
          EnumAniSubsetSD <-c()                                                         # Creates a blank vector to add stuff to add standard deviations to to that also gets reset after each run
          EnumAniParallel <- c()
          EnumAniParallel <- foreach(n=1:10000, .combine = "c", .packages = c("stats","fitdistrplus")) %dopar%{     # 100 replicate, output list
            EnumAniSubsetData <- c()                                                        # Creates a blank vector that randomly generated data is added to.
            for (q in 1:EnumAniRowSub){                                                      # Runs the forloop for each row of the subset
              EnumRandn <-EnumAniSubsetSub$WValues[q]                                        # Sets number of random numbers generated based on distribution
              EnumRandmean <- EnumAniSubsetSub$lnLognormalMean[q]                            # Sets the lognormal mean for the distribution
              EnumRandsd <- EnumAniSubsetSub$lnLognormalSD[q]                                # Sets the lognormal standard deviation for the distribution
              EnumRand <- rlnorm(EnumRandn, meanlog = EnumRandmean, sdlog = EnumRandsd)   # Creates the random numbers
              EnumAniSubsetData <- c(EnumAniSubsetData,EnumRand)                          # Creates a vector that combines all the random numbers generated from the subset
            }
            EnumAniSubsteFit <- fitdist(EnumAniSubsetData, "lnorm")                       # Fits a curve to the data
            EnumAniSubsetMeans <- c(EnumAniSubsteFit$estimate[1], EnumAniSubsteFit$estimate[2])      # Creates a mean with the log mean and standard devation.
          }
          EnumMeanName <- c(EnumMeanName, mean(EnumAniParallel[names(EnumAniParallel)=="meanlog"]))   # Makes a vector of the mean of the means from the replicates
          EnumSDName <- c(EnumSDName, mean(EnumAniParallel[names(EnumAniParallel)=="sdlog"]))
          AggRunMean <- EnumAniParallel[names(EnumAniParallel)=="meanlog"]                  # All the Means calculated for each of the 10000 runs
          AggRunSD <- EnumAniParallel[names(EnumAniParallel)=="sdlog"]                      # All of the Standard deviation for each of the 10000 runs
          AggRunSERow <- AggRunSD/sqrt(10000)                                               # Standard error for each 10000 runs
          AggRunSERow2 <- AggRunSERow^2                                                     # Squared errorfor each row
          AggRunMeanDif2 <- (AggRunMean-mean(AggRunMean))^2                                  # The square difference between the mean of the distribution and mean of all means
          AggRubSE <- sqrt((1/10000)*sum(AggRunSERow2) + (1+(1/10000))*(1/(10000-1))*sum(AggRunMeanDif2)) # Calcualtes the Rubin Standard error
          EnumRubSE <- c(EnumRubSE, AggRubSE)                                               #Collects the Rubin Standard error
          EnumAniName <- c(EnumAniName, EnumAniSubset$Animal[1])                                 # Makes a vector of Animal names 
          EnumMethName <- c(EnumMethName, EnumAniSubset$Method[1])                               # Makes a vector of Method names 
          EnumPathName <- c(EnumPathName, EnumAniSubset$Pathogen[1])                             # Makes a vector of Pathogen names 
          EnumStudyName <- c(EnumStudyName, length(unique(EnumAniSubset$StudyNum)))              # Makes a vector of Number of Studies
          EnumDatasetName <- c(EnumDatasetName, EnumAniRow)                                      # Makes a vector of Number of Datasets
          EnumSampleName <- c(EnumSampleName, sum(EnumAniSubset$`Sample(n)`, na.rm = TRUE))      # Makes a vector of Total Number of Samples
          SampleNA <- 
            if(is.na(sum(EnumAniSubset$`Sample(n)`))==TRUE){"Yes"}else{"No"}
          EnumSampleUN <- c(EnumSampleUN, SampleNA)
          GRADES <- min(EnumAniSubset$GRADE_Score, na.rm = FALSE)                                # Takes the lowest GRADE
          GRADES <- GRADEfunct(GRADES)                                                           # Gives a GRADE rating to the score
          EnumGRADE <- c(EnumGRADE, GRADES)
          EnumStudiesSubUN <- unique(EnumAniSubsetSub$Study)
          EnumStudiesUN <- unique(EnumAniSubset$Study)
          EnumStudiesSub <- paste(unlist(EnumStudiesSubUN), collapse = ", ")
          EnumStudies <- paste(unlist(EnumStudiesUN), collapse = ", ")
          EnumStudyAgg <- c(EnumStudyAgg, EnumStudies)                                                     # Makes a vector of the GRADES
          print(paste(i,MethNumEO1,"_",j,AniNumEO1,"_",k,PathNumEO1))
        } else{next}
      }
    }
  }
}

#####
# Synthesizing total

if (nrow(EnumOut1) > 0){
  for (i in 1:MethNumEO1) {
    for (k in 1:PathNumEO1) {                                                             # filter by animal
      EnumAniSubset <- EnumOut1[(EnumOut1$Meth.numsub==i &
                                   EnumOut1$Path.numsub==k),]                                # Creates subsets
      if (nrow(EnumAniSubset) > 0){                                                   # Skips subsets that have no data
        EnumSum <- sum(EnumAniSubset$Weight)                                          # The total of the subset weights
        EnumAniSubset$WValues <- EnumAniSubset$Weight/EnumSum*10000                   # The number of weighted values that are going to be randomly generated
        EnumAniSubset$WValues <- round(EnumAniSubset$WValues, digits=0)               # This is rounded since only a whole number value can be generated.  There might be an error here since studies with larger than 10000 difference will have no representation in the final plot.
        EnumAniSubsetSub <- EnumAniSubset[(EnumAniSubset$WValues >0),]                  # Further subsetting to increase speed, since not all of the rows have data that is going to be used since many are knocked out by havng high
        
        EnumAniRow <- nrow(EnumAniSubset)                                             # Number of rows in the subset
        EnumAniRowSub <- nrow(EnumAniSubsetSub) 
        EnumAniSubsetMeans <- c()                                                     # Creates a blank vector to add stuff to add means to that also gets reset after each run
        EnumAniSubsetSD <-c()                                                         # Creates a blank vector to add stuff to add standard deviations to to that also gets reset after each run
        EnumAniParallel <- c()
        EnumAniParallel <- foreach(n=1:10000, .combine = "c", .packages = c("stats","fitdistrplus")) %dopar%{     # 100 replicate, output list
          EnumAniSubsetData <- c()                                                        # Creates a blank vector that randomly generated data is added to.
          for (q in 1:EnumAniRowSub){                                                      # Runs the forloop for each row of the subset
            EnumRandn <-EnumAniSubsetSub$WValues[q]                                        # Sets number of random numbers generated based on distribution
            EnumRandmean <- EnumAniSubsetSub$lnLognormalMean[q]                            # Sets the lognormal mean for the distribution
            EnumRandsd <- EnumAniSubsetSub$lnLognormalSD[q]                                # Sets the lognormal standard deviation for the distribution
            EnumRand <- rlnorm(EnumRandn, meanlog = EnumRandmean, sdlog = EnumRandsd)   # Creates the random numbers
            EnumAniSubsetData <- c(EnumAniSubsetData,EnumRand)                          # Creates a vector that combines all the random numbers generated from the subset
          }
          EnumAniSubsteFit <- fitdist(EnumAniSubsetData, "lnorm")                       # Fits a curve to the data
          EnumAniSubsetMeans <- c(EnumAniSubsteFit$estimate[1], EnumAniSubsteFit$estimate[2])      # Creates a mean with the log mean and standard devation.
        }
        EnumMeanName <- c(EnumMeanName, mean(EnumAniParallel[names(EnumAniParallel)=="meanlog"]))   # Makes a vector of the mean of the means from the replicates
        EnumSDName <- c(EnumSDName, mean(EnumAniParallel[names(EnumAniParallel)=="sdlog"]))
        AggRunMean <- EnumAniParallel[names(EnumAniParallel)=="meanlog"]                  # All the Means calculated for each of the 10000 runs
        AggRunSD <- EnumAniParallel[names(EnumAniParallel)=="sdlog"]                      # All of the Standard deviation for each of the 10000 runs
        AggRunSERow <- AggRunSD/sqrt(10000)                                               # Standard error for each 10000 runs
        AggRunSERow2 <- AggRunSERow^2                                                     # Squared errorfor each row
        AggRunMeanDif2 <- (AggRunMean-mean(AggRunMean))^2                                  # The square difference between the mean of the distribution and mean of all means
        AggRubSE <- sqrt((1/10000)*sum(AggRunSERow2) + (1+(1/10000))*(1/(10000-1))*sum(AggRunMeanDif2)) # Calcualtes the Rubin Standard error
        EnumRubSE <- c(EnumRubSE, AggRubSE)                                               #Collects the Rubin Standard error
        EnumAniName <- c(EnumAniName, "Total")                                              # Makes a vector of Animal names 
        EnumMethName <- c(EnumMethName, EnumAniSubset$Method[1])                               # Makes a vector of Method names 
        EnumPathName <- c(EnumPathName, EnumAniSubset$Pathogen[1])                             # Makes a vector of Pathogen names 
        EnumStudyName <- c(EnumStudyName, length(unique(EnumAniSubset$StudyNum)))              # Makes a vector of Number of Studies
        EnumDatasetName <- c(EnumDatasetName, EnumAniRow)                                      # Makes a vector of Number of Datasets
        EnumSampleName <- c(EnumSampleName, sum(EnumAniSubset$`Sample(n)`, na.rm = TRUE))      # Makes a vector of Total Number of Samples
        SampleNA <- 
          if(is.na(sum(EnumAniSubset$`Sample(n)`))==TRUE){"Yes"}else{"No"}
        EnumSampleUN <- c(EnumSampleUN, SampleNA)
        GRADES <- min(EnumAniSubset$GRADE_Score, na.rm = FALSE)                                # Takes the lowest GRADE
        GRADES <- GRADEfunct(GRADES)                                                           # Gives a GRADE rating to the score
        EnumGRADE <- c(EnumGRADE, GRADES)
        EnumStudiesSubUN <- unique(EnumAniSubsetSub$Study)
        EnumStudiesUN <- unique(EnumAniSubset$Study)
        EnumStudiesSub <- paste(unlist(EnumStudiesSubUN), collapse = ", ")
        EnumStudies <- paste(unlist(EnumStudiesUN), collapse = ", ")
        EnumStudyAgg <- c(EnumStudyAgg, EnumStudies)                                                     # Makes a vector of the GRADES
        print(paste(i,MethNumEO1,"_",j,AniNumEO1,"_",k,PathNumEO1))
      } else{next}
    }
  }
}

# Synthesizing by North America

EnumOutNA <- EnumOut1[EnumOut1$Global.Region == 'North America',]

if (nrow(EnumOutNA) > 0){
  MethNumEO1 <- max(EnumOutNA$Meth.numsub, na.rm = TRUE)
  PathNumEO1 <- max(EnumOutNA$Path.numsub, na.rm = TRUE)
  for (i in 1:MethNumEO1) {
    for (k in 1:PathNumEO1) {                                                             # filter by animal
      EnumAniSubset <- EnumOutNA[(EnumOutNA$Meth.numsub==i &
                                    EnumOutNA$Path.numsub==k),]                                # Creates subsets
      if (nrow(EnumAniSubset) > 0){                                                   # Skips subsets that have no data
        EnumSum <- sum(EnumAniSubset$Weight)                                          # The total of the subset weights
        EnumAniSubset$WValues <- EnumAniSubset$Weight/EnumSum*10000                   # The number of weighted values that are going to be randomly generated
        EnumAniSubset$WValues <- round(EnumAniSubset$WValues, digits=0)               # This is rounded since only a whole number value can be generated.  There might be an error here since studies with larger than 10000 difference will have no representation in the final plot.
        EnumAniSubsetSub <- EnumAniSubset[(EnumAniSubset$WValues >0),]                  # Further subsetting to increase speed, since not all of the rows have data that is going to be used since many are knocked out by havng high
        EnumAniRow <- nrow(EnumAniSubset)                                             # Number of rows in the subset
        EnumAniRowSub <- nrow(EnumAniSubsetSub) 
        EnumAniSubsetMeans <- c()                                                     # Creates a blank vector to add stuff to add means to that also gets reset after each run
        EnumAniSubsetSD <-c()                                                         # Creates a blank vector to add stuff to add standard deviations to to that also gets reset after each run
        EnumAniParallel <- c()
        EnumAniParallel <- foreach(n=1:10000, .combine = "c", .packages = c("stats","fitdistrplus")) %dopar%{     # 100 replicate, output list
          EnumAniSubsetData <- c()                                                        # Creates a blank vector that randomly generated data is added to.
          for (q in 1:EnumAniRowSub){                                                      # Runs the forloop for each row of the subset
            EnumRandn <-EnumAniSubsetSub$WValues[q]                                        # Sets number of random numbers generated based on distribution
            EnumRandmean <- EnumAniSubsetSub$lnLognormalMean[q]                            # Sets the lognormal mean for the distribution
            EnumRandsd <- EnumAniSubsetSub$lnLognormalSD[q]                                # Sets the lognormal standard deviation for the distribution
            EnumRand <- rlnorm(EnumRandn, meanlog = EnumRandmean, sdlog = EnumRandsd)   # Creates the random numbers
            EnumAniSubsetData <- c(EnumAniSubsetData,EnumRand)                          # Creates a vector that combines all the random numbers generated from the subset
          }
          EnumAniSubsteFit <- fitdist(EnumAniSubsetData, "lnorm")                       # Fits a curve to the data
          EnumAniSubsetMeans <- c(EnumAniSubsteFit$estimate[1], EnumAniSubsteFit$estimate[2])      # Creates a mean with the log mean and standard devation.
        }
        EnumMeanName <- c(EnumMeanName, mean(EnumAniParallel[names(EnumAniParallel)=="meanlog"]))   # Makes a vector of the mean of the means from the replicates
        EnumSDName <- c(EnumSDName, mean(EnumAniParallel[names(EnumAniParallel)=="sdlog"]))
        AggRunMean <- EnumAniParallel[names(EnumAniParallel)=="meanlog"]                  # All the Means calculated for each of the 10000 runs
        AggRunSD <- EnumAniParallel[names(EnumAniParallel)=="sdlog"]                      # All of the Standard deviation for each of the 10000 runs
        AggRunSERow <- AggRunSD/sqrt(10000)                                               # Standard error for each 10000 runs
        AggRunSERow2 <- AggRunSERow^2                                                     # Squared errorfor each row
        AggRunMeanDif2 <- (AggRunMean-mean(AggRunMean))^2                                  # The square difference between the mean of the distribution and mean of all means
        AggRubSE <- sqrt((1/10000)*sum(AggRunSERow2) + (1+(1/10000))*(1/(10000-1))*sum(AggRunMeanDif2)) # Calcualtes the Rubin Standard error
        EnumRubSE <- c(EnumRubSE, AggRubSE)                                               #Collects the Rubin Standard error
        EnumAniName <- c(EnumAniName, "North America")                                              # Makes a vector of Animal names 
        EnumMethName <- c(EnumMethName, EnumAniSubset$Method[1])                               # Makes a vector of Method names 
        EnumPathName <- c(EnumPathName, EnumAniSubset$Pathogen[1])                             # Makes a vector of Pathogen names 
        EnumStudyName <- c(EnumStudyName, length(unique(EnumAniSubset$StudyNum)))              # Makes a vector of Number of Studies
        EnumDatasetName <- c(EnumDatasetName, EnumAniRow)                                      # Makes a vector of Number of Datasets
        EnumSampleName <- c(EnumSampleName, sum(EnumAniSubset$`Sample(n)`, na.rm = TRUE))      # Makes a vector of Total Number of Samples
        SampleNA <- 
          if(is.na(sum(EnumAniSubset$`Sample(n)`))==TRUE){"Yes"}else{"No"}
        EnumSampleUN <- c(EnumSampleUN, SampleNA)
        GRADES <- min(EnumAniSubset$GRADE_Score, na.rm = FALSE)                                # Takes the lowest GRADE
        GRADES <- GRADEfunct(GRADES)                                                           # Gives a GRADE rating to the score
        EnumGRADE <- c(EnumGRADE, GRADES)
        EnumStudiesSubUN <- unique(EnumAniSubsetSub$Study)
        EnumStudiesUN <- unique(EnumAniSubset$Study)
        EnumStudiesSub <- paste(unlist(EnumStudiesSubUN), collapse = ", ")
        EnumStudies <- paste(unlist(EnumStudiesUN), collapse = ", ")
        EnumStudyAgg <- c(EnumStudyAgg, EnumStudies)                                                     # Makes a vector of the GRADES
        print(paste(i,MethNumEO1,"_",j,AniNumEO1,"_",k,PathNumEO1))
      } else{next}
    }
  }
}


# Next we are going to take all the vectors that were created in the previous
# step and add them to the final list that is going to be created.

EnumTableLenght <- length(EnumMeanName)                                                       # length of the forloop based on the elements in the list created earlier
for (s in 1:EnumTableLenght){                                                                 # the Forlop
  FinalEnumRow <- c(EnumAniName[s],                                                           # Create the row that is added to the list.
                    EnumPathName[s], 
                    EnumMethName[s],
                    EnumGRADE[s],
                    EnumStudyName[s],
                    EnumDatasetName[s],
                    EnumSampleName[s],
                    EnumSampleUN[s],
                    exp(EnumMeanName[s]),                                                      # Creates a row with the median, ln mean, and ln standard deviation
                    EnumMeanName[s], 
                    EnumSDName[s],
                    EnumRubSE[s],
                    EnumStudyAgg[s])                                                            # Adds the name to the rows
  FinalEnumList[[length(FinalEnumList)+1]] <- FinalEnumRow                                    # Adds the row to the list
}


####### Stop here and repeat ######################

# After running previous section run
stopCluster(cl)

# Creates final Enumerated Dataframe

FinalEnumStatsTable <- as.data.frame(do.call(rbind, FinalEnumList))                       # Creates dataframe from the list
colnames(FinalEnumStatsTable) <- c("Subgroup",
                                   "Microbe",
                                   "Method",
                                   "GRADE",
                                   "Report#",
                                   "Datasets",
                                   "Samples",
                                   "Uncounted_Samples",
                                   "Median",
                                   "LN_Mean",
                                   "LN_SD",
                                   "SE",
                                   "Report Name")
FinalEnumStatsTable.num <- c("Report#",                                                    # Changes rows to numeric
                             "Datasets", 
                             "Samples", 
                             "Median",
                             "LN_Mean",
                             "LN_SD",
                             "SE")
FinalEnumStatsTable[FinalEnumStatsTable.num] <- sapply(FinalEnumStatsTable[FinalEnumStatsTable.num],as.numeric)


# Creates final Prevalence Dataframe
FinalPrevalenceTable <- as.data.frame(do.call(rbind, FinalPrevalencelist))
colnames(FinalPrevalenceTable) <- c("Subgroup", 
                                    "Microbe", 
                                    "Method",
                                    "GRADE", 
                                    "Datasets",
                                    "Report#",
                                    "Positive",
                                    "Sample(n)", 
                                    "Fixed.effects.mean",
                                    "Fixed.effects.95low",
                                    "Fixed.effects.95high",
                                    "Random.effects.mean",
                                    "Random.effects.95low",
                                    "Random.effects.95high",
                                    "AniNum",
                                    "MethNum",
                                    "PathNum",
                                    "Report Name") 

FinPrevcol <- ncol(FinalPrevalenceTable)
FinalPrevalenceTable[ , 5:(FinPrevcol-1)] <- apply(FinalPrevalenceTable[ , 5:FinPrevcol], 2,         
                                                   function(x) as.numeric(as.character(x)))
if(!file.exists("Data_out")){
  dir.create("Data_out")
}

write.csv(FinalEnumStatsTable, "Data_out/Total_Enumerated_Dataset.csv", row.names = FALSE)
write.csv(FinalPrevalenceTable, "Data_out/Total_Prevalence_Dataset1.csv", row.names = FALSE)

##### This will create the plots for the articles #####
# Prevalence forest plot
PrevalenceData <- read.csv("Data_out/Total_Prevalence_Dataset1.csv")
EnumData <- read.csv("Data_out/Total_Enumerated_Dataset.csv")

# Next I only want to use the conventional data from both datasets

PrevalenceCon <- PrevalenceData %>% filter(Method == "Conventional")
EnumCon <- EnumData %>% filter(Method == "Conventional")

# Since we can only display a limited amount of information in the script
# I am only going to display data for the largest pathogen categories

PrevalencePath <- PrevalenceCon %>% filter(Microbe == "Fecal coliform"|
                                             Microbe =="Escherichia coli"|
                                             Microbe =="Enterococci"|
                                             Microbe =="Salmonella"|
                                             Microbe =="Giardia"|
                                             Microbe =="Campylobacter"|
                                             Microbe =="Cryptosporidium")

EnumPath <- EnumCon  %>% filter(Microbe == "Fecal coliform"|
                                  Microbe =="Escherichia coli"|
                                  Microbe =="Enterococci"|
                                  Microbe =="Salmonella"|
                                  Microbe =="Giardia"|
                                  Microbe =="Campylobacter"|
                                  Microbe =="Cryptosporidium")

# I am going to remove plots with one report since there is not enought data to give
# an unbiased meta analysis also it is not providing any new information from viewing
# the original papers and it already represented in the total and maybe the USA data

PrevalencePath <- PrevalencePath %>% filter(Report. > 1 & !Subgroup=="USA")


# Removes studies with less than 1 study and USA data
EnumPath1 <- EnumPath %>% filter(Report. > 1 & !Subgroup=="USA")
# Makes another data frame with North america to keep north america with <1 sample
EnumPath2 <- EnumPath %>% filter(Subgroup == "North America")
# Joings the two table
EnumPath3 <- full_join(EnumPath1,EnumPath2)
EnumPath3 <- arrange(EnumPath3, EnumPath3$Microbe)
EnumPath <- EnumPath3

# Since I am adding the total and the USA data I am just going to keep everything at the same
# everything the same since it is already arranged with the total and USA last and I am
# not going to do the tornado plot thing of arranging from highest to lowest since the
# gets weird.


ForestSubSort <- arrange(PrevalencePath, PrevalencePath$Microbe)
# Structures the data in a way that forestplot can use to display.

# Changes prevalence from a ratio to a percentage
ForestSubSort$Random.effects.mean <- 100 * ForestSubSort$Random.effects.mean
ForestSubSort$Random.effects.95low <- 100 * ForestSubSort$Random.effects.95low
ForestSubSort$Random.effects.95high <- 100 * ForestSubSort$Random.effects.95high


anum <- nrow(ForestSubSort)+1
ForestDat <- structure(list(mean  = c(NA, ForestSubSort$Random.effects.mean), 
                            lower = c(NA, ForestSubSort$Random.effects.95low),
                            upper = c(NA, ForestSubSort$Random.effects.95high)),
                       .Names = c("mean", "lower", "upper"), 
                       row.names = c(NA, -anum),
                       class = "data.frame")
Foresttext <- cbind(c("Microbe", ForestSubSort$Microbe),
                    c("Subgroup", ForestSubSort$Subgroup),
                    c("Report#", ForestSubSort$Report.),
                    c("Positive", ForestSubSort$Positive),
                    c("Sample(n)", ForestSubSort$Sample.n.),
                    c("GRADE", ForestSubSort$GRADE))
# Sets the Dimensions of the graph
fheight <- 15*anum+100

# Line assignments for horizontal lines which you assign by changing the number
# in fp_add_lines(h_) The first number is always 2 since that is the line under the
# the header.
line1 <- 2+nrow(ForestSubSort[ForestSubSort$Microbe == "Campylobacter",])
line2 <- line1+nrow(ForestSubSort[ForestSubSort$Microbe == "Cryptosporidium",])
line3 <- line2+nrow(ForestSubSort[ForestSubSort$Microbe == "Enterococci",])
line4 <- line3+nrow(ForestSubSort[ForestSubSort$Microbe == "Escherichia coli",])
line5 <- line4+nrow(ForestSubSort[ForestSubSort$Microbe == "Fecal coliform",])
line6 <- line5+nrow(ForestSubSort[ForestSubSort$Microbe == "Giardia",])
line7 <- line6+nrow(ForestSubSort[ForestSubSort$Microbe == "Salmonella",])

# Creates the Forest plot 
png(paste('Forest/Journal_forestplot.png'),
    # These adjust the size of the figure.
    width = 850, height = fheight , units = "px") 
print(ForestDat |> 
        forestplot(labeltext = Foresttext, 
                   is.summary = c(rep(TRUE, 1), rep(FALSE, anum - 1)),
                   boxsize = .1,
                   xticks = c(0,20,40,60,80,100),
                   xlab = "% Prevalence",
                   # title = "Prevalence of pathogens in animal feces",
                   graphwidth = unit (5, "cm"),) |>
        fp_add_lines(h_2 = gpar(lty = 1),
                     h_6 = gpar(lty = 1),
                     h_10 = gpar(lty = 1),
                     h_13 = gpar(lty = 1),
                     h_17 = gpar(lty = 1),
                     h_20 = gpar(lty = 1),
                     h_24 = gpar(lty = 1),
                     h_30 = gpar(lty = 1))|>
        fp_set_style(box = "royalblue",
                     line = "lightblue",
                     summary = "royalblue",
                     txt_gp = fpTxtGp(ticks = gpar(fontfamily = "",  cex = 1),
                                      xlab = gpar(fontfamily = "", cex =1))
                     
        )
)

dev.off()

# Creating a boxplot using the data that was generated before.
# EnumPath

# The first thing that I need to do is change the data so that it is in a format 
# that a boxplot can be made out of which I am going to do by creating a thousand
# random points for each distribution.

# Creates the data frame for the boxplot

AniboxC <- c()
PathboxC <- c()
EnumRandC <- c()

EnumBNum <- nrow(EnumPath)
nruns <- 1000
for (i in 1:EnumBNum) {                                                           # Go row by row
  Enummeanbox <- EnumPath$LN_Mean[i]                                              # Sets the lognormal mean for the distribution
  EnumSDbox <- EnumPath$LN_SD[i]                                                  # Sets the lognormal standard deviation for the distribution
  EnumRand <- rlnorm(nruns, meanlog = Enummeanbox, sdlog = EnumSDbox)             # Creates the random numbers
  Anibox <- rep(EnumPath$Subgroup[i],nruns)
  Pathbox <-rep(EnumPath$Microbe[i],nruns)
  AniboxC <- c(AniboxC, Anibox)
  PathboxC <- c(PathboxC, Pathbox)
  EnumRandC <- c(EnumRandC,EnumRand)
} 

EnumboxDF <- data.frame(PathboxC, AniboxC, EnumRandC)

#Creates the boxplot


png(paste('Forest/Journal_EnumBoxplot.png'),
    # These adjust the size of the figure.
    width = 700, height = 700 , units = "px") 
print(Enumboxplot <- ggplot(EnumboxDF, aes(x=PathboxC, y=EnumRandC, fill=AniboxC)) + 
        geom_boxplot(position = position_dodge2(preserve = "single")) +
        geom_hline(yintercept=1)+
        facet_wrap(~PathboxC, scale="free") + 
        scale_y_continuous(trans='log10')+ 
        xlab("Subgroup")+
        ylab("Microbe/cyst/oocyst concentration per gram of wet feces")+
        theme_classic()+
        theme(legend.position = "bottom",
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              text = element_text(size = 15),
              legend.title = element_blank())+
        scale_fill_grey(start = 0, end = 1)+
        labs(color=NULL)
)

dev.off()

# Creates a table of the boxplot data
write.csv(EnumPath, "Data_out/Enumerated_boxplot.csv", row.names = FALSE)

# Converts to nomral median and standard devation

Enummeanboxtable <- aggregate(EnumRandC~PathboxC+AniboxC, median, data=EnumboxDF)
EnumIQRboxtable <- aggregate(EnumRandC~PathboxC+AniboxC, IQR, data=EnumboxDF)
Enumboxtable <- cbind(Enummeanboxtable,EnumIQRboxtable$EnumRandC)
colnames(Enumboxtable) <- c("Microbe", "Subgroup", "Mean","IQR")

EnumDis <- EnumPath
EnumDisp <- EnumDis %>% left_join(Enumboxtable)
EnumDispTable <- EnumDisp[,c("Microbe", "Subgroup", "Report.", "Samples", "GRADE", "LN_Mean", "LN_SD", "Median", "IQR")]
write.csv(EnumDispTable, "Data_out/Enumerated_box_display.csv", row.names = FALSE)

# Saves the workspace
save.image('workspaceNA.RData')
