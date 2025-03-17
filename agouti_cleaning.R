# CLEANING AGOUTI OUTPUT
require(dplyr)
require(lubridate)
require(stringr)
require(ggplot2)
require(DescTools)
require(magrittr)
library(tidyverse)

# check if working directory is set okay
setwd("~/GitHub/camtrap_coiba") # change to Git 

#### GENERAL CLEANING ####
### SOME NOTES ON UNDERSTANDING THE OUTPUT FROM AGOUTI:
## This data is downloaded by pressing "export" on the agouti website
## It gives the following things: all sequences in all deployments that are uploaded to agouti, coded or not
## For coded sequences, if the sequence was coded as "blank" it is not included in this output
## For not finished deployments, all the uncoded sequences are also included (so this also includes 'blank' sequences that haven't been marked as such yet)
## Agouti also provides a multimedia file and deployment keys

## Update March 2025. Agouti changed their output, and they now changed "sequenceID" to "eventID". They also have seqend and seqstart in as event end and event start
# individual ID is now real name not key
# camera set up is now "setup" or nothing instead of TRUE/FALSE


### TO DO's:
## - for known juvenile individuals, need to code which year/deployment they shift to (sub)adult
## - for affiliative interactions, how do we know who is grooming (or being friendly) with whom if we've coded more than 2 individuals as grooming?
## - can also look at autocorrelation by for instance making a day-hour variable and doing ftable of day-hour and all Jicaron tool use locations. 
#   then see if you have activity at several camera traps within the same hour (hour might be too large of a timescale). Need spatial depth for this. 

# open Agouti output file (observations) that you have downloaded from the agouti website. Use most recent version
agoutigross <- read.csv("agouti_output/coiba-national-park-tool-use-20250228090345/observations.csv", header = TRUE)
# add sequenceID variable (replace eventID)
agoutigross$sequenceID <- agoutigross$eventID
# remake camera set up
agoutigross$cameraSetup <- ifelse(agoutigross$cameraSetupType == "setup", "True", "False")

# open the associated deployment keys (also downloaded from agouti.eu)
depl_keys <- read.csv("agouti_output/coiba-national-park-tool-use-20250228090345/deployments.csv", header = TRUE)
depl_keys$tags <- depl_keys$deploymentTags

# filter out test deployments/not relevant ones (so create variable to filter test ones)
## THIS WILL NEED TO BE MORE FINETUNED LATER. THERE ARE SOME TRIAL/WRONG DATA ON THERE THAT MAY NOT BE CAPTURED NOW.
# also have "caution" flag for deployments that might need some modification/fixing. 
depl_keys$flag <- ifelse(grepl("Flag", depl_keys$deploymentTags) | grepl("Test", depl_keys$deploymentTags) | depl_keys$deploymentTags == "", 1, 0 )
agoutigross <- left_join(agoutigross, depl_keys, "deploymentID")
agoutigross <- agoutigross[agoutigross$flag == 0,]

#### Timestamps & Sequence Info #####
# turn sequence start time into POSIXct format and replace T with a space
agoutigross$time <- str_replace(agoutigross$eventStart, "T", " ")
agoutigross$time <- as.POSIXct(agoutigross$time, tz = "America/Panama", format = "%Y-%m-%d %H:%M:%S")

# identify and correct wrong timestamps
# open the multimedia csv containing the correct timestamps (also from agouti)
multimedia <- read.csv("agouti_output/coiba-national-park-tool-use-20250228090345/media.csv", header = TRUE)
# now the sequenceID is buried in the mediaComments
multimedia$sequenceID <- sapply(str_split(multimedia$mediaComments, ":"), '[', 2)

# have both timestamps we entered incorrectly (e.g. 1970) and those that shifted 5 hours by accident
multimedia$time <- str_replace(multimedia$timestamp, "T", " ")
multimedia$time <- as.POSIXct(multimedia$time, tz = "America/Panama", format = "%Y-%m-%d %H:%M:%S")

# obtain correct timestamp from file name
multimedia$time_file <- paste(sapply(str_split(multimedia$fileName, "__"), "[", 2), sapply(str_split(multimedia$fileName, "__"), "[", 3), " ")
# mp4 and jpg strings are extracted differently
multimedia$time_file[which(str_detect(multimedia$time_file, "mp4"))] <- substr(multimedia$time_file[which(str_detect(multimedia$time_file, "mp4"))], 1,
                                                                               nchar(multimedia$time_file[which(str_detect(multimedia$time_file, "mp4"))])-6)

multimedia$time_file <- sapply(str_split(multimedia$time_file, "\\("), '[', 1)
multimedia$time_file <- as.POSIXct(multimedia$time_file, tz = "America/Panama", format = "%Y-%m-%d %H-%M-%S")

# now flag those where the timestamp does not match the file time
multimedia$timeflag <- ifelse(multimedia$time == multimedia$time_file, 0, 1)
ftable(multimedia$deploymentID, multimedia$timeflag)

# the file time is (nearly) always correct, so can just use that one for the one's that are flagged
# covers both the entry mistakes and the ones that drifted 5 hours
# create timestamp_correct 
multimedia$timestamp_correct <- multimedia$timestamp
multimedia$timestamp_correct[which(multimedia$timeflag == 1)] <- as.character(multimedia$time_file[which(multimedia$timeflag == 1)])

# EXCEPTIONS: In at least one deployment (CEBUS-09-R3) the filetime is one year off. Also the R9 deployment has this consistently
# deployment key of CEBUS-09-R3: 09bcd191-78ec-4ded-a438-cffaa14fa55a. Fix this
# deployment CEBUS-03-R8: 82cd9f63-8756-4987-bc32-8082090996a3 is 12 hours late
multimedia$timestamp_correct[which(multimedia$deploymentID == "09bcd191-78ec-4ded-a438-cffaa14fa55a")] <- gsub("2017", "2018", multimedia$timestamp_correct[which(multimedia$deploymentID == "09bcd191-78ec-4ded-a438-cffaa14fa55a")]) 
multimedia$timestamp_correct <- multimedia$timestamp_correct %>%
  as.character(.) %>%
  str_replace(., "T", " ") %>%
  as.POSIXct(., tz = "America/Panama", format = "%Y-%m-%d %H:%M:%S")
multimedia$timestamp_correct[which(multimedia$deploymentID == "82cd9f63-8756-4987-bc32-8082090996a3")] <- multimedia$timestamp_correct[which(multimedia$deploymentID == "82cd9f63-8756-4987-bc32-8082090996a3")] - 43200

# create sequence_start, sequence_end and duration column (so this is per sequence)
startseq <- aggregate(x = list(seq_start = multimedia$timestamp_correct), by = list(sequenceID = multimedia$sequenceID), FUN = min)
endseq <- aggregate(x = list(seq_end = multimedia$timestamp_correct), by = list(sequenceID = multimedia$sequenceID), FUN = max)

multimedia <- left_join(multimedia, startseq, "sequenceID")
multimedia <- left_join(multimedia, endseq, "sequenceID")
multimedia$seq_lengthstill <-as.numeric(difftime(multimedia$seq_end,multimedia$seq_start,units="secs"))

#### Add temperature and video sequence length ####
# first run script "exiftempseq_cleaning.R" in the exiftool_metadata folder
## clean filename in multimedia agouti (take off the upload time, which is always 15 characters)
multimedia$FileName <- substring(multimedia$fileName, 16)

multimedia <- left_join(multimedia, exifstillsclean, "FileName")
multimedia <- left_join(multimedia, exifvideosclean, "FileName")

multimedia$temperature <- ifelse(is.na(multimedia$temperature.x), multimedia$temperature.y, multimedia$temperature.x)
multimedia$seq_length <- ifelse((multimedia$fileMediatype == "video/mp4"), multimedia$seq_lengthv, multimedia$seq_lengthstill) 
multimedia$mediatype <- ifelse((multimedia$fileMediatype == "video/mp4"), "video", "still") 

# left join does not work if the y dataframe (so multimedia2) has duplicated values in the column you're matching by (so sequenceID). Need to drop duplicates before the left_join
multimedia2 <- multimedia[,c("sequenceID", "seq_start", "seq_end", "seq_length", "temperature", "mediatype", "captureMethod")]
multimedia2 <- multimedia2[!duplicated(multimedia2$sequenceID),]
agoutigross <- left_join(agoutigross, multimedia2, "sequenceID")

# add correct deployment start and end time (can still double check, now took the minimum and maximum sequence time within deployment id)
# need to do this per camera (because camera could stop before deployment end)
# rename SURVEY-CEBUS-24-01 to CEBUS-04
agoutigross$locationName[which(agoutigross$locationName == "SURVEY-CEBUS-24-01")] <- "CEBUS-04"
# create unique variable (like deployment ID) that is location name + tag
agoutigross$uniqueloctag <- paste(agoutigross$locationName, agoutigross$tags, sep = "-")
# remove "caution" from the name and instead create a caution column
agoutigross$caution <- ifelse(str_detect(agoutigross$uniqueloctag, "-Caution \\| ") == TRUE, 1, 0)
agoutigross$uniqueloctag <- str_replace(agoutigross$uniqueloctag, "-Caution \\| R", "-")
startdep <- aggregate(x = list(dep_start = agoutigross$seq_start), by = list(uniqueloctag = agoutigross$uniqueloctag), FUN = min)
enddep <- aggregate(x = list(dep_end = agoutigross$seq_start), by = list(uniqueloctag = agoutigross$uniqueloctag), FUN = max)

agoutigross <- left_join(agoutigross, startdep, "uniqueloctag")
agoutigross <- left_join(agoutigross, enddep, "uniqueloctag")
# add deployment length in hours, this is an exposure
agoutigross$dep_length_hours <-as.numeric(difftime(agoutigross$dep_end,agoutigross$dep_start,units="hours"))

## Corrections: TU-155 was turned on two days prior to deployment, real deployment start is 2022-07-15 11:12:50
agoutigross$dep_start[which(agoutigross$locationName == "TU-155")] <- "2022-07-15 11:12:50 EST"

#### ADDING EXTRA VARIABLES ####
# make month column
agoutigross$month <- month(agoutigross$seq_start)

# make wet season may-nov and dry season dec-april
agoutigross$season <- ifelse(agoutigross$month == 12 | agoutigross$month == 1 | agoutigross$month == 2 | agoutigross$month == 3 | 
                         agoutigross$month == 4, "Dry", "Wet") 

# pull island location and tool use/non tool use from coiba_camtrap_ids_gps.csv
# take long-lat from coiba_camtrap_ids_gps instead of what we entered on agouti
setwd("~/GitHub/camtrap_coiba") # change to Git if need be
deployment_info <- read.csv("coiba_camtrap_ids_gps.csv")
deployment_info$locationName <- deployment_info$camera_id
deployment_info <- deployment_info[!deployment_info$locationName == "SURVEY-CEBUS-24-01",] 

# drop columns we don't want to attach
deployment_info2 <- deployment_info[, !names(deployment_info) %in% c("camera_id", "number")]

# FYI, long/lats in agouti are sometimes different from ones in csv. assume that csv is correct
# compare longitudes and latitudes from agouti with ones from csv
test <- agoutigross[!duplicated(agoutigross$locationName), c("locationName", "longitude", "latitude")]
test <- left_join(test, deployment_info2, "locationName")
test <- test[,c("locationName", "longitude.x", "longitude.y", "latitude.x", "latitude.y")]
test[,2:5] <- round(test[,2:5], 5)

test$latitudematch <- ifelse(test$latitude.x == test$latitude.y, 1, 0)
test$longitudematch <- ifelse(test$longitude.x == test$longitude.y, 1, 0)
test$mismatch <- ifelse(test$latitudematch == 0 | test$longitudematch == 0, 1, 0)
mismatches <- test[test$latitudematch == 0 | test$longitudematch == 0,]

#ggplot(test[test$mismatch == 1,], aes(label = locationName)) + geom_point(aes(x= longitude.x, y = latitude.x), col = "black") + 
# geom_text(aes(x = longitude.x, y = latitude.x), hjust = 0, vjust = 0, col = "black") +
#  coord_cartesian(xlim = c(-81.8225, -81.816), ylim = c(7.266,7.273)) +
#  geom_point(aes(x = longitude.y, y = latitude.y), col = "red", alpha = 0.5) + geom_text(aes(x = longitude.y, y = latitude.y), hjust = 0, vjust = 0, col = "red") 

# overwrite those in agoutigross with ones from csv
agoutigross <- agoutigross[,!names(agoutigross) %in% c("longitude", "latitude")]

agoutigross <- left_join(agoutigross, deployment_info2, "locationName")

#### Capuchin detections #####
# create column per sequence if a capuchin was present 1/0 and if capuchins, how many in the sequence (a count)
agoutigross$capuchin <- ifelse(agoutigross$scientificName == "Cebus imitator", 1, 0)

# make blank sex or life stage put as unknown when there is a capuchin
agoutigross$sex[which(agoutigross$capuchin == TRUE & agoutigross$sex == "")] <- "unknown"
agoutigross$lifeStage[which(agoutigross$capuchin == TRUE & agoutigross$lifeStage == "")] <- "unknown"

## for ID'd juvenile and subadult capuchins, code when their ages change 
# for now only did this for the ones that we have seen at various lifestages
# SPT, RIC, INK, LAR, YOD (not coded later yet)
# need to find their agouti codes
# with march updates agouti key codes are gone and IDs are in main output
# SPT = 66437816-fdf4-4fde-9dbd-8591bfe8cb1c
# LAR = c6114eb1-b53e-4494-9565-0dd621ff88b9
# YOD = cbf8f2a0-f1ed-425b-9a71-8b59227b5546
# RIC = 108f3544-b14a-43d9-a3a0-b805f35ae2db
# INK =  db68538d-ff60-433a-bf94-e9a886a6d23c
agoutigross$lifeStage <- ifelse(agoutigross$individualID == "SPT (Spot)" & year(agoutigross$seq_start) > 2018 & year(agoutigross$seq_start) < 2022, "subadult", 
                                ifelse(agoutigross$individualID == "SPT (Spot)" & year(agoutigross$seq_start) > 2021, "adult",
                                       ifelse(agoutigross$individualID == "INK (Ink)" & year(agoutigross$seq_start) > 2020, "adult", 
                                              ifelse(agoutigross$individualID == "LAR (Larry)" & year(agoutigross$seq_start) > 2019 & year(agoutigross$seq_start) < 2022, "subadult",
                                                     ifelse(agoutigross$individualID == "LAR (Larry)" & year(agoutigross$seq_start) > 2021, "adult", agoutigross$lifeStage)))))

## CHECK THIS!!

# make agesex variable that is paste of sex and life stage
agoutigross$agesex <- ifelse(agoutigross$capuchin == TRUE, 
                             paste(agoutigross$lifeStage, agoutigross$sex, sep = " "), NA)
as.matrix(ftable(agoutigross$agesex))
agoutigross$agesexF <- as.factor(agoutigross$agesex)

agoutigross_cap <- agoutigross[agoutigross$capuchin == 1, ]
cap_numbers <- agoutigross_cap %>%
  dplyr::count(sequenceID)

agoutigross <- left_join(agoutigross, cap_numbers, "sequenceID")
# replace NAs with 0 for the capuchin count
agoutigross$n[is.na(agoutigross$n)] <- 0

## First add ID string instead of code
# no longer necessary with March 2025 update
#IDS <- read.csv("agouti_output/individual_age_sex_species_202306131450.csv")
#IDS <- IDS[,c("id", "name")]

#agoutigross <- left_join(agoutigross, IDS, by = c("individualID" = "id"))

##### Age/sex classes #####
# make frequency table of agesex per sequence
cap_agesex <- as.data.frame(as.matrix(ftable(agoutigross_cap$sequenceID, agoutigross_cap$agesex)))
colnames(cap_agesex) <- c("nAF", "nAM", "nAU", "nJF", "nJM", "nJU", "nSF", "nSM", "nSU", 
                          "nUF", "nUM", "nUU")
cap_agesex$sequenceID <- rownames(cap_agesex)

# number of adult females carrying infant
cap_agesex_infant <- as.data.frame(as.matrix(ftable(agoutigross_cap$sequenceID[which(str_detect(agoutigross_cap$behavior, "Infant Care") == TRUE)], 
                                             agoutigross_cap$agesex[which(str_detect(agoutigross_cap$behavior, "Infant Care") == TRUE)])))
colnames(cap_agesex_infant) <- c("nAF_infant", "nAM_infant", "nAU_infant", "nJF_infant", "nJM_infant", "nJU_infant", "nSF_infant", "nSM_infant", "nSU_infant", 
                          "nUF_infant", "nUM_infant", "nUU_infant")
cap_agesex_infant$sequenceID <- rownames(cap_agesex_infant)
## there is a bit of an error here, where sometimes for the dorsal infant themselves "infant care" was coded, so they get counted as an unknown juvenile carrying an infant
## THIS IS NOT YET FIXED, SEEMS TO ALSO BE FOR CASES WHERE INFANT WAS NOT CODED AS INFANT CARE!! 

# number of neck infants
cap_neck_infant <- as.data.frame(ftable(agoutigross$sequenceID[which(str_detect(agoutigross$observationComments, "neck infant") == TRUE)]))
colnames(cap_neck_infant) <- c("sequenceID", "n_neckinfant")

# add this to agoutigross dataframe
agoutigross <- left_join(agoutigross, cap_agesex, "sequenceID")
agoutigross <- left_join(agoutigross, cap_agesex_infant, "sequenceID")
agoutigross <- left_join(agoutigross, cap_neck_infant, "sequenceID")

# make new colum for total number of males, total number of females,
# total adults, total subadults, total juveniles per sequence
agoutigross$nAdult <- agoutigross$nAF + agoutigross$nAM + agoutigross$nAU
agoutigross$nSubadult <- agoutigross$nSF + agoutigross$nSM + agoutigross$nSU
agoutigross$nJuvenile <- agoutigross$nJF + agoutigross$nJM + agoutigross$nJU
agoutigross$nFemales <- agoutigross$nAF + agoutigross$nSF + agoutigross$nJF + agoutigross$nUF
agoutigross$nMales <- agoutigross$nAM + agoutigross$nSM + agoutigross$nJM + agoutigross$nUM

##### Inspections #####
# add nr or capuchins inspecting camera trap
agoutigross_capinsp <- agoutigross[agoutigross$capuchin == 1 & str_detect(agoutigross$behavior, "Inspecting"), ]
# nr of adults inspecting
cap_numbers_insp_ad <- agoutigross_capinsp[agoutigross_capinsp$lifeStage == "adult",] %>%
  dplyr::count(sequenceID)
colnames(cap_numbers_insp_ad) <- c("sequenceID", "n_inspect_adult")
head(cap_numbers_insp_ad)
# nr of non-adults (juveniles + subadults)
cap_numbers_insp_nonad <- agoutigross_capinsp[agoutigross_capinsp$lifeStage == "juvenile" | agoutigross_capinsp$lifeStage == "subadult",] %>%
  dplyr::count(sequenceID)
colnames(cap_numbers_insp_nonad) <- c("sequenceID", "n_inspect_nonadult")
# nr of unknowns
cap_numbers_insp_unknown <- agoutigross_capinsp[agoutigross_capinsp$lifeStage == "unknown",] %>%
  dplyr::count(sequenceID)
colnames(cap_numbers_insp_unknown) <- c("sequenceID", "n_inspect_unknown")

agoutigross <- left_join(agoutigross, cap_numbers_insp_ad, "sequenceID")
agoutigross <- left_join(agoutigross, cap_numbers_insp_nonad, "sequenceID")
agoutigross <- left_join(agoutigross, cap_numbers_insp_unknown, "sequenceID")

##### Displacements #####
# info we'd like to have for analyses and need to extract now
# very broadly could do whether displacement occurs in a sequence 1/0 and then see if there are adult females present
agoutigross_disp <- agoutigross[str_detect(agoutigross$observationComments, "displace"), c("agesex", "lifeStage", "sex", "sequenceID", "observationComments")]

agoutigross$displacement <- ifelse(agoutigross$sequenceID %in% agoutigross_disp$sequenceID, 1, 0)

# agesex of displacer based on comments
# agesex of displacee based on comments
agoutigross_disp$victim_as <- ifelse(str_detect(agoutigross_disp$observationComments, "displaced"), agoutigross_disp$agesex, NA)
agoutigross_disp$aggressor_as <- ifelse(str_detect(agoutigross_disp$observationComments, "displaced") == FALSE, agoutigross_disp$agesex, NA)
victims <- agoutigross_disp[!is.na(agoutigross_disp$victim_as),c("sequenceID", "victim_as")]
aggressors <- agoutigross_disp[!is.na(agoutigross_disp$aggressor_as),c("sequenceID", "aggressor_as")]
displacements <- full_join(victims, aggressors, "sequenceID")
# so don't have all victims/aggressors but something is better than nothing!

displacements$victim_as[is.na(displacements$victim_as) | displacements$victim_as == "unknown unknown"] <- "unknown"
displacements$aggressor_as[is.na(displacements$aggressor_as)| displacements$aggressor_as == "unknown unknown"] <- "unknown"

agoutigross <- left_join(agoutigross, displacements, "sequenceID")

ftable(agoutigross$victim_as)
# so get some information from this! (this is still multiple times per sequence, need to drop down to sequence level ofc)

#### Foraging items ####
# what is being foraged on per sequence
# check if one individual forages on more things
multipleitems <- agoutigross[which(str_count(agoutigross$behavior, "F:")> 1),]
# not that many, and usually they eat the same thing. So I think it's fine to only record one item per individual
# follows the order specified below

# what they are foraging using tools
agoutigross$tool_item <- ifelse(str_detect(agoutigross$behavior, "TAF: Palm"), "palm", 
                           ifelse(str_detect(agoutigross$behavior, "TAF: Coconut"), "coconut", 
                                  ifelse(str_detect(agoutigross$behavior, "TAF: Embedded"), "insect",
                                         ifelse(str_detect(agoutigross$behavior, "TAF: Halloween"), "hwcrab",
                                                ifelse(str_detect(agoutigross$behavior, "TAF: Hermit"), "hecrab",
                                                       ifelse(str_detect(agoutigross$behavior, "TAF: Snail"), "snail",
                                                              ifelse(str_detect(agoutigross$behavior, "TAF: Almendra"), "almendra",
                                                                     ifelse(str_detect(agoutigross$behavior, "TAF: Other"), "other", 
                                                                            ifelse(str_detect(agoutigross$behavior, "TAF: Unknown"), "unknown", NA)))))))))                                        



ftable(agoutigross$tool_item)

# what they are foraging on without tools
agoutigross$normal_item <- ifelse(str_detect(agoutigross$behavior,  "F: Fruit"), "fruit", 
                                  ifelse(str_detect(agoutigross$behavior, "F: Coco flesh|F: Coco water"), "coconut", 
                                         ifelse(str_detect(agoutigross$behavior, "F: Insect"), "insect",
                                                                     ifelse(str_detect(agoutigross$behavior, "F: Almendra flesh"), "almendra",
                                                                            ifelse(str_detect(agoutigross$behavior, "\\bF: Other"), "other", 
                                                                                   ifelse(str_detect(agoutigross$behavior, "\\bF: Unknown"), "unknown", NA))))))                                      

ftable(agoutigross$normal_item)

# if 'other', see if we can get information from the comments
agoutigross$comment_item <- NA
agoutigross$comment_item[which(str_detect(agoutigross$behavior, "Other"))] <- 
  agoutigross$observationComments[which(str_detect(agoutigross$behavior, "Other"))]

# include at least fruit and crab and snail
agoutigross$normal_item[which(agoutigross$normal_item == "other")] <- 
  ifelse(str_detect(agoutigross$comment_item[which(agoutigross$normal_item == "other")], "crab"), "crab",
         ifelse(str_detect(agoutigross$comment_item[which(agoutigross$normal_item == "other")], "fruit|berr"), "fruit",
                ifelse(str_detect(agoutigross$comment_item[which(agoutigross$normal_item == "other")], "snail"), "snail", "other")))

# also for tool using 
agoutigross$comment_item[which(agoutigross$tool_item == "other")]

agoutigross$tool_item[which(agoutigross$tool_item == "other")] <- 
  ifelse(str_detect(agoutigross$comment_item[which(agoutigross$tool_item == "other")], "Halloween crab"), "hwcrab",
         ifelse(str_detect(agoutigross$comment_item[which(agoutigross$tool_item == "other")], "fruit"), "fruit",
                ifelse(str_detect(agoutigross$comment_item[which(agoutigross$tool_item == "other")], "snail"), "snail", "other")))

# whether tool-using occurred in the sequence or not
agoutigross$tooluse <- str_detect(agoutigross$behavior, "TAF") # now just takes all types of tool use, also unknown
agoutigross_tools <- agoutigross[agoutigross$tooluse == TRUE,]
# amount of individuals using tools per sequence
tooluse_count <- agoutigross_tools %>%
  dplyr::count(sequenceID)
colnames(tooluse_count)[2] <- "n_tooluse"

agoutigross <- left_join(agoutigross, tooluse_count, "sequenceID")
# replace NAs with 0 for the tool use
agoutigross$n_tooluse[is.na(agoutigross$n_tooluse)] <- 0
# make sure only the tool use sequences have items (don't want to get the normal foraging items)
# not really necessary but a good check
agoutigross$tool_item[agoutigross$tooluse == "FALSE"] <- NA

# also do tool using identity (adult/subadult/juvenile or even age/sex (though no females really so eh))
head(agoutigross_tools)

# make frequency table of agesex per sequence
cap_toolsagesex <- as.data.frame(as.matrix(ftable(agoutigross_tools$sequenceID, agoutigross_tools$agesexF)))

colnames(cap_toolsagesex) <- c("tu_nAF", "tu_nAM", "tu_nAU", "tu_nJF", "tu_nJM", "tu_nJU", "tu_nSF", "tu_nSM", "tu_nSU", 
                          "tu_nUF", "tu_nUM", "tu_nUU")
cap_toolsagesex$sequenceID <- rownames(cap_toolsagesex)

# add this to agoutigross dataframe
agoutigross <- left_join(agoutigross, cap_toolsagesex, "sequenceID")

# make new colum for total adults, total subadults, total juveniles per sequence
agoutigross$tu_nAdult <- agoutigross$tu_nAF + agoutigross$tu_nAM + agoutigross$tu_nAU
agoutigross$tu_nSubadult <- agoutigross$tu_nSF + agoutigross$tu_nSM + agoutigross$tu_nSU
agoutigross$tu_nJuvenile <- agoutigross$tu_nJF + agoutigross$tu_nJM + agoutigross$tu_nJU

### APPROACH 1: get sequence-level variable of what is mostly being processed in that sequence
# also get number of capuchins consuming each item per sequence 
# tool items
tool_items <- as.data.frame(as.matrix(ftable(agoutigross_tools$sequenceID, agoutigross_tools$tool_item)))
# set unknown tool use to 1 for all rows, so that when there are two things being processed (e.g. one unknown, one almendra) you keep the almendra
# this is crude solution but works
tool_items$unknown <- 1
tool_items$seq_toolitem <- colnames(tool_items)[max.col(tool_items,ties.method="first")]  
# now gets the most frequently processed item in that sequence. If they tie it takes the first one first (so almendra, then coconut etc)
tool_items$sequenceID <- rownames(tool_items)
tool_items2 <- tool_items[,c("sequenceID", "seq_toolitem")]

agoutigross <- left_join(agoutigross, tool_items2, "sequenceID")

# normal items
agoutigross_foraging <- agoutigross[str_detect(agoutigross$behavior, "\\bF:"),]
normal_items <- as.data.frame(as.matrix(ftable(agoutigross_foraging$sequenceID, agoutigross_foraging$normal_item)))
# set unknown foraging to 1 for all rows, so that when there are two things being processed (e.g. one unknown, one almendra) you keep the almendra
# this is crude solution but works
normal_items$unknown <- 1
normal_items$seq_normalitem <- colnames(normal_items)[max.col(normal_items,ties.method="first")]  
# now gets the most frequently processed item in that sequence. If they tie it takes the first one first (so almendra, then coconut etc)
normal_items$sequenceID <- rownames(normal_items)
normal_items2 <- normal_items[,c("sequenceID", "seq_normalitem")]

agoutigross <- left_join(agoutigross, normal_items2, "sequenceID")

# see whether they both forage with and without tools and create variable with most information
# so what's being foraged irrespective of with or without tools
agoutigross$bothforage <- ifelse(is.na(agoutigross$seq_toolitem) == FALSE & is.na(agoutigross$seq_normalitem) == FALSE, 1, 0)

# how many cases in which they forage with and without tools
unique(agoutigross$behavior[which(agoutigross$bothforage == 1)])

# make general item one, if they both use tools and forage normally, tools take precedent
agoutigross$seq_item <- ifelse(agoutigross$bothforage == 0 & is.na(agoutigross$seq_toolitem) == FALSE, agoutigross$seq_toolitem,
                           ifelse(agoutigross$bothforage == 0 & is.na(agoutigross$seq_normalitem) == FALSE, agoutigross$seq_normalitem,
                                  ifelse(agoutigross$bothforage == 1 & agoutigross$seq_toolitem == "unknown", agoutigross$seq_normalitem,
                                         ifelse(agoutigross$bothforage == 1 & agoutigross$seq_toolitem != "unknown", agoutigross$seq_toolitem, NA))))

# simplify to categories that overlap
agoutigross$seq_item <- ifelse(str_detect(agoutigross$seq_item, "crab"), "crab", agoutigross$seq_item)

### APPROACH 2: get number of capuchins foraging certain resource per sequence for each itemtype
# for how many capuchins are normal item and tool item the same?
agoutigross[which(agoutigross$normal_item == agoutigross$tool_item),]
# wouldn't want to inflate this
agoutigross$foraging_item1 <- agoutigross$normal_item
agoutigross$foraging_item2 <- agoutigross$tool_item
agoutigross$foraging_item2[which(agoutigross$normal_item == agoutigross$tool_item)] <- NA

agoutigross[which(is.na(agoutigross$foraging_item1) == FALSE & is.na(agoutigross$foraging_item2) == FALSE),]

# combine both tool_items and normal_items (if they are not the same)
# sum normal and tool items. This means that if a capuchin eats two different things they are counted twice
# so if 4 out of 5 capuchins eat insects in a sequence, still 3 could be eating almendras (cause they eat several things at once)
# this seems fair. Can also do this otherwise. 

agoutigross_allforaging <- agoutigross[str_detect(agoutigross$behavior, "F:"),]

d1 <- as.data.frame(as.matrix(ftable(agoutigross_allforaging$sequenceID, agoutigross_allforaging$foraging_item1)))
d2 <- as.data.frame(as.matrix(ftable(agoutigross_allforaging$sequenceID, agoutigross_allforaging$foraging_item2)))
d1$sequenceID <- rownames(d1)
d2$sequenceID <- rownames(d2)

all_items <- full_join(d1, d2, by = "sequenceID")

# only make the four categories now
all_items$nr_almendra <- all_items$almendra.x + all_items$almendra.y
all_items$nr_coconut <- all_items$coconut.x + all_items$coconut.y
all_items$nr_fruit <- all_items$fruit.x + all_items$fruit.y + all_items$palm
all_items$nr_invertebrate <- all_items$hecrab + all_items$hwcrab + all_items$insect.x + 
  all_items$insect.y + all_items$snail.x + all_items$snail.y
all_items$nr_other <- all_items$other.x + all_items$other.y
all_items$nr_unknown <- all_items$unknown.x + all_items$unknown.y

all_items2 <- all_items[,c("sequenceID", "nr_almendra", "nr_coconut", "nr_fruit", "nr_invertebrate",
                          "nr_other", "nr_unknown")]

agoutigross <- left_join(agoutigross, all_items2, "sequenceID")

##### Tolerated scrounging and foraging on anvil debris ######
cap_scrounge <- agoutigross[str_detect(agoutigross$behavior, "scroung") & agoutigross$n_tooluse >0,]
cap_scrounge_as <- as.data.frame(as.matrix(ftable(cap_scrounge$sequenceID, cap_scrounge$agesexF)))
colnames(cap_scrounge_as) <- c("sc_nAF", "sc_nAM", "sc_nAU", "sc_nJF", "sc_nJM", "sc_nJU", "sc_nSF", "sc_nSM", "sc_nSU", 
                               "sc_nUF", "sc_nUM", "sc_nUU")
cap_scrounge_as$sequenceID <- rownames(cap_scrounge_as)

agoutigross <- left_join(agoutigross, cap_scrounge_as, "sequenceID")

# anvil debris
cap_debris <- agoutigross[str_detect(agoutigross$behavior, "debris") & agoutigross$capuchin == 1,]
cap_debris_as <- as.data.frame(as.matrix(ftable(cap_debris$sequenceID, cap_debris$agesexF)))
colnames(cap_debris_as) <- c("ad_nAF", "ad_nAM", "ad_nAU", "ad_nJF", "ad_nJM", "ad_nJU", "ad_nSF", "ad_nSM", "ad_nSU", 
                               "ad_nUF", "ad_nUM", "ad_nUU")
cap_debris_as$sequenceID <- rownames(cap_debris_as)

agoutigross <- left_join(agoutigross, cap_debris_as, "sequenceID")

#### Cleaning and to sequence-level ##### 
# set all NAs to 0 except for categorical variables
agoutigross <- agoutigross %>%
  mutate_at(vars("n", "nAF", "nAM", "nAU", "nJF", "nJM", "nJU", "nSF", "nSM", "nSU", "nUF", "nUM", "nUU", "nAdult",
                 "nSubadult", "nJuvenile", "nFemales", "nMales", "tu_nAF", "tu_nAM", "tu_nAU", "tu_nJF", "tu_nJM", 
                 "tu_nJU", "tu_nSF", "tu_nSM", "tu_nSU", "tu_nUF", "tu_nUM", "tu_nUU", "tu_nAdult", "tu_nSubadult", 
                 "tu_nJuvenile", "nr_almendra", "nr_coconut", "nr_fruit", "nr_invertebrate", "nr_other", "nr_unknown",
                 "sc_nAF", "sc_nAM", "sc_nAU", "sc_nJF", "sc_nJM", "sc_nJU", "sc_nSF", "sc_nSM", "sc_nSU", 
                 "sc_nUF", "sc_nUM", "sc_nUU", "ad_nAF", "ad_nAM", "ad_nAU", "ad_nJF", "ad_nJM", "ad_nJU", "ad_nSF", 
                 "ad_nSM", "ad_nSU", "ad_nUF", "ad_nUM", "ad_nUU", "nAF_infant", "nAM_infant", "nAU_infant", "nJF_infant", "nJM_infant", "nJU_infant", "nSF_infant", "nSM_infant", "nSU_infant", 
                 "nUF_infant", "nUM_infant", "nUU_infant", "n_neckinfant", "n_inspect_adult", "n_inspect_nonadult", "n_inspect_unknown"), ~replace_na(.,0))

# make total inspect variable
agoutigross$n_inspect <- agoutigross$n_inspect_adult + agoutigross$n_inspect_nonadult + agoutigross$n_inspect_unknown

### at this point can FILTER OUR THE BLANKS
# remove all timelapse triggers
# flag the uncoded sequences in here that are not camera setup
agoutigross$uncoded <- ifelse(agoutigross$observationType == "unclassified" & agoutigross$cameraSetup == "False" & agoutigross$classificationTimestamp == "" & agoutigross$captureMethod == "activityDetection", 1, 0)

# exclude blanks, unclassified and unknown
ftable(agoutigross$observationType)
# exclude all timelapse triggers, all blanks. Keep the unclassified (uncoded) ones in cause those allow you to filter out which deployments are coded or not
agoutigross2 <- agoutigross[which(agoutigross$observationType != "blank" & agoutigross$captureMethod == "activityDetection"),]

# can still clean up by removing unnecessary columns
# keep checking if these are the right ones to remove
agouticlean <- agoutigross2[, !names(agoutigross2) %in% c("eventStart", "eventID", "eventEnd",  "mediaID", "observationLevel", "cameraSetupType", "classificationProbability", "coordinateUncertainty",
                                                          "start", "end", "cameraID", "cameraModel", "baitUse", "featureType", "individualPositionRadius", "individualPositionAngle", "individualSpeed",
                                                          "bboxX", "bboxY", "bboxWidth", "bboxHeight", "cameraID", "cameraModel", "cameraDelay", "cameraHeight", "cameraDepth", "cameraTilt", "cameraHeading",
                                                          "detectionDistance", "timestampIssues", "featureType", "deploymentGroups", "deploymentTags", "baitUse",
                                                          "habitat", "time", "comment_item", "foraging_item1", "X", "foraging_item2", "countNew", "cameraHeading")]

# drop duplicated sequences, necessary form for tidal analyses and data exploration
# lose individual info etc so to analyze like that need to use agoutigross
agoutisequence <- agouticlean[!duplicated(agouticlean$sequenceID),]
agoutisequence$tooluse <- agoutisequence$n_tooluse > 0
length(unique(agoutisequence$sequenceID)) 

#### Add tide times ####
# Run tide_cleaning.R script first!!

# for each sequence get time to nearest low tide (need to match day and get low tide times then)
# two options, either get absolute difference (so always positive)
# or before/after difference
# Absolute
agoutisequence$tidedifabs <- NA
which(is.na(agoutisequence$time))

for (i in 1:nrow(agoutisequence)) {
  agoutisequence$tidedifabs[i] <- min(abs(difftime(agoutisequence$seq_start[i], TidesLow$TIDE_TIME, units = "hours")))
}

# Both positive and negative

agoutisequence$tidedif <- NA

for (i in 1:nrow(agoutisequence)) {
  agoutisequence$tidedif[i] <- Closest((as.vector(difftime(agoutisequence$seq_start[i], TidesLow$TIDE_TIME,   units = "hours"))), 0)
}

# get an error but it does seem to work... 

## time to high tide instead of low
agoutisequence$tidedif2 <- NA
for (i in 1:nrow(agoutisequence)) {
  agoutisequence$tidedif2[i] <- Closest((as.vector(difftime(agoutisequence$seq_start[i], TidesHigh$TIDE_TIME,   units = "hours"))), 0)
}

#### Further prep and filter uncoded out #####
# make temperature numerical
agoutisequence$temperature <- as.numeric(agoutisequence$temperature)
# add seqday variable (RDate format)
agoutisequence$seqday <- as.Date(format(agoutisequence$seq_start, "%Y-%m-%d"))

## EXPOSURE
# how long the camera was running on a day
## NOTE: in later deployments (R10 on) to conserve batteries we have cameras on a schedule to be off at night, but before run for 24 hours
# extract the day of the sequence and also the day of the deployment start and end
# if day of sequence is in day of deployment start or end, then it's time for that. otherwise it's 24 hours
agoutisequence$dep_startday <- format(agoutisequence$dep_start, "%Y-%m-%d")
agoutisequence$dep_endday <- format(agoutisequence$dep_end, "%Y-%m-%d")
agoutisequence$seq_startday <- format(agoutisequence$seq_start, "%Y-%m-%d")
agoutisequence$dep_starttime <- hour(agoutisequence$dep_start)
agoutisequence$dep_endtime <- hour(agoutisequence$dep_end)

# maybe too coarse, only on level of hours now 
agoutisequence$exposure <- ifelse(agoutisequence$seq_startday == agoutisequence$dep_startday, 
                                  24-agoutisequence$dep_starttime, 
                                  ifelse(agoutisequence$seq_startday == agoutisequence$dep_endday,
                                         24 - agoutisequence$dep_endtime, 24))

# need to filter out deployments that are not fully coded
# have a lot of unclassified sequences in --> these seem to largely be the 00:00:00 automated timecapture moments, that weren't coded
# so have two variables for this, whether a sequence is a timelapse sequence (1 yes, 0 no) and whether it is uncoded (1 yes, 0 no)
# want to exclude deployments that have uncoded sequences (that are not timelapse)
# I already took the timelapse ones out, so all unclassified left are "real" uncoded sequences
# can either very strictly subset on only 100% coded deployments or less strictly on all that have less than 5 uncoded sequences or something
cd <- as.data.frame((ftable(agoutisequence$uniqueloctag, agoutisequence$uncoded)))
codeddeployments_total <- as.character(cd$Var1[cd$Var2 == 1 & cd$Freq < 5]) # deployments that miss less than 5 sequences

# subset only fully coded deployments 
agoutisequence_c <- agoutisequence[(agoutisequence$uniqueloctag %in% codeddeployments_total),]
agoutisequence_c <- droplevels.data.frame(agoutisequence_c)

agoutisequence_c$hour <- hour(agoutisequence_c$seq_start)
agoutisequence_c$toolusers <- factor(agoutisequence_c$tool_site, levels = c(0,1), labels = c("Non-tool-users", "Tool-users"))
agoutisequence_c$locationfactor <- as.factor(agoutisequence_c$locationName)
agoutisequence_c$depdays <- as.numeric(difftime(agoutisequence_c$seq_startday, agoutisequence_c$dep_startday, units = "days"))
agoutisequence_c$depnr <- as.numeric(unlist(regmatches(agoutisequence_c$tags, gregexpr("[[:digit:]]+", agoutisequence_c$tags))))

## also want to flag human presence, which is when cameras were deployed or picked up. 
# here try to be specific, e.g. if one camera is placed in the tool-user range on Jicaron on a specific date, all cameras running there at that time might be affected. 
# so use variable that specifies which island and which "side" (north/west) etc. 
# and then identify days that a camera was picked up or removed
agoutisequence_c$islandside <- paste(agoutisequence_c$island, agoutisequence_c$side, sep = "-")
picksetupdays <- agoutisequence_c[agoutisequence_c$cameraSetup ==  "True", c("seqday", "islandside")]
picksetupdays <- picksetupdays[!duplicated(picksetupdays),]

# create flag for which ones need to be removed
islandsides <- unique(agoutisequence_c$islandside)
agoutisequence_c$picksetup <- NA
for (i in 1:length(islandsides)) {
  agoutisequence_c$picksetup[which(agoutisequence_c$islandside == islandsides[i])] <- ifelse(agoutisequence_c$seqday[which(agoutisequence_c$islandside == islandsides[i])] %in% picksetupdays$seqday[which(picksetupdays$islandside == islandsides[i])], 1, 0)
}

# now removing these is as simple as filtering out all those where picksetup == 1

##### generating 0's ####
## NOTE: this does not yet take into consideration that some cameras may not be active at night (this is only later deployments).

# for some purposes we want to know when the camera was active but not triggered, these 0's are not included yet
# currently if we look at 0's it's only triggers by other species
# so need to make this explicit to be able to distinguish between absence of data because camera was not triggered vs camera was not deployed

# do this on the big full agoutisequence_c dataset so that you can easily subset from it

## need to know start and end date of each unique deployment
locations2 <- data.frame(uniqueloctag = unique(agoutisequence_c$uniqueloctag)) 
locations2 <- left_join(locations2, agoutisequence_c[,c("uniqueloctag", "dep_start", "dep_end")], by = "uniqueloctag")
locations2 <- locations2[!duplicated(locations2$uniqueloctag),]
# take time off and keep just date variable
locations2$dep_startday <- as.Date(format(locations2$dep_start, "%Y-%m-%d"))
locations2$dep_endday <- as.Date(format(locations2$dep_end, "%Y-%m-%d"))

###  generate all the days that should be present within each deployment
# first create dataframe for first one
depldays3 <<- data.frame(uniqueloctag = locations2$uniqueloctag[1], seqday = seq(locations2$dep_startday[1], locations2$dep_endday[1], by = "days"))
# now iterate over all the other ones and append those to the dataframe
for (i in 2:nrow(locations2)) {
  depldays4 = data.frame(uniqueloctag = locations2$uniqueloctag[i], seqday = seq(locations2$dep_startday[i], locations2$dep_endday[i], by = "days"))
  depldays3 <<- rbind(depldays3, depldays4)
} 

# need to expand this dataframe to have an hour a day for each deployment.
depldayhour <- data.frame(uniqueloctag = rep(depldays3$uniqueloctag, 24), seqday = rep(depldays3$seqday, 24))
depldayhour <- depldayhour[order(depldayhour$uniqueloctag, depldayhour$seqday),]
depldayhour$hour <- rep(1:24, (nrow(depldayhour)/24))

agoutiselect2 <- left_join(depldayhour[depldayhour$uniqueloctag %in% agoutisequence_c$uniqueloctag,], agoutisequence_c, by = c("uniqueloctag", "seqday", "hour"), relationship = "many-to-many")
agoutiselect2$noanimal <- ifelse(is.na(agoutiselect2$sequenceID), 1, 0)

## fill in NAs like above, using a metadata file
metadata2 <- agoutisequence_c[!duplicated(agoutisequence_c$uniqueloctag), c("uniqueloctag", "deploymentID", "locationName", "tags", "dep_start", "dep_end", "dep_length_hours",
                                                                            "island", "tool_anvil", "tool_site", "toolusers", "streambed", "mediatype", "depnr", "depdays")]

for (i in 1:nrow(metadata2)) {
  agoutiselect2$locationName[agoutiselect2$uniqueloctag == metadata2$uniqueloctag[i]] <- metadata2$locationName[metadata2$uniqueloctag == metadata2$uniqueloctag[i]]
  agoutiselect2$tags[agoutiselect2$uniqueloctag == metadata2$uniqueloctag[i]] <- metadata2$tags[metadata2$uniqueloctag == metadata2$uniqueloctag[i]]
  agoutiselect2$dep_start[agoutiselect2$uniqueloctag == metadata2$uniqueloctag[i]] <- metadata2$dep_start[metadata2$uniqueloctag == metadata2$uniqueloctag[i]] 
  agoutiselect2$dep_end[agoutiselect2$uniqueloctag == metadata2$uniqueloctag[i]] <- metadata2$dep_end[metadata2$uniqueloctag == metadata2$uniqueloctag[i]] 
  agoutiselect2$dep_length_hours[agoutiselect2$uniqueloctag == metadata2$uniqueloctag[i]] <- metadata2$dep_length_hours[metadata2$uniqueloctag == metadata2$uniqueloctag[i]] 
  agoutiselect2$island[agoutiselect2$uniqueloctag == metadata2$uniqueloctag[i]] <- metadata2$island[metadata2$uniqueloctag == metadata2$uniqueloctag[i]] 
  agoutiselect2$tool_anvil[agoutiselect2$uniqueloctag == metadata2$uniqueloctag[i]] <- metadata2$tool_anvil[metadata2$uniqueloctag == metadata2$uniqueloctag[i]]
  agoutiselect2$deploymentID[agoutiselect2$uniqueloctag == metadata2$uniqueloctag[i]] <- metadata2$deploymentID[metadata2$uniqueloctag == metadata2$uniqueloctag[i]]
  agoutiselect2$tool_site[agoutiselect2$uniqueloctag == metadata2$uniqueloctag[i]] <- metadata2$tool_site[metadata2$uniqueloctag == metadata2$uniqueloctag[i]]
  agoutiselect2$streambed[agoutiselect2$uniqueloctag == metadata2$uniqueloctag[i]] <- metadata2$streambed[metadata2$uniqueloctag == metadata2$uniqueloctag[i]]
  agoutiselect2$mediatype[agoutiselect2$uniqueloctag == metadata2$uniqueloctag[i]] <- metadata2$mediatype[metadata2$uniqueloctag == metadata2$uniqueloctag[i]]
  agoutiselect2$toolusers[agoutiselect2$uniqueloctag == metadata2$uniqueloctag[i]] <- metadata2$toolusers[metadata2$uniqueloctag == metadata2$uniqueloctag[i]]
  agoutiselect2$depnr[agoutiselect2$uniqueloctag == metadata2$uniqueloctag[i]] <- metadata2$depnr[metadata2$uniqueloctag == metadata2$uniqueloctag[i]]
  agoutiselect2$depdays[agoutiselect2$uniqueloctag == metadata2$uniqueloctag[i]] <- metadata2$depdays[metadata2$uniqueloctag == metadata2$uniqueloctag[i]]
  
  
} 

agoutiselect2$n[agoutiselect2$noanimal == 1] <- 0
agoutiselect2$capuchin[agoutiselect2$noanimal == 1] <- 0
# not sure if we should change seq_length to 0
agoutiselect2$seq_length[agoutiselect2$noanimal == 1] <- 0
agoutiselect2[,c( "n", "nAF", "nAM", "nAU","nJF", "nJM","nJU","nSF","nSM","nSU","nUF","nUM","nUU","nAF_infant","nAM_infant","nAU_infant","nJF_infant","nJM_infant","nJU_infant",
                  "nSF_infant","nSM_infant","nSU_infant", "nUF_infant","nUM_infant","nUU_infant", "n_neckinfant","nAdult","nSubadult","nJuvenile","nFemales",  "nMales","n_inspect",
                  "n_tooluse","tu_nAF","tu_nAM","tu_nAU", "tu_nJF","tu_nJM","tu_nJU","tu_nSF","tu_nSM","tu_nSU","tu_nUF","tu_nUM","tu_nUU" , "tu_nAdult","tu_nSubadult", "tu_nJuvenile",
                  "nr_almendra","nr_coconut","nr_fruit","nr_invertebrate","nr_other" , "nr_unknown", "sc_nAF","sc_nAM","sc_nAU","sc_nJF","sc_nJM","sc_nJU","sc_nSF","sc_nSM","sc_nSU",
                  "sc_nUF","sc_nUM","sc_nUU","ad_nAF","ad_nAM","ad_nAU","ad_nJF","ad_nJM","ad_nJU", "ad_nSF","ad_nSM","ad_nSU","ad_nUF", "ad_nUM")][is.na(agoutiselect2[,c(
                    "n", "nAF", "nAM", "nAU","nJF", "nJM","nJU","nSF","nSM","nSU","nUF","nUM","nUU","nAF_infant","nAM_infant","nAU_infant","nJF_infant","nJM_infant","nJU_infant",
                    "nSF_infant","nSM_infant","nSU_infant", "nUF_infant","nUM_infant","nUU_infant", "n_neckinfant","nAdult","nSubadult","nJuvenile","nFemales",  "nMales","n_inspect",
                    "n_tooluse","tu_nAF","tu_nAM","tu_nAU", "tu_nJF","tu_nJM","tu_nJU","tu_nSF","tu_nSM","tu_nSU","tu_nUF","tu_nUM","tu_nUU" , "tu_nAdult","tu_nSubadult", "tu_nJuvenile",
                    "nr_almendra","nr_coconut","nr_fruit","nr_invertebrate","nr_other" , "nr_unknown", "sc_nAF","sc_nAM","sc_nAU","sc_nJF","sc_nJM","sc_nJU","sc_nSF","sc_nSM","sc_nSU",
                    "sc_nUF","sc_nUM","sc_nUU","ad_nAF","ad_nAM","ad_nAU","ad_nJF","ad_nJM","ad_nJU", "ad_nSF","ad_nSM","ad_nSU","ad_nUF", "ad_nUM")])] <- 0

# fill in important variables
agoutiselect2$month <- month(agoutiselect2$seqday)

# make wet season may-nov and dry season dec-april
agoutiselect2$season <- ifelse(agoutiselect2$month == 12 | agoutiselect2$month == 1 | agoutiselect2$month == 2 | agoutiselect2$month == 3 | 
                               agoutiselect2$month == 4, "Dry", "Wet") 

agoutiselect2$locationfactor <- as.factor(agoutiselect2$locationName)
agoutiselect2$tooluse[agoutiselect2$noanimal == 1] <- 0

# Save dataframes as RDS so you don't always have to re-run this script but can just load these in
saveRDS(agoutiselect2, "agouti_output/agoutiselect2.rds")
saveRDS(agouticlean, "agouti_output/agouticlean.rds")
saveRDS(agoutisequence, "agouti_output/agoutisequence.rds")
saveRDS(agoutisequence_c, "agouti_output/agoutisequence_c.rds")

agoutiselect2 <- readRDS("agouti_output/agoutiselect2.rds")
agouticlean <- readRDS("agouti_output/agouticlean.rds")
agoutisequence <- readRDS("agouti_output/agoutisequence.rds")
agoutisequence_c <- readRDS("agouti_output/agoutisequence_c.rds")
