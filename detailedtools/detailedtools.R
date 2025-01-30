## Detailed tool use analyses -- Cleaning Script
## MPI-AB; Z Goldsborough

## Script how to clean BORIS output data

## packages needed
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
require(DescTools)
library(reshape2)
library(data.table)

## NOTE: add code about catching when anvil type switches (e.g. if it's wooden anvil at comment seq_start, or when they switch back to stone/wood)
## NOTE: check how to extract end location of second hammer when they switch (code hammerend loc again?)

### Loading dataset ####
# load csv files with aggregated BORIS output (from BORIScoding/ExportedData Google Drive)
# Zoë's csv
dettools1 <- read.csv("detailedtools/ZGdetailedtoolscoding.csv")
# Meredith's csv
dettools2 <- read.csv("detailedtools/EXP-ANV-01-R11_MC.csv")
dettools2$Coder.ID <- "MC"
# Leonie's csv
dettools3 <- read.csv("detailedtools/CEBUS-02-R11_R12_2022_EXP-ANV-01_R12_LRdetailedtoolscoding.csv")

# bind all three datasets together (after making sure they have the same number and order of columns)
# TEMPORARY until we are all on same version
new_cols <- names(dettools1)[!names(dettools1) %in% names(dettools2)]
old_cols <- names(dettools2)[!names(dettools2) %in% names(dettools1)]

dettools1 <- dettools1[, !names(dettools1) %in% new_cols] 
dettools2 <- dettools2[, !names(dettools2) %in% old_cols]
dettools3 <- dettools3[, !names(dettools3) %in% old_cols]
## END OF TEMPORARY

# for now not 2b in as there are still errors/is incomplete
dettools <- rbind(dettools1, dettools2, dettools3)
# sort so that observations from the same video are clustered together and it's chronological
dettools <- dettools[order(dettools$Observation.id),]

# remove unnecessary columns and rename the ones we keep
dettools_r <- data.frame("videoID" = dettools$Observation.id, "codingdate" = dettools$Observation.date,
                         "medianame" = dettools$Media.file.name, "videolength" = dettools$Media.duration..s., "coder" = 
                           dettools$Coder.ID, "subjectID" = dettools$Subject, "behavior" = dettools$Behavior,
                         "modifier1" = dettools$Modifier..1,  "modifier2" = dettools$Modifier..2,  "modifier3" = dettools$Modifier..3,  "modifier4" = dettools$Modifier..4, 
                         "starttime" = dettools$Start..s., "comment" = dettools$Comment.start)

# take out Zoë's coding of "female" tool use for other project
# note: if we have any other test sequences, we can filter them out here
# for now also filtering out Leonie's sequence where BAL was processing two almendras at once alternating between them
flags <- dettools_r$videoID[which(str_detect(dettools_r$videoID, "femaletooluse1|test|double") == TRUE)]
dettools_r <- dettools_r[!dettools_r$videoID %in% flags,]

### Create unique sequence ID #### 
# Sequence ID that is same for sequences continuing across multiple videos
# and different for 2 or 3 sequences in the same video
### NOTE: it is crucial for this to work that seq_start is the first thing and seq_end the last thing in each sequence!!!! 

# change seq_end to seq_cont if it has the continue modifier
dettools_r$behavior[which(str_detect(dettools_r$modifier1, "cont") == TRUE)] <- "seqcont"

# create ascending number for each sequence
curseq <- 1
cache <- 0
dettools_r$seqnumber <- NA

for (i in 1:nrow(dettools_r)) {
  dettools_r$seqnumber[i] <- curseq
  
  if(dettools_r$behavior[i] == "seqend") {
    cache <- curseq
    dettools_r$seqnumber[i] <- curseq
    curseq <- NA
  }
  if(dettools_r$behavior[i] == "seqstart") {
  curseq <- cache + 1
  dettools_r$seqnumber[i] <- curseq
  cache <- NA
  }

}

# check for NAs
dettools_r[is.na(dettools_r$seqnumber) == TRUE,]
# check for focalsubjects not coded or unknown
unknownIDs <- unique(dettools_r$videoID[which(dettools_r$subjectID == "No focal subject" | dettools_r$subjectID == "unknown")])
mistakenIDs <- dettools_r[dettools_r$videoID %in% unknownIDs,]

# combine with location and date to get unique seq_ID
dettools_r$location <- ifelse(str_detect(dettools_r$videoID, "EXP-ANV") == TRUE, "EXP-ANV-01", "CEBUS-02")
dettools_r$mediadate <- sapply(str_split(dettools_r$videoID, "__"), '[', 2)
dettools_r$sequenceID <- paste(dettools_r$location, dettools_r$mediadate, dettools_r$seqnumber, sep = "_" )

head(dettools_r$sequenceID)

# if there was a comment that the ending was missed or unknown, then change outcome from "none"  to "unknown"
dettools_r$modifier1[which(dettools_r$behavior == "seqend" & dettools_r$modifier1 == "None" &
                             str_detect(dettools_r$comment, "unknown|missed|missing|not done|not continued|not on video|not sure|unkown|ends before|not on camera") == TRUE)] <- "Unknown"

### Extract modifiers per sequence ####
# make dataframe with sequence_level information and populate it, then left_join at the end
seqdat <- data.frame(sequenceID = unique(dettools_r$sequenceID))
# here every row is a sequence

## what item is being consumed
seqdat$item[which(seqdat$sequenceID == dettools_r[dettools_r$behavior == "seqstart",]$sequenceID)] <-  dettools_r[dettools_r$behavior == "seqstart",]$modifier1[which(seqdat$sequenceID == dettools_r[dettools_r$behavior == "seqstart",]$sequenceID)]
ftable(seqdat$item)
# see overwhelming majority of consumed items are almendras

## hammerstone information
# note: anvil location can be in modifier 3 or modifier 4
hammers <- dettools_r[dettools_r$behavior == "hammerstone",]
hammers$h_startloc <- hammers$modifier1
hammers$h_endloc <- ifelse(hammers$modifier4 == "", hammers$modifier3, hammers$modifier4)
hammers$h_endloc <- ifelse(hammers$h_endloc == "None",
                           lead(hammers$h_endloc, order_by = hammers$sequenceID),
                           hammers$h_endloc)
hammers <- hammers[which(hammers$h_startloc != "None"),]
# check that the next line works correctly to identify unmarked and unknown hammerstones
hammers$hammerID <- ifelse(hammers$modifier2 != "None", hammers$modifier2, hammers$comment)
hammers <- hammers[,c("sequenceID", "h_startloc", "h_endloc", "hammerID")]
# sometimes there is additional whitespace or a comment attached to the hammerstone ID, pull those out
# first exclude the ones where hammerstone ID is as it should be (or has check attached to it)
hammers$hammerID[! hammers$hammerID %in% c("FRE", "unmarked", "BAM", "PEB", "unknown", "WIL", "DWA", "BCH", "LCH", 
                                           "DWA_A", "DWA_B", "DPL", "DPL_A", "check ID", "DPL_A CHECK", "")] <-
  str_trim(str_extract(hammers$hammerID[! hammers$hammerID %in% c("FRE", "unmarked", "BAM", "PEB", "unknown", "WIL", "DWA", "BCH", "LCH", 
                                                                  "DWA_A", "DWA_B", "DPL", "DPL_A", "check ID", "DPL_A CHECK", "")],
                       "([:upper:]|[:space:]){2,}"))

# generate lists of blanks that we missed so we can correct them
blank <- hammers$sequenceID[which(hammers$hammerID == "")]
blank_videonames_ZG <- unique(dettools_r$videoID[which(dettools_r$sequenceID %in% blank & dettools_r$coder == "ZG")])
blank_videonames_MKWC <- unique(dettools_r$videoID[which(dettools_r$sequenceID %in% blank & dettools_r$coder == "MKWC")])
blank_videonames_LR <- unique(dettools_r$videoID[which(dettools_r$sequenceID %in% blank & dettools_r$coder == "LR")])

seqdat <- left_join(seqdat, hammers, "sequenceID")
ftable(seqdat$hammerID)

## sequence end information
seqendings <- dettools_r[dettools_r$behavior == "seqend" | dettools_r$behavior == "seqcont",]
# outcome
seqendings$outcome <- seqendings$modifier1
ftable(seqendings$outcome)
# displacement
seqendings$displacement <- seqendings$modifier3
ftable(seqendings$displacement)
# social attention
seqendings$socatt <- seqendings$modifier4
ftable(seqendings$socatt)
# scrounging
seqendings$scrounging <- seqendings$modifier2
ftable(seqendings$scrounging)

# filter out seqcont ones after making sure their information is included
seqendings <- seqendings[!seqendings$outcome == "seqcont", c("sequenceID", "outcome", "displacement", "socatt", "scrounging")]

seqdat <- left_join(seqdat, seqendings, "sequenceID")

# I think that is all the sequence specific info, attach this to the main dataframe
dettools_r2 <- left_join(dettools_r, seqdat, "sequenceID")

### Extract modifiers on behavior level  ####
# (so not aggregated to sequence)
# specify pound type
dettools_r2$poundtype <- ifelse(dettools_r2$behavior == "pound", dettools_r2$modifier1, NA)
# specify one-footed yes/no
dettools_r2$onefoot <- ifelse(dettools_r2$behavior == "pound", str_detect(dettools_r2$modifier2, "1foot"), NA)
# specify overhead yes/no
dettools_r2$overhead <- ifelse(dettools_r2$behavior == "pound", str_detect(dettools_r2$modifier3, "overhead"), NA)
# specify one-handed yes/no
dettools_r2$onehand <- ifelse(dettools_r2$behavior == "pound", str_detect(dettools_r2$modifier2, "1hand"), NA)
# specify tail-support yes/no
dettools_r2$tailsupport <- ifelse(dettools_r2$behavior == "pound", str_detect(dettools_r2$modifier2, "tailsupport"), NA)

# specify mistake type
dettools_r2$mistaketype <- ifelse(dettools_r2$behavior == "misstrike", dettools_r2$modifier1, NA)
# specify repositioning type
dettools_r2$repostype <- ifelse(dettools_r2$behavior == "reposit", dettools_r2$modifier1, NA)

# will have to do when there's a hammerswitch that hammerID changes until seq_end (use for loop for it like above)
currenthammerID <- dettools_r2$hammerID[1]

for (i in 1:nrow(dettools_r2)) {
  dettools_r2$hammerID2[i] <- currenthammerID
   if(dettools_r2$behavior[i] == "hammerswitch") {
    currenthammerID <- ifelse(dettools_r2$modifier2[i] != "None", dettools_r2$modifier2[i], dettools_r2$comment[i])
    dettools_r2$hammerID2[i] <- currenthammerID
  }

  if(dettools_r2$behavior[i] == "seqstart") {
    currenthammerID <- dettools_r2$hammerID[i]
    dettools_r2$hammerID2[i] <- currenthammerID
  }
}

# extract location of hammerstone they switched to 
dettools_r2$h_switchloc <- ifelse(dettools_r2$behavior == "hammerswitch", dettools_r2$modifier1, NA)

# add information on hammerstone switch to sequence-level dataframe
# there are sometimes many hammerstone switches in one sequence, so it is not straightforward to add all the IDs to the sequence level dataframe
# for now add the number of switches per sequence as a sequence-level variable
switchsequences <- dettools_r2$sequenceID[which(dettools_r2$behavior == "hammerswitch")]
switches <- aggregate(dettools_r2$behavior[dettools_r2$sequenceID %in% switchsequences & dettools_r2$behavior == "hammerswitch"], by = list(sequenceID = dettools_r2$sequenceID[dettools_r2$sequenceID %in% switchsequences & dettools_r2$behavior == "hammerswitch"]), FUN = length)
seqdat$hammerswitches <- 0
seqdat$hammerswitches[seqdat$sequenceID %in% switchsequences] <- switches$x
ftable(seqdat$hammerswitches)

## incorporate what type of anvil it is
# at EXP-ANV is stone, at CEBUS-02 is wood
seqdat$anviltype <- ifelse(str_detect(seqdat$sequenceID, "CEBUS-02") == TRUE, "wood", "stone")
# if it differed from the anviltype of the main anvil we made a comment at the seqstart
seqdat$anviltype <- ifelse(seqdat$sequenceID %in% dettools_r2$sequenceID[which(dettools_r2$behavior == "seqstart" & str_detect(dettools_r2$comment, "wood") == TRUE)], "wood", ifelse(
  seqdat$sequenceID %in% dettools_r2$sequenceID[which(dettools_r2$behavior == "seqstart" & str_detect(dettools_r2$comment, "stone") == TRUE)], "stone", seqdat$anviltype))
# attach this to main dataframe
dettools_r2 <- left_join(dettools_r2, seqdat[,c("sequenceID", "anviltype")], "sequenceID")

# when there's an anvilswitch then anviltype changes until seq_end (use for loop for it like above)
currentanviltype <- dettools_r2$anviltype[1]

for (i in 1:nrow(dettools_r2)) {
  dettools_r2$anviltype2[i] <- currentanviltype
  if(dettools_r2$behavior[i] == "anvilswitch") {
    currentanviltype <- dettools_r2$modifier1[i]
    dettools_r2$anviltype2[i] <- currentanviltype
  }
  
  if(dettools_r2$behavior[i] == "seqstart") {
    currentanviltype <- dettools_r2$anviltype[i]
    dettools_r2$anviltype2[i] <- currentanviltype
  }
}

# add number of anvilswitches to seqdat
aswitchsequences <- dettools_r2$sequenceID[which(dettools_r2$behavior == "anvilswitch")]
aswitches <- aggregate(dettools_r2$behavior[dettools_r2$sequenceID %in% aswitchsequences & dettools_r2$behavior == "anvilswitch"], by = list(sequenceID = dettools_r2$sequenceID[dettools_r2$sequenceID %in% aswitchsequences & dettools_r2$behavior == "anvilswitch"]), FUN = length)
seqdat$anvilswitches <- 0
seqdat$anvilswitches[seqdat$sequenceID %in% aswitchsequences] <- aswitches$x
ftable(seqdat$anvilswitches)

## so in dettools_r2 anviltype2 and hammerID2 are the correct variables

### Get metrics of efficiency aggregated per sequence #### 
# calculate duration 
# need to extract the start and end time of the videos
dettools_r2$videostart <- as.POSIXct(paste(sapply(str_split(dettools_r2$videoID, "__"), '[', 2), sapply(str_split(dettools_r2$videoID, "__"), '[', 3), sep = " "), tz = "America/Panama", format = "%Y-%m-%d %H-%M-%S") 
dettools_r2$videoend <- dettools_r2$videostart + 60

# make flag if it's a split sequence or not
# then calculate duration and  start/end time differently depending on whether it is a split sequence or not. 
for (i in 1:nrow(seqdat)) {
  seqdat$split[i] <- ifelse("seqcont" %in% dettools_r2$behavior[which(dettools_r2$sequenceID== seqdat$sequenceID[i])], TRUE, FALSE)
  seqdat$seqstart[i] <- dettools_r2$starttime[which(dettools_r2$sequenceID == seqdat$sequenceID[i] & dettools_r2$behavior == "seqstart")]
  
  # if sequence is within one video it is just the time of seqstart and seqend
  if( seqdat$split[i] == FALSE) { 
  seqdat$seqend[i] <- dettools_r2$starttime[which(dettools_r2$sequenceID == seqdat$sequenceID[i] & dettools_r2$behavior == "seqend")]
  seqdat$seqduration[i] <- seqdat$seqend[i] - seqdat$seqstart[i]
  }
  
  # if sequence is split over consecutive videos, we use the real time (the start time of the video + timestamp of seqstart and seqend) to calculate the duration
  # this also includes the time when the camera was not active but retriggering! (usually around 10-30 seconds)
  if (seqdat$split[i] == TRUE) { 
  realstart <- unique(dettools_r2$videostart[which(dettools_r2$sequenceID == seqdat$sequenceID[i] & dettools_r2$behavior == "seqstart")] + seqdat$seqstart[i])  
  realend <- unique(dettools_r2$videostart[which(dettools_r2$sequenceID == seqdat$sequenceID[i] & dettools_r2$behavior == "seqend")] + 
                        dettools_r2$starttime[which(dettools_r2$sequenceID == seqdat$sequenceID[i] & dettools_r2$behavior == "seqend")])
  seqdat$seqduration[i] <- as.numeric(difftime(realend, realstart, unit = "secs"))
  seqdat$seqend[i] <- seqdat$seqstart[i] + seqdat$seqduration[i]
  }
}

# left join again to get this sequence information in 
dettools_r2 <- left_join(dettools_r2, seqdat[, c("sequenceID", "seqduration", "seqstart", "seqend")], "sequenceID")
# sequence duration only makes sense as a metric of efficiency if the outcome was indeed that the item was opened, so will need to filter on that when doing these analyses (via the outcome variable)

# nr of pounds
poundsonly <- dettools_r2[dettools_r2$behavior == "pound",]
nr_pounds <- poundsonly %>%
  dplyr::count(sequenceID)

colnames(nr_pounds) <- c("sequenceID", "n_pounds")

dettools_r2 <- left_join(dettools_r2, nr_pounds, "sequenceID")
# set the NA's to 0
dettools_r2$n_pounds[which(is.na(dettools_r2$n_pounds) == TRUE)] <- 0

# nr of misstrikes (true miss)
missonly <- dettools_r2[dettools_r2$behavior == "misstrike" & str_detect(dettools_r2$mistaketype, "None") == TRUE,]
nr_miss <- missonly %>% 
  dplyr::count(sequenceID)
colnames(nr_miss) <- c("sequenceID", "n_miss")

dettools_r2 <- left_join(dettools_r2, nr_miss, "sequenceID")
dettools_r2$n_miss[which(is.na(dettools_r2$n_miss) == TRUE)] <- 0

# nr of misstrikes (item flies)
fliesonly <- dettools_r2[dettools_r2$behavior == "misstrike" & str_detect(dettools_r2$mistaketype, "itemflies") == TRUE,]
nr_flies <- fliesonly %>% 
  dplyr::count(sequenceID)
colnames(nr_flies) <- c("sequenceID", "n_flies")

dettools_r2 <- left_join(dettools_r2, nr_flies, "sequenceID")
dettools_r2$n_flies[which(is.na(dettools_r2$n_flies) == TRUE)] <- 0

# nr of misstrikes (hammer lost)
hlossonly <- dettools_r2[dettools_r2$behavior == "misstrike" & str_detect(dettools_r2$mistaketype, "hammerlost") == TRUE,]
nr_hloss <- hlossonly %>% 
  dplyr::count(sequenceID)
colnames(nr_hloss) <- c("sequenceID", "n_hloss")

dettools_r2 <- left_join(dettools_r2, nr_hloss, "sequenceID")
dettools_r2$n_hloss[which(is.na(dettools_r2$n_hloss) == TRUE)] <- 0

# make combined metric
dettools_r2$n_misstotal <- dettools_r2$n_miss + dettools_r2$n_flies + dettools_r2$n_hloss

# nr of repositions (hammer)
hreponly <- dettools_r2[dettools_r2$behavior == "reposit" & dettools_r2$repostype == "hammer",]
nr_hrepos <- hreponly %>%
  dplyr::count(sequenceID)
colnames(nr_hrepos) <- c("sequenceID", "n_hamreposit")

dettools_r2 <- left_join(dettools_r2, nr_hrepos, "sequenceID")
dettools_r2$n_hamreposit[which(is.na(dettools_r2$n_hamreposit) == TRUE)] <- 0

# nr of repositions (item)
ireponly <- dettools_r2[dettools_r2$behavior == "reposit" & dettools_r2$repostype == "item",]
nr_irepos <- ireponly %>%
  dplyr::count(sequenceID)
colnames(nr_irepos) <- c("sequenceID", "n_itemreposit")

dettools_r2 <- left_join(dettools_r2, nr_irepos, "sequenceID")
dettools_r2$n_itemreposit[which(is.na(dettools_r2$n_itemreposit) == TRUE)] <- 0

# nr of repositions (peel)
peelonly <- dettools_r2[dettools_r2$behavior == "reposit" & dettools_r2$repostype == "peel",]
nr_peel <- peelonly %>%
  dplyr::count(sequenceID)
colnames(nr_peel) <- c("sequenceID", "n_peel")

dettools_r2 <- left_join(dettools_r2, nr_peel, "sequenceID")
dettools_r2$n_peel[which(is.na(dettools_r2$n_peel) == TRUE)] <- 0

# combined metric (hammer and item repositioning, not peeling)
dettools_r2$n_reposit <- dettools_r2$n_hamreposit + dettools_r2$n_itemreposit

#### add age/sex to IDS ####
# load in file with capuchin age and sexes. KEEP THIS UP TO DATE
capID <- read.csv("detailedtools/capuchinIDs.csv", sep = ",")
dettools_r2 <- left_join(dettools_r2, capID, by = c("subjectID" = "ID"))

# fill in age sex for unidentified individuals
dettools_r2$Age[which(is.na(dettools_r2$Age) == TRUE)] <- ifelse(str_detect(dettools_r2$subjectID[which(is.na(dettools_r2$Age) == TRUE)], "juvenile") == TRUE, "Juvenile", 
                                                                 ifelse(str_detect(dettools_r2$subjectID[which(is.na(dettools_r2$Age) == TRUE)], "adult") == TRUE, "Adult",
                                                                 ifelse(str_detect(dettools_r2$subjectID[which(is.na(dettools_r2$Age) == TRUE)], "subadult") == TRUE, "Subadult", "Unknown"))) 
dettools_r2$Sex[which(is.na(dettools_r2$Sex) == TRUE)] <- ifelse(str_detect(dettools_r2$subjectID[which(is.na(dettools_r2$Sex) == TRUE)], "male") == TRUE, "Male", 
                                                                 ifelse(str_detect(dettools_r2$subjectID[which(is.na(dettools_r2$Sex) == TRUE)], "female") == TRUE, "Female", "Unknown")) 
# have some extra  spaces that snuck in
dettools_r2$Age <- str_trim(dettools_r2$Age)
# make age ordered factor for easy plotting
dettools_r2$age_of <- factor(dettools_r2$Age, ordered = TRUE, levels = c("Juvenile", "Subadult", "Adult"))
dettools_r2$hammerID2 <- str_trim(dettools_r2$hammerID2)
# extract deployment number
dettools_r2$deployment <- ifelse(str_detect(dettools_r2$videoID, "R11") == TRUE, "R11", 
                                 ifelse(str_detect(dettools_r2$videoID, "R12") == TRUE, "R12", "R13"))
# make variable types correct
dettools_r2$mediadate <- as.POSIXct(dettools_r2$mediadate)
# add split to dataframe
dettools_r2 <- left_join(dettools_r2, seqdat[,c("sequenceID", "split")], by = "sequenceID")

### Check for coding errors and fix minor mistakes #### 
# fix occurrences of DWA after fracture (2022-03-11), should be DWA_A
dettools_r2$hammerID[which(dettools_r2$mediadate > "2022-03-11" & dettools_r2$hammerID == "DWA")] <- "DWA_A"
dettools_r2$hammerID2[which(dettools_r2$mediadate > "2022-03-11" & dettools_r2$hammerID2 == "DWA")] <- "DWA_A"

# some histograms to look for mistakes
hist(dettools_r2$n_pounds)
hist(dettools_r2$seqduration)
max(dettools_r2$seqduration)
ftable(dettools_r2$subjectID)
table(dettools_r2$location, dettools_r2$deployment, dettools_r2$coder)
hist(dettools_r2$n_reposit)
ftable(dettools_r2$h_startloc)
ftable(dettools_r2$h_endloc)
unique(dettools_r2$videoID[which(dettools_r2$h_endloc == "None")])
ftable(dettools_r2$Sex)

# sequence level dataframe (only for analyzing things like nr of pounds/duration. things that are fixed per sequence)
detseq <- dettools_r2[!duplicated(dettools_r2$sequenceID),]
detseq <- left_join(detseq, seqdat[,c("sequenceID", "hammerswitches", "anvilswitches")], by = "sequenceID")

# add a pound if the hammerstone was "inhand" or if there is a split
detseq$n_pounds[which(detseq$h_startloc == "inhand" | detseq$split == TRUE)] <- detseq$n_pounds[which(detseq$h_startloc == "inhand" | detseq$split == TRUE)] + 1

### Integrate social attention coding ####
# coded in more detail who is present, who displaces, who scrounges, and who pays social attention (and when)
# can do more with this, but for now will extract per sequence this basic additional info
# can have both extra info on the sequence level, and more detailed on the social attention level
# think about it

# first filter to relevant part of dataset 
# which is all sequences where there were capuchins present (so not "none" for social attention)
soc_att <- detseq[detseq$socatt != "None",]
nas <- soc_att[!complete.cases(soc_att),]

# maybe work with behavior-level dataframe?
socatt_l <- dettools_r2[which(dettools_r2$sequenceID %in% soc_att$sequenceID),]

# load in social attention coding of these sequences
socatt_c <- read.csv("detailedtools/socialattentioncoding.csv")
# exclude "double" coded sequences where BAL processed two at the same time
socatt_c <- socatt_c[!str_detect(socatt_c$Observation.id, "double") == TRUE,]
head(socatt_c)
# clean similarly to detailed tool output
socatt_c <- socatt_c[order(socatt_c$Observation.id),]
socatt_c <- data.frame("videoID" = socatt_c$Observation.id, "codingdate" = socatt_c$Observation.date,
                       "medianame" = socatt_c$Media.file.name, "videolength" = socatt_c$Media.duration..s., "coder" = 
                         socatt_c$Coder.ID, "subjectID" = socatt_c$Subject, "behavior" = socatt_c$Behavior,
                       "modifier1" = socatt_c$Modifier..1,  
                       "starttime" = socatt_c$Start..s., "comment" = socatt_c$Comment.start)

# need to assign correct sequence ID to each line.
socatt_c$sequenceID <- NA
# it's easy to match up when there is just one sequence in the video, so flag it and do those first
vidfreq <- as.data.frame(ftable(soc_att$videoID)) # see how many times each videoID occurs
# filter just to multiples
multiples <- as.character(vidfreq$Var1[vidfreq$Freq > 1])
soc_att$multiple <- ifelse(soc_att$videoID %in% multiples, 1, 0)
socatt_c$multiple <- ifelse(socatt_c$videoID %in% multiples, 1, 0)
#dataset with singles only
socatt_cs <- socatt_c[socatt_c$multiple == 0,]

for(i in 1:nrow(socatt_cs)) {
    socatt_cs$sequenceID[i] <- soc_att$sequenceID[soc_att$videoID == socatt_cs$videoID[i]]
}

# now all that's left is to match the remaining lines to a sequenceID
# I think this will be easiest using the behavior level coding
# filter to just the NAs
socatt_cm <- socatt_c[socatt_c$multiple == 1,]

for(i in 1:nrow(socatt_cm)) {
  # pull for this video the other coded behaviors from the detailed tool use coding
  thisvideo <- dettools_r2[dettools_r2$videoID == socatt_cm$videoID[i],]
  # find behavior with closest starttime to each NA row
  thisvideo$difftimes <- abs(thisvideo$starttime - socatt_cm$starttime[i])
  # then take the sequenceID of that behavior
  socatt_cm$sequenceID[i] <- thisvideo$sequenceID[thisvideo$difftime == min(thisvideo$difftimes)]
}

# CHECK IS THIS WORKED OK OR IF THERE ARE OBVIOUS MISTAKES TO FIX
# get behavior-level dataframe to same columns as soc_att_c
# filter out ones we already have sequence ID for 
socatt_l <- socatt_l[socatt_l$videoID %in% socatt_cm$videoID, names(socatt_l) %in% names(socatt_cm)]

# just attach and then reorder
socatt_lc <- rbind(socatt_l, socatt_cm[, -12])
socatt_lc <- socatt_lc[order(socatt_lc$videoID, socatt_lc$starttime),]

#write.csv(socatt_lc, "detailedtools/socialattentioncheck.csv")

# later do more sanity checks
# for now make one dataframe with all the social attention coding with sequence IDs
socatt_ct <- rbind(socatt_cm, socatt_cs)
sum(is.na(socatt_ct$sequenceID)) # no NAs left

# add age sex of named individuals
socatt_ct <- left_join(socatt_ct, capID, by = c("subjectID" = "ID"))
socatt_ct$agesex <- ifelse(is.na(socatt_ct$Age) == TRUE, socatt_ct$subjectID, 
                          paste(socatt_ct$Age, socatt_ct$Sex, ""))
socatt_ct$agesex <- tolower(gsub(" ", "", socatt_ct$agesex))

# this is a nice overview of all the exact timings of the behaviors etc.
# could envision adding this into the main coding so it goes in between the tool use behavior
# for now will extract the important information per sequence to use in analyses
# make sequence level dataframe
socatt_seq <- socatt_ct[!duplicated(socatt_ct$sequenceID), !names(socatt_ct) %in% c("subjectID","behavior", 
                                                                                     "modifier1", "starttime",
                                                                                     "comment")]
head(socatt_seq)

# Clean to information we want
# who is paying social attention?
socatt_total <- dcast(setDT(socatt_ct[socatt_ct$behavior == "socialattention", c("sequenceID", "subjectID")]), 
      sequenceID ~ rowid(sequenceID, prefix = "socatt_ID"),
      value.var = "subjectID")

# how many individuals are paying social attention
socatt_total$n_socatt <- rowSums(!is.na(socatt_total[,2:4]))

socatt_seq <- left_join(socatt_seq, socatt_total, by = "sequenceID")
# timing of social attention (during tool use or after, falling within sequence start-end or after)
# this is a possibility for the future, but I am leaving this out for now

# also have per age sex class how many were paying social attention (to be able to divide social attention by opportunity)
socatt_as <- as.data.frame(as.matrix(ftable(socatt_ct$sequenceID[which(socatt_ct$behavior == "socialattention")], 
                                            socatt_ct$agesex[which(socatt_ct$behavior == "socialattention")])))
colnames(socatt_as) <- c("sa_nAF", "sa_nAM", "sa_nJM", "sa_nJU", "sa_nSM")
socatt_as$sequenceID <- rownames(socatt_as)
socatt_as$sa_nJuveniles <- socatt_as$sa_nJM + socatt_as$sa_nJU
socatt_as$sa_nAdults <- socatt_as$sa_nAF + socatt_as$sa_nAM
socatt_as$sa_nSubadults <- socatt_as$sa_nSM
head(socatt_as)

socatt_seq <- left_join(socatt_seq, socatt_as, by = "sequenceID")

# who is displacing
disp_only <- dcast(setDT(socatt_ct[socatt_ct$behavior == "displace", c("sequenceID", "subjectID")]), 
                      sequenceID ~ rowid(sequenceID, prefix = "disp_ID"),
                      value.var = "subjectID")
disp_only$n_disp <- rowSums(!is.na(disp_only[,2:3]))

socatt_seq <- left_join(socatt_seq, disp_only, by = "sequenceID")

# who is scrounging
scr_IDs <- dcast(setDT(socatt_ct[socatt_ct$behavior == "scrounge", c("sequenceID", "subjectID")]), 
                   sequenceID ~ rowid(sequenceID, prefix = "scr_ID"),
                   value.var = "subjectID")
scr_IDs$n_scr <- rowSums(!is.na(scr_IDs[,2:5]))

socatt_seq <- left_join(socatt_seq, scr_IDs, by = "sequenceID")

# possibly can also differentiate the scrounging type and add this information
# differentiate tolerated and afterwards
socatt_ct$scrounging <- paste(socatt_ct$modifier1, socatt_ct$behavior)
# scrounging that is neither tolerated nor afterwards is stealing
socatt_ct$scrounging <- ifelse(socatt_ct$scrounging == "tolerated scrounge", "tol_scrounge",
                               ifelse(socatt_ct$scrounging == "afterwards scrounge", "aft_scrounge",
                                      ifelse(socatt_ct$scrounging == "None scrounge", "steal_scrounge", socatt_ct$scrounging)))
scrounger_mod <- as.data.frame(as.matrix(ftable(socatt_ct[socatt_ct$behavior == "scrounge", c("sequenceID", "scrounging")])))
head(scrounger_mod)
scrounger_mod$sequenceID <- rownames(scrounger_mod)

socatt_seq <- left_join(socatt_seq, scrounger_mod, by = "sequenceID")

# who is present
presence <- socatt_ct[socatt_ct$behavior == "present",]

present_as <- as.data.frame(as.matrix(ftable(presence$sequenceID, presence$agesex)))
colnames(present_as) <- c("p_nAF", "p_nAM", "p_nJM", "p_nJU", "p_nSM", "p_nUU")
present_as$sequenceID <- rownames(present_as)
present_as$p_total <- rowSums(present_as[,1:6])
present_as$p_nJuveniles <- present_as$p_nJM + present_as$p_nJU
present_as$p_nAdults <- present_as$p_nAF + present_as$p_nAM
present_as$p_nSubadults <- present_as$p_nSM
head(present_as)

socatt_seq <- left_join(socatt_seq, present_as, by = "sequenceID")

# set all NA's to 0 except for categorical variables
socatt_seq <- socatt_seq %>%
  mutate_at(vars("n_socatt", "sa_nAF", "sa_nAM", "sa_nJM", "sa_nJU", "sa_nSM", 
                 "sa_nJuveniles", "sa_nAdults", "sa_nSubadults", "n_disp", "n_scr", "p_nAF",
                 "p_nAM", "p_nJM", "p_nJU", "p_nSM", "p_nUU", "p_total", "p_nJuveniles",
                 "p_nAdults", "p_nSubadults", "aft_scrounge", "tol_scrounge", "steal_scrounge"), ~replace_na(.,0))

### Datasets we are now left with ####
# detseq #
# aggregated to one row per sequence, contains all information on efficiency etc
head(detseq)
detseq <- detseq[,c("videoID", "codingdate", "medianame", "videolength", "coder", "subjectID", "seqnumber", "location",
                    "mediadate", "sequenceID", "item", "h_startloc", "h_endloc", "hammerID", "outcome", "displacement", 
                    "socatt", "scrounging", "anviltype", "videostart", "videoend", "seqduration", "n_pounds", "n_miss",
                    "n_flies", "n_hloss", "n_misstotal", "n_itemreposit", "n_hamreposit", "n_peel",
                    "n_reposit", "Age", "Sex", "split", "hammerswitches", "anvilswitches", "age_of", "deployment")]
#saveRDS(detseq, "detailedtools/RDS/detseq.rds")

# dettools_r2 #
# not aggregated, every row is a behavior, for detailed looks at the behavior in the sequences
head(dettools_r2)
dettools_r2 <- dettools_r2[,c("videoID", "codingdate", "medianame", "videolength", "coder", "starttime", "subjectID", "behavior", "comment", 
                              "seqnumber", "location", "mediadate", "sequenceID", "item", "h_startloc", "h_endloc", 
                              "outcome", "displacement", "socatt", "scrounging", "poundtype", "onefoot", "overhead", "onehand", "tailsupport",
                              "mistaketype", "repostype", "hammerID2", "h_switchloc", "anviltype2", 
                              "videostart", "videoend", "seqduration", "seqstart", "seqend", "n_pounds", "n_miss", 
                              "n_flies", "n_hloss", "n_misstotal", "n_itemreposit", "n_hamreposit", "n_peel",
                              "n_reposit", "Age", "Sex", "split", "age_of", "deployment")]
#saveRDS(dettools_r2, "detailedtools/RDS/dettools_r2.rds"

# socatt_ct #
# coding only of sequences with capuchins present
# not aggregated, every row is a behavior, detailing scrounging, displacement, avoidance and social attention
# also coding which age/sex classes and individuals are present during the sequence
head(socatt_ct)
socatt_ct <- socatt_ct[,c("videoID", "codingdate", "medianame", "videolength", "coder", "subjectID","behavior", 
                          "modifier1", "starttime", "comment",  "sequenceID", "multiple", "agesex")]
#saveRDS(socatt_ct, "detailedtools/RDS/socatt_ct.rds")


# socatt_seq #
# aggregated to one row per sequence, summarized information on how many individuals scrounged, paid social attention
# which age sex classes were present, who displaced, etc
head(socatt_seq)
socatt_seq <- socatt_seq[,c("videoID", "codingdate", "medianame", "videolength", "coder","sequenceID", "n_socatt",
                            "n_disp", "n_scr", "p_total", "socatt_ID1", "socatt_ID2", "socatt_ID3", "sa_nAF",
                            "sa_nAM", "sa_nJU", "sa_nSM", "sa_nJuveniles","sa_nAdults", "sa_nSubadults", "disp_ID1",
                            "disp_ID2", "scr_ID1", "scr_ID2", "scr_ID3", "scr_ID4", "p_nAF", "p_nAM", "p_nJM",
                            "p_nJU", "p_nSM", "p_nUU", "p_nJuveniles","p_nJuveniles", "p_nAdults", "p_nSubadults",
                            "aft_scrounge", "tol_scrounge", "steal_scrounge")]
#saveRDS(socatt_seq, "detailedtools/RDS/socatt_seq.rds")

# Generate sample of 100 sequences for interrater reliability
set.seed(22)
IRseqs <- sample_n(detseq[which(detseq$coder != "LR"),], 150)
IRvideos <- unique(dettools_r$videoID[which(dettools_r$sequenceID %in% IRseqs$sequenceID)])
IRvideos

