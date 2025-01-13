## Detailed tool use analyses -- Cleaning Script
## MPI-AB; Z Goldsborough

## Script how to clean BORIS output data

## packages needed
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)

## NOTE: add code about catching when anvil type switches (e.g. if it's wooden anvil at comment seq_start, or when they switch back to stone/wood)
## NOTE: check how to extract end location of second hammer when they switch (code hammerend loc again?)

### Loading dataset ####
# load csv files with aggregated BORIS output (from BORIScoding/ExportedData Google Drive)
# Zoë's csv
dettools1 <- read.csv("detailedtools/ZGdetailedtoolscoding.csv")
# Meredith's csv
dettools2 <- read.csv("detailedtools/EXP-ANV-01-R11_MC.csv")
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

## add age/sex to IDS
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
# fix occurrences of DWA after fracture (2022-03-11), should be DWA_A
dettools_r2$hammerID[which(dettools_r2$mediadate > "2022-03-11" & dettools_r2$hammerID == "DWA")] <- "DWA_A"
dettools_r2$hammerID2[which(dettools_r2$mediadate > "2022-03-11" & dettools_r2$hammerID2 == "DWA")] <- "DWA_A"

dettools_r2 <- left_join(dettools_r2, seqdat[,c("sequenceID", "split")], by = "sequenceID")

# sequence level dataframe (only for analyzing things like nr of pounds/duration. things that are fixed per sequence)
detseq <- dettools_r2[!duplicated(dettools_r2$sequenceID),]
detseq <- left_join(detseq, seqdat[,c("sequenceID", "hammerswitches", "anvilswitches")], by = "sequenceID")

# add a pound if the hammerstone was "inhand" or if there is a split
detseq$n_pounds[which(detseq$h_startloc == "inhand" | detseq$split == TRUE)] <- detseq$n_pounds[which(detseq$h_startloc == "inhand" | detseq$split == TRUE)] + 1

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

# Generate sample of 100 sequences for interrater reliability
set.seed(22)
IRseqs <- sample_n(detseq[which(detseq$coder != "LR"),], 150)
IRvideos <- unique(dettools_r$videoID[which(dettools_r$sequenceID %in% IRseqs$sequenceID)])
IRvideos

