## Detailed tool use -- Analyses not used in PhD chapter
## MPI-AB; Z Goldsborough

## Analysing other aspects of tool use
# E.g., technique, seasonality etc

## packages needed
library(stringr)
library(viridis)
library(dplyr)
library(tidyr)
library(brms)
library(ggplot2)
library(tidybayes)
library(fitdistrplus)
library(reshape2)
library(lubridate)
library(mgcv)
library(gratia)

### Load datasets, cleaned in detailedtools.R script ####
#detseq <- readRDS("detailedtools/RDS/detseq.rds")
head(detseq)
# every row is a tool use sequence, all the information is aggregated to the level of a single sequence
#dettools_r2 <- readRDS("detailedtools/RDS/dettools_r2.rds")
head(dettools_r2)
# every row is a behavior, so this still contains every coded behavior with timestamp, no aggregation


### Technique ####
## Individual variation in technique
# very basic, bar chart of types of pounds

# work not with detseq_o2 but non-aggregated dataset
str(dettools_r2)
# filter to known individuals and opened, not split, hammerstone not already in hand and only brown almendras? trying to keep most things consistent
dettools_r2oi <- dettools_r2[dettools_r2$subjectID %in% knownids$ID & dettools_r2$outcome == "opened" & dettools_r2$item == "almendrabrown" & dettools_r2$split == FALSE & !dettools_r2$h_startloc == "inhand",]
# also filter out to only the behaviors (so remove hammerstone for sure)
dettools_r2oi <- dettools_r2oi[!dettools_r2oi$behavior %in% c("hammerstone"), ]

ggplot(data = dettools_r2oi[dettools_r2oi$behavior == "pound" & dettools_r2oi$poundtype %in% c("crouch", "stand", "jump"),], aes(x = poundtype, fill = subjectID)) + geom_histogram(stat = "count") + facet_wrap(~ subjectID)
ggplot(data = dettools_r2oi[dettools_r2oi$behavior == "pound" & dettools_r2oi$poundtype %in% c("crouch", "stand", "jump"),], aes(x = onefoot, fill = subjectID)) + geom_histogram(stat = "count") + facet_wrap(~ subjectID)

# make new behavior comment that also includes what type of pound it is 
dettools_r2oi$behavior[which(dettools_r2oi$behavior == "pound")] <- ifelse(dettools_r2oi$poundtype[which(dettools_r2oi$behavior == "pound")] =="stand", "standpound", 
                                                                           ifelse(dettools_r2oi$poundtype[which(dettools_r2oi$behavior == "pound")] == "crouch", "crouchpound", "jumppound"))
dettools_r2oi$behavior[which(dettools_r2oi$behavior == "reposit")] <- ifelse(dettools_r2oi$repostype[which(dettools_r2oi$behavior == "reposit")] == "peel", "peel", "reposit")

ftable(dettools_r2oi$behavior)

## sunburst
require(sunburstR)
#function to be able to suppress the NA's
paste5 <- function(..., sep = " ", collapse = NULL, na.rm = F) {
  if (na.rm == F)
    paste(..., sep = sep, collapse = collapse)
  else
    if (na.rm == T) {
      paste.na <- function(x, sep) {
        x <- gsub("^\\s+|\\s+$", "", x)
        ret <- paste(na.omit(x), collapse = sep)
        is.na(ret) <- ret == ""
        return(ret)
      }
      df <- data.frame(..., stringsAsFactors = F)
      ret <- apply(df, 1, FUN = function(x) paste.na(x, sep))
      
      if (is.null(collapse))
        ret
      else {
        paste.na(ret, sep = collapse)
      }
    }
}

#make sure it doesn't round to the nearest 10
custom.message = "function (d) {
  root = d;
  while (root.parent) {
    root = root.parent
  }
  p = (100*d.value/root.value).toPrecision(3);
  msg = p+' %<br/>'+d.value+' of '+root.value;
  return msg;
}"

cols <- data.frame(col = c("purple",  # anvilswitch
                           "#990F0F",  # crouchpound
                           "#005000", # jumppound
                           "yellow", #misstrike
                           "#008600", # peel
                           "#009292", # reposit
                           "darkblue", # standpound
                           "#e0e2e3", #seqend
                           "orange"), behavior = c("anvilswitch", "crouchpound", "jumppound", "misstrike", "peel", "reposit", "standpound", "seqend", "hammerswitch"))

## make a sunburst, would want to be able to split by adult, subadult, juvenile (I think?)
# need columns: sequenceID, individual, age, behavior, get number of what number of behavior it is in the sequence
sunbursttools <- dettools_r2oi[,c("sequenceID", "subjectID", "Age", "behavior")]

sunbursttools$nr <- ave(sunbursttools$sequenceID, sunbursttools$sequenceID, FUN = seq_along)
hist(as.numeric(sunbursttools$nr))
sunbursttools$nr_n <- as.numeric(sunbursttools$nr)
# so okay some go really long.. Potentially could consider only taking the ones that don't go beyond 10
long_sequences <- unique(sunbursttools$sequenceID[which(sunbursttools$nr_n > 10)])
sunbursttools_s <- subset(sunbursttools, ! sunbursttools$sequenceID %in% long_sequences)
sunbursttools_s <- sunbursttools_s[complete.cases(sunbursttools_s),]

suntoolsC <- dcast(sunbursttools_s, sequenceID ~ nr, value.var = "behavior")

#steps to make sunburst
suntoolsC$First <- suntoolsC$`1`
suntoolsC$Second <- suntoolsC$`2`
suntoolsC$Third <- suntoolsC$`3`
suntoolsC$Fourth <- suntoolsC$`4`
suntoolsC$Fifth <- suntoolsC$`5`
suntoolsC$Sixth <- suntoolsC$`6`
suntoolsC$Seventh <- suntoolsC$`7`
suntoolsC$Eighth <- suntoolsC$`8`
suntoolsC$Ninth <- suntoolsC$`9`
suntoolsC$Tenth <- suntoolsC$`10`

suntoolsC$sequence <- with(suntoolsC, paste5(Second, Third, Fourth, Fifth, Sixth, Seventh, Eighth, Ninth, Tenth, sep = "-", na.rm = T))

# general sunburst
sunburst(data.frame(table(suntoolsC$sequence)), colors = cols, explanation = custom.message)

# still need to make colors the same across sunbursts!! 
# https://stackoverflow.com/questions/49993198/how-to-specify-the-colors-and-toggle-labels-for-each-category-in-r-sunburst
# https://community.plotly.com/t/sunburst-color-levels-r/33253/5
# https://stackoverflow.com/questions/70246083/how-to-setup-a-color-per-category-accross-all-layers-of-a-sunburst-plotly-graph

# adults
adultsequences <- unique(sunbursttools_s$sequenceID[which(sunbursttools_s$Age == "Adult")])
suntoolsC_adult <- suntoolsC[suntoolsC$sequenceID %in% adultsequences,]

sunburst(data.frame(table(suntoolsC_adult$sequence)), explanation = custom.message)

# subadults
subadultsequences <- unique(sunbursttools_s$sequenceID[which(sunbursttools_s$Age == "Subadult")])
suntoolsC_subadult <- suntoolsC[suntoolsC$sequenceID %in% subadultsequences,]

sunburst(data.frame(table(suntoolsC_subadult$sequence)), explanation = custom.message)

# juveniles
juvsequences <- unique(sunbursttools_s$sequenceID[which(sunbursttools_s$Age == "Juvenile")])
suntoolsC_juv <- suntoolsC[suntoolsC$sequenceID %in% juvsequences,]

sunburst(data.frame(table(suntoolsC_juv$sequence)), explanation = custom.message)

## potentially could look only at the pound progression (so only have jumppound, standpound, crouchpound)
# but this does look like a good visual tool for comparison (once I have the same behaviors in the same colors across sunbursts)

### Hammerstones #####

# how often do they transport hammerstone in?
# on all sequences, including those that didnt finish
ftable(detseq$h_startloc)
ggplot(detseq, aes(x = h_startloc, fill = h_startloc)) + geom_histogram(stat = "count") + theme_bw() + facet_wrap(~age_of)
ggplot(detseq, aes(x = h_endloc, fill = h_endloc)) + geom_histogram(stat = "count") + theme_bw() + facet_wrap(~age_of)

# so mostly juveniles transport hammers in, but they are also the ones who usually leave the hammerstone not on the anvil
# so could be artefact of repeat sequences where they dont put the hammerstone back properly so have to fetch it for the new sequence

## Hammerstone identities ##
# what hammerstone is used to process what item?
# what is the average number of pounds per hammerstone?

# filter to opened sequences, with hammerstone IDs known
detseq_oh <- detseq_o[detseq_o$hammerID %in% c("BAM", "BCH", "DPL", "DWA", "DWA_A", "DWA_B", "FRE", "LCH", "PEB", "BRK", "BOA", "BOA_A") & detseq_o$split == FALSE,]

# number of pounds per hammerstone
ftable(detseq_oh$hammerID)
ggplot(detseq_oh[detseq_oh$Age == "Adult" | detseq_oh$Age == "Subadult",], aes(x=hammerID, y=n_pounds, fill = location)) + 
  geom_violin() + theme_bw() + facet_wrap(~item, scales = "free_x")

# what items are being opened with what hammerstone
ggplot(detseq[detseq$hammerID %in% c("BAM", "BCH", "DPL", "DWA", "DWA_A", "DWA_B", "FRE", "LCH", "PEB", "BRK", "BOA"),], aes(x = item, fill = item)) + geom_histogram(stat = "count") + theme_bw() + facet_wrap(~hammerID, scales = "free_x")

## hammer timeline
ggplot(detseq_oh[detseq_oh$location == "EXP-ANV-01",], 
       aes(x = mediadate, fill = hammerID)) + geom_histogram() + theme_bw() 
ggplot(detseq_oh[detseq_oh$location == "CEBUS-02",], 
       aes(x = mediadate, fill = hammerID)) + geom_histogram() + theme_bw() 

# are some individuals contributing a lot to the accumulation at a site
ggplot(detseq[detseq$subjectID %in% unique(detseq_oi$subjectID),], 
       aes(x = subjectID, fill = Age)) + geom_histogram(stat = "count") +
  theme_bw() + facet_grid(cols = vars(location), rows = vars(deployment))

# individuals preferences for hammerstones
ggplot(detseq[detseq$subjectID %in% unique(detseq_oi$subjectID) & detseq$hammerID %in% c("BAM", "BCH", "DPL", "DWA", "DWA_A", "DWA_B", "FRE", "LCH", "PEB", "BRK", "BOA", "BOA_A") & detseq$location == "EXP-ANV-01" & detseq$deployment == "R11",], 
       aes(x = subjectID, fill = Age)) + geom_histogram(stat = "count") +
  theme_bw() + facet_wrap(~hammerID)

ggplot(detseq[detseq$subjectID %in% unique(detseq_oi$subjectID) & detseq$hammerID %in% c("BAM", "BCH", "DPL", "DWA", "DWA_A", "DWA_B", "FRE", "LCH", "PEB", "BRK", "BOA", "BOA_A") & detseq$location == "CEBUS-02",], 
       aes(x = subjectID, fill = Age)) + geom_histogram(stat = "count") +
  theme_bw() + facet_wrap(~hammerID)

# descriptives
ftable(detseq_o$socatt)
ftable(detseq_o$displacement)
ftable(detseq_o$scrounging)

ftable(dettools_r2$mistaketype)
ftable(dettools_r2$repostype)

##### seasonality #######
### What is being processed when? ####

# now looking only at EXP-ANV-01-R11 as that is the only one fully coded
# later could look at all sites/times
# for this plot exclude really rare items
ggplot(detseq[detseq$location == "CEBUS-02" & 
                detseq$item %in% c("almendrabrown", "almendragreen", "almendraunknown", "almendrared"),], 
       aes(x = mediadate, fill = item)) + geom_histogram() + theme_bw() + facet_wrap(~item)

## GAM model
# multinomial model (?) item being processed depending on day of the year*location interaction
# make month variable
detseq$month <- month(detseq$mediadate)
# for now don't have enough years, but later could include year
detseq$year <- year(detseq$mediadate)

detseq_gam <- detseq[detseq$item %in% c("almendrabrown", "almendragreen", "almendraunknown", "almendrared"),]
detseq_gam$itemF <- as.factor(detseq_gam$item)
detseq_gam$locationF <- as.factor(detseq_gam$location)

# brms
# smooth version 
alm_bm1 <- brm(itemF ~ s(month, bs ="cc", k = 11, by = locationF) + locationF, data=detseq_gam, family="categorical", 
               knots = list(month = c(0.5,12.5)), chains=2, cores = 4, backend = "cmdstanr", save_pars = save_pars(all = TRUE),
               iter = 1000)

# saveRDS(alm_bm1, file = "detailedtools/RDS/alm_bm1.rds")
# alm_bm1 <- readRDS("detailedtools/RDS/alm_bm1.rds")

summary(alm_bm1)
mcmc_plot(alm_bm1)
plot(conditional_smooths(alm_bm1, categorical = TRUE))
plot(conditional_effects(alm_bm1, categorical = TRUE))

# best plot
conditions <- make_conditions(alm_bm1, "locationF")
alm_plot <- plot(conditional_effects(alm_bm1, categorical = TRUE, conditions = conditions), plot = FALSE)[[2]]
alm_plot + theme_bw()

