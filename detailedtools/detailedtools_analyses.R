## Detailed tool use -- Analyses
## MPI-AB; Z Goldsborough

## Analysing efficiency of tool use and social attention paid to tool use

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
library(easystats)
library(scales)
library(ggh4x)
library(ggeffects)

## NOTES ###
# If we have enough data, ideally exclude sequences that are split across multiple videos (split == TRUE) and of which outcome is unknown
# Similarly, if we have enough data with identified individuals, only work on that dataset so individual can be a random effect

### Load datasets, cleaned in detailedtools.R script ####
#detseq <- readRDS("detailedtools/RDS/detseq.rds")
head(detseq)
# every row is a tool use sequence, all the information is aggregated to the level of a single sequence
#dettools_r2 <- readRDS("detailedtools/RDS/dettools_r2.rds")
head(dettools_r2)
# every row is a behavior, so this still contains every coded behavior with timestamp, no aggregation

### Filter and general diagnostics ####
## Efficiency coding ##
# how much data do we have?
nrow(detseq)
ftable(detseq$item)
ftable(detseq$outcome)
ftable(detseq$subjectID)
mean(ftable(detseq$mediadate))
length(unique(detseq$mediadate))
ftable(detseq$location)
ftable(detseq$scrounging)
ftable(detseq$displacement)
ftable(detseq$socatt)
ftable(detseq$hammerID, detseq$location)

# plot of sequences per location
# png("detailedtools/RDS/locationscoverage.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(detseq, 
       aes(x = mediadate)) + geom_histogram() + facet_wrap(~location) + theme_bw() +
  labs(x = "Date", y = "Number of tool use sequences") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) + theme(axis.text.x = element_text(angle =45, hjust = 1))
# dev.off()

# only include tool use sequences that were successful (so outcome = opened) to be able to measure efficiency
# consider variability in the items processed
ftable(detseq$item)
# see that anything other than almendra (sea almond) is very rare
# also likely need different techniques and nr of pounds for different items
# so subset to almendra's only (all levels of ripeness)
detseq_o <- detseq[detseq$outcome == "opened" & detseq$item %in% c("almendrabrown", "almendragreen", "almendraunknown", "almendrared"),]

# assess quantity and quality of data
ggplot(detseq_o, aes(x=coder, y=n_pounds)) + geom_violin()
table(detseq_o$location, detseq_o$deployment)
ftable(detseq_o$subjectID) # for vast majority of sequences we can identify the tool user
ftable(detseq_o$Age) # most sequences are by subadults, but good range of the different age classes
nrow(detseq_o)
# so at the moment have 3203 sequences

# subset to the best of the best: known individuals, and no split videos
# Split videos (spanning several videos) likely miss pounds, and as such are not reliable measures of efficiency
ftable(detseq_o$split) # exclude 254 split sequences 
ftable(detseq_o$subjectID) # exclude 59 not ID'd sequences 
detseq_oi <- detseq_o[!detseq_o$subjectID %in% c("adultmale", "subadultmale", "juvenileunknown") & detseq_o$split == FALSE,]
nrow(detseq_oi)
ftable(detseq_oi$Age)
# then end up with 2890 sequences

# make age a factor in the right order (for plotting)
detseq_oi$Age <- factor(detseq_oi$Age, levels = c("Juvenile", "Subadult", "Adult"))


## Social attention ##
## TEMPORARY CHECK WHAT STILL TO CODE
soc_att <- detseq[!detseq$socatt == "None",]
socatt_vidnames <- soc_att[,c("videoID", "coder", "subjectID", "socatt", "scrounging", "displacement")]
# filter to ones not coded yet
# load in coding
socatt_c <- read.csv("detailedtools/socialattentioncoding.csv")
tocode <- socatt_vidnames[!socatt_vidnames$videoID %in% socatt_c$Observation.id,]
tocode[0:nrow(tocode),]

## social attention coding dataframes, two levels
head(socatt_seq) # aggregated to sequence
head(socatt_ct) # every row a behavior occurrence

# add variable whether there is social attention in the sequence
socsequences <- unique(socatt_ct$sequenceID[which(socatt_ct$behavior == "socialattention")])
socatt_ct$socatt <- ifelse(socatt_ct$sequenceID %in% socsequences, 1, 0)
socatt_cts <- socatt_ct[socatt_ct$behavior == "socialattention",]

#make dataframe with already all presence without social attention
socatt_final <- socatt_ct[socatt_ct$behavior == "present" & socatt_ct$socatt == 0,]

# work on the sequence level
for(i in 1:length(socsequences)){
  observers <- socatt_cts$subjectID[which(socatt_cts$sequenceID == socsequences[i])]
  seq_present <- socatt_ct[socatt_ct$sequenceID == socsequences[i] & socatt_ct$behavior == "present",]
  # set all socatt to 0 except for the individual paying social attention 
  seq_present$socatt <- 0
  #easiest to match for each observer
  # for first observer, just take first line with matching presence 
  # (e.g. if there are multiple juveniles just take the first, and make it a 1)
  seq_present$socatt[which(seq_present$subjectID == observers[1])][1] <- 1
  if(length(observers) > 1) {
    # to avoid making the same one 1 again, subset to those still 0
    seq_present$socatt[which(seq_present$subjectID == observers[2] & seq_present$socatt == 0)][1] <- 1
  }
  if(length(observers) > 2) {
    seq_present$socatt[which(seq_present$subjectID == observers[3] & seq_present$socatt == 0)][1] <- 1
  }
  socatt_final <- rbind(socatt_final, seq_present)
}

socatt_final <- socatt_final[order(socatt_final$sequenceID),]
# in this sample, we have:
length(unique(socatt_final$sequenceID))
# 992 sequences with capuchins present
# now attach the relevant information from the main dataframe (information on the tool user etc)
head(detseq)
detseq$socialattention <- detseq$socatt

socatt_final <- left_join(socatt_final, detseq[,c("sequenceID", "subjectID", "coder", "location", "item",
                                                  "outcome", "displacement", "scrounging", "socialattention",
                                                  "anviltype", "seqduration", "n_pounds", "n_misstotal",
                                                  "Age", "deployment", "split", "hammerswitches", "anvilswitches")],
                          by = "sequenceID")

# also add more detailed information on number of capuchins present, number scrounging etc, from the social coding
head(socatt_seq)
socatt_final <- left_join(socatt_final, socatt_seq[,c("sequenceID", "n_socatt", "n_disp", "n_scr", "p_total")],
                          by = "sequenceID")

######### TEMPORARY UNTIL ALL CODING IS DONE
# CHECK IF SCORUNGING/DISPLACEMENT/SOCATT TRACKED
table(socatt_final$socialattention, socatt_final$n_socatt)
# generate list of incongruities to fix
checklist <- unique(socatt_final$sequenceID[socatt_final$n_disp == 0 & socatt_final$displacement == "fulldisp"| 
                                             socatt_final$n_disp > 0 & socatt_final$displacement == "None" |
                                             socatt_final$n_disp >0 & socatt_final$displacement == "nodisplacement"|
                                       socatt_final$n_scr == 0 & socatt_final$scrounging == "scrounging" |
                                       socatt_final$n_scr > 0 & socatt_final$scrounging == "None" |
                                       socatt_final$n_scr > 0 & socatt_final$scrounging == "noscrounging" |
                                       socatt_final$n_socatt >0 & socatt_final$socialattention == "noattention"|
                                         socatt_final$n_socatt ==0 & socatt_final$socialattention == "socialattention"])
# check the checklist files, first in the social attention coding to see if that is correct
# if that is correct, then these need to be changed in the original BORIS coding
# change the ones I can, and have Leonie changes the ones in her and Meredith's coding
#############################

# clean up the final dataframe
head(socatt_final)
# subjectID.x = observer, coder.x = coder of social attention
socatt_final$observerID <- socatt_final$subjectID.x
socatt_final$tooluserID <- socatt_final$subjectID.y
socatt_final$coder_socatt <- socatt_final$coder.x
socatt_final$coder_tooluse <- socatt_final$coder.y
socatt_final$observer_agesex <- ifelse(str_detect(socatt_final$agesex, "juvenile"), "juvenile", socatt_final$agesex)
socatt_final$tooluser_age <- socatt_final$Age
# variable for if there is social attention in sequence
socatt_final$socialattention <- ifelse(socatt_final$sequenceID %in% socsequences, "socialattention", "nosocialattention")

socatt_final <- socatt_final[,c("sequenceID", "videoID", "deployment", "location", "anviltype", "seqduration", "coder_socatt", 
                                "coder_tooluse", "split", "item", "observerID", "observer_agesex", "socatt","tooluserID", "tooluser_age",
                                "n_pounds", "n_misstotal", "hammerswitches", "anvilswitches", "n_socatt", "n_disp", "n_scr", "p_total",
                                "outcome", "displacement", "scrounging", "socialattention")]

# some descriptives
# how many of these sequences are split across various videos (so information missing)
ftable(socatt_final$split[!duplicated(socatt_final$sequenceID)])
table(socatt_final$socialattention[!duplicated(socatt_final$sequenceID)], 
      socatt_final$split[!duplicated(socatt_final$sequenceID)])
# how many sequences have unknown observers
unique(socatt_final$videoID[which(socatt_final$observer_agesex == "unknown")])
# what are outcomes depending on social attention yes/no
table(socatt_final$socialattention[!duplicated(socatt_final$sequenceID)], 
      socatt_final$outcome[!duplicated(socatt_final$sequenceID)])
# unknown means that we do not know the ending because it is not on camera 
# none could be unknown or another reason
# either way, good to exclude these because we missed part of the sequence and 
# therefore are not sure who was present/who paid attention

# exclude these sequences
socatt_final <- socatt_final[socatt_final$split == FALSE & !socatt_final$observer_agesex == "unknown" &
                               !socatt_final$outcome == "None" & !socatt_final$outcome == "Unknown",]
# then have 839 sequences
length(unique(socatt_final$sequenceID))



###
### EFFICIENCY ####
###

## Comparing efficiency between age classes on opened sequences, multiple measures
## Other factors that likely affect efficiency are:
# The ripeness of the sea almond
# The type of anvil (wood or stone)
# The identity of the tool user

##### 1. Sequence duration (seconds) #####
## determine distribution
descdist(detseq_oi$seqduration)
hist(detseq_oi$seqduration)

testdist1.1 <- fitdist(detseq_oi$seqduration, "norm")
plot(testdist1.1)

testdist1.2 <- fitdist(detseq_oi$seqduration, "gamma")
plot(testdist1.2)
# Gamma appears best 

### Model_e1 ###
# Outcome: sequence duration (seconds)
# Fixed effects: age, item (ripeness of sea almond) and anviltype (wood or stone) 
# Random effects: subjectID (identity of tool user)
m_e1 <- brm(seqduration ~ Age + item + anviltype + (1|subjectID), data = detseq_oi, iter = 2000, 
            save_pars = save_pars(all = TRUE), chain = 2, core = 2, backend = "cmdstanr", 
            control = list(adapt_delta = 0.99), family = "gamma")
# m_e1 <- add_criterion(m_e1, c("loo", "loo_R2", "bayes_R2"), reloo = TRUE, backend = "cmdstanr", ndraws = 2000) 

# saving model to not have to run it again
# saveRDS(m_e1, "detailedtools/RDS/m_e1.rds")
# m_e1 <- readRDS("detailedtools/RDS/m_e1.rds")

# Diagnostics
summary(m_e1)
mcmc_plot(m_e1)
pp_check(m_e1)
plot(conditional_effects(m_e1))

loo(m_e1) # had one influential case (which is observation 265, Ink who is adult but bad)
loo_R2(m_e1) # 0.22
round(bayes_R2(m_e1),2) # 0.23

plot(m_e1$criteria$loo, label_points = TRUE)

# Interpretation
round(exp(2.99),2)

hypothesis(m_e1, "Intercept  > Intercept + AgeSubadult", alpha = 0.05)
hypothesis(m_e1, "Intercept  > Intercept + AgeAdult", alpha = 0.05)
hypothesis(m_e1, "Intercept + AgeAdult  < Intercept + AgeSubadult", alpha = 0.05)
hypothesis(m_e1, "Intercept + itemalmendragreen > Intercept", alpha = 0.05)

# report(m_e1)

# Visualization
m_type_pred <- m_e1 %>% 
  epred_draws(newdata = tibble(item = detseq_oi$item,
                               Age = detseq_oi$Age,
                               anviltype = detseq_oi$anviltype,
                               subjectID = detseq_oi$subjectID))

# Age difference in duration to open item
# png("detailedtools/RDS/m_e1_age.png", width = 8, height = 7, units = 'in', res = 300)
ggplot(data = m_type_pred, aes(x = Age, y = .epred)) + geom_violin(aes(color = Age, fill = Age), alpha = 0.4)  +
  stat_summary(detseq_oi, inherit.aes = FALSE, mapping=aes(x = Age, y = seqduration, color = Age), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Age", y = "Seconds required to open item") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14))
# dev.off()

# Sea almond ripeness difference in duration to open item
# png("detailedtools/RDS/m_e1_item.png", width = 8, height = 7, units = 'in', res = 300)
ggplot(data = m_type_pred, aes(x = item, y = .epred)) + geom_violin(aes(color = item, fill = item), alpha = 0.4) + ylim(0,100) +
  stat_summary(detseq_oi, inherit.aes = FALSE, mapping=aes(x = item, y = seqduration, color = item), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Item type", y = "Seconds required to open item") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 
# dev.off()


##### 2. Number of pounds ######
descdist(detseq_oi$n_pounds)

testdist2.1 <- fitdist(detseq_oi$n_pounds, "pois")
plot(testdist2.1)
# Poisson is best family

### Model_e2 ###
# Outcome: number of pounds
# Fixed effects: age, item, and anviltype
# Random effects: subjectID
m_e2 <- brm(n_pounds ~ Age + item + anviltype + (1|subjectID), data = detseq_oi, family = "poisson", 
            iter = 2000, chain = 2, core = 2, save_pars = save_pars(all = TRUE), 
            control = list(adapt_delta = 0.99), backend = "cmdstanr")
# m_e2 <- add_criterion(m_e2, c("loo", "loo_R2", "bayes_R2"), reloo = TRUE, backend = "cmdstanr", ndraws = 2000) 

# save and load model
# saveRDS(m_e2, "detailedtools/RDS/m_e2.rds")
# m_e2 <- readRDS("detailedtools/RDS/m_e2.rds")

# Diagnostics
summary(m_e2)
pp_check(m_e2)
mcmc_plot(m_e2)
plot(conditional_effects(m_e2))

loo(m_e2) # all cases good
loo_R2(m_e2) # 0.12
round(bayes_R2(m_e2),2) # 0.12

# Interpretation
round(exp(1.64),2)
hypothesis(m_e2, "Intercept + itemalmendrared > Intercept", alpha = 0.05)
hypothesis(m_e2, "Intercept > Intercept + AgeAdult", alpha = 0.05)
hypothesis(m_e2, "Intercept > Intercept + anviltypewood", alpha = 0.05)


# Visualization
m_type_pred2 <- m_e2 %>% 
  epred_draws(newdata = tibble(item = detseq_oi$item,
                               Age = detseq_oi$Age,
                               anviltype = detseq_oi$anviltype,
                               subjectID = detseq_oi$subjectID))

# Age difference in number of pounds to open item
# png("detailedtools/RDS/m_e2_pound.png", width = 8, height = 7, units = 'in', res = 300)
ggplot(data = m_type_pred2, aes(x = Age, y = .epred)) + geom_violin(aes(color = Age, fill = Age), alpha = 0.4) +
  stat_summary(detseq_oi, inherit.aes = FALSE, mapping=aes(x = Age, y = n_pounds, color = Age), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Age", y = "Number of pounds required to open item") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14))
# dev.off()

# Item difference in number of pounds to open item
ggplot(data = m_type_pred2, aes(x = item, y = .epred)) + geom_violin(aes(color = item, fill = item), alpha = 0.4) +
  stat_summary(detseq_oi[which(detseq_oi$item %in% c("almendrabrown", "almendragreen")),], inherit.aes = FALSE, mapping=aes(x = item, y = n_pounds, color = item), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Item type", y = "Number of pounds required to open item") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

# relationship between nr of pounds and sequence duration
# colors for juvenile, subadult, adult
cols <- viridis(3, option = "plasma", end = 0.8)

# png("detailedtools/RDS/duration_pound.png", width = 8, height = 7, units = 'in', res = 300)
ggplot(detseq_oi, aes(y = seqduration, x = n_pounds, color = Age, shape = Age)) + geom_point(size = 3, alpha = 0.3) + geom_smooth() + 
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  labs(y = "Seconds needed to open item", x = "Number of pounds needed to open item") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 
# dev.off()

##### 3. Number of repositions ######

## 3a: repositions of item

### Model e3a ### 
# Outcome: number of repositions
# Fixed effects: age, item, anviltype
# Random effects: subjectID
m_e3a <- brm(n_itemreposit ~ Age + item + anviltype + (1|subjectID), data = detseq_oi, family = "poisson", 
             iter = 2000, chain = 2, core = 2, control = list(adapt_delta = 0.99),
             save_pars = save_pars(all = TRUE), backend = "cmdstanr")
# m_e3a <- add_criterion(m_e3a, c("loo", "loo_R2", "bayes_R2"), reloo = TRUE, backend = "cmdstanr", ndraws = 2000) 

# saving and loading model
# saveRDS(m_e3a, "detailedtools/RDS/m_e3a.rds")
# m_e3a <- readRDS("detailedtools/RDS/m_e3a.rds")

# diagnostics
summary(m_e3a)
mcmc_plot(m_e3a)
pp_check(m_e3a)
plot(conditional_effects(m_e3a))

loo(m_e3a) # all cases good
round(loo_R2(m_e3a),2) # 0.15
round(bayes_R2(m_e3a),2) # 0.17 

# Interpretation
round(exp(0.20),2)
hypothesis(m_e3a, "Intercept > Intercept + AgeAdult", alpha = 0.05)
hypothesis(m_e3a, "Intercept > Intercept + AgeSubadult", alpha = 0.05)

# make violin plot
m_type_pred3 <- m_e3a %>% 
  epred_draws(newdata = tibble(Age = detseq_oi$Age,
                               item = detseq_oi$item,
                               anviltype = detseq_oi$anviltype,
                               subjectID = detseq_oi$subjectID))

# age difference in number of item repositions
ggplot(data = m_type_pred3, aes(x = Age, y = .epred)) + geom_violin(aes(color = Age, fill = Age), alpha = 0.4) +
  stat_summary(detseq_oi, inherit.aes = FALSE, mapping=aes(x = Age, y = n_itemreposit, color = Age), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Age", y = "Average number of repositions per sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

## 3b: peeling

### Model_3b ### 
# Outcome: number of peels
# Fixed effects: age, item, anviltype
# Random effects: subjectID
m_e3b <- brm(n_peel ~ Age + item + anviltype + (1|subjectID), data = detseq_oi,
             save_pars = save_pars(all = TRUE), family = "poisson", iter = 2000,
             chain = 2, core = 2, backend = "cmdstanr")
# m_e3b <- add_criterion(m_e3b, c("loo", "loo_R2", "bayes_R2"), reloo = TRUE, backend = "cmdstanr", ndraws = 2000) 

# saving and loading model
# saveRDS(m_e3b, "detailedtools/RDS/m_e3b.rds")
# m_e3b <- readRDS("detailedtools/RDS/m_e4b.rds")

# diagnostics
summary(m_e3b)
mcmc_plot(m_e3b)
pp_check(m_e3b)
plot(conditional_effects(m_e3b))

loo(m_e3b) # all cases good
round(loo_R2(m_e3b),2) # 0.10
round(bayes_R2(m_e3b),2) # 0.12 

# Interpretation
round(exp(0.35),2)
hypothesis(m_e3b, "Intercept > Intercept + AgeSubadult", alpha = 0.05)
hypothesis(m_e3b, "Intercept + itemalmendragreen > Intercept", alpha = 0.05)

# make violin plot
m_type_pred3b <- m_e3b %>% 
  epred_draws(newdata = tibble(Age = detseq_oi$Age,
                               item = detseq_oi$item,
                               anviltype = detseq_oi$anviltype,
                               subjectID = detseq_oi$subjectID))

# age difference in number of peels
ggplot(data = m_type_pred3b, aes(x = Age, y = .epred)) + geom_violin(aes(color = Age, fill = Age), alpha = 0.4) +
  stat_summary(detseq_oi, inherit.aes = FALSE, mapping=aes(x = Age, y = n_peel, color = Age), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Age", y = "Average number of peels per sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

##### 4. Number of mistakes ######
# we coded three different "types" of mistakes: the item flying off, dropping the hammer, and general misstrikes (not hitting the item)
table(detseq_oi$n_miss, detseq_oi$Age)
table(detseq_oi$n_flies, detseq_oi$Age)
ftable(detseq_oi$n_hloss, detseq_oi$Age)
# in all cases mistakes are very rare, so these models are heavily zero-inflated
ftable(detseq_oi$n_misstotal)
sum(detseq_oi$n_misstotal)

# anvil type could play a big role here, with items being more likely to fly off on some material
t.test(detseq_oi$n_flies ~ as.factor(detseq_oi$anviltype))
# see more flying on stone than on wood

# do not pool the different mistakes together, since different processes might underlie them

## 4a: True misses (n_miss)
descdist(detseq_oi$n_miss)
testdist3.1 <- fitdist(detseq_oi$n_miss, "pois")
plot(testdist3.1)
# use a zero-inflated poisson

### Model_e4a ### 
# Outcome: number of misstrikes (true misses)
# Fixed effects: age, item, and anviltype
# Random effects: subjectID
m_e4a <- brm(n_miss ~ Age + item + anviltype + (1|subjectID), data = detseq_oi, family = zero_inflated_poisson, 
             iter = 2000, chain = 2, core = 2, save_pars = save_pars(all = TRUE), backend = "cmdstanr", control = list(adapt_delta = 0.99))
# m_e4a <- add_criterion(m_e4a, c("loo", "loo_R2", "bayes_R2"), reloo = TRUE, backend = "cmdstanr", ndraws = 2000) 

# saving and loading model
# saveRDS(m_e4a, "detailedtools/RDS/m_e4a.rds")
# m_e4a <- readRDS("detailedtools/RDS/m_e4a.rds")

# diagnostics
summary(m_e4a)
mcmc_plot(m_e4a)
pp_check(m_e4a)
plot(conditional_effects(m_e4a))

loo(m_e4a) # all cases good
loo_R2(m_e4a) # 0.17
round(bayes_R2(m_e4a),2) # 0.16

# Interpretation
round(exp(-0.77-3.63),2) * 0.48
hypothesis(m_e4a, "Intercept > Intercept + AgeSubadult", alpha = 0.05)

# make violin plot
m_type_pred4 <- m_e4a %>% 
  epred_draws(newdata = tibble(Age = detseq_oi$Age,
                               item = detseq_oi$item,
                               anviltype = detseq_oi$anviltype,
                               subjectID = detseq_oi$subjectID))

# age difference in misstrikes
ggplot(data = m_type_pred4, aes(x = Age, y = .epred)) + geom_boxplot(aes(color = Age, fill = Age), alpha = 0.4) +
  stat_summary(detseq_oi, inherit.aes = FALSE, mapping=aes(x = Age, y = n_miss, color = Age), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Age", y = "Average number of mistakes per tool use sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

## 4b: Item flying (itemflies)

### Model_e4b ### 
# Outcome: number of itemflies (item flying off)
# Fixed effects: age, item, and anviltype
# Random effects: subjectID
m_e4b <- brm(n_flies ~ Age + item + anviltype + (1|subjectID), data = detseq_oi, 
             family = zero_inflated_poisson, iter = 2000, chain = 2, core = 2, 
             save_pars = save_pars(all = TRUE),
             backend = "cmdstanr", control = list(adapt_delta = 0.99))
# m_e4b <- add_criterion(m_e4b, c("loo", "loo_R2", "bayes_R2"), reloo = TRUE, backend = "cmdstanr", ndraws = 2000) 

# saving and loading model
# saveRDS(m_e4b, "detailedtools/RDS/m_e4b.rds")
# m_e4b <- readRDS("detailedtools/RDS/m_e4b.rds")

# diagnostics
summary(m_e4b)
mcmc_plot(m_e4b)
pp_check(m_e4b)
plot(conditional_effects(m_e4b))

loo(m_e4b) # all cases good
loo_R2(m_e4b) # 0.10
round(bayes_R2(m_e4b),2) # 0.13

# Interpretation
round(exp(-0.54),2) * 0.40
hypothesis(m_e4b, "Intercept > Intercept + AgeAdult", alpha = 0.05)
hypothesis(m_e4b, "Intercept > Intercept + anviltypewood", alpha = 0.05)
hypothesis(m_e4b, "Intercept + itemalmendragreen > Intercept", alpha = 0.05)
hypothesis(m_e4b, "Intercept  + itemalmendrared > Intercept", alpha = 0.05)

# make violin plot
m_type_pred4b <- m_e4b %>% 
  epred_draws(newdata = tibble(Age = detseq_oi$Age,
                               item = detseq_oi$item,
                               anviltype = detseq_oi$anviltype,
                               subjectID = detseq_oi$subjectID))

# age difference in items flying
ggplot(data = m_type_pred4b, aes(x = Age, y = .epred)) + geom_boxplot(aes(color = Age, fill = Age), alpha = 0.4) +
  stat_summary(detseq_oi, inherit.aes = FALSE, mapping=aes(x = Age, y = n_flies, color = Age), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Age", y = "Average number of items flying per tool use sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

##### Visualizations ####
# combine dataframes
m_type_predtotal <- m_type_pred
m_type_predtotal$.epred_pound <- m_type_pred2$.epred
m_type_predtotal$.epred_reposition <- m_type_pred3$.epred
m_type_predtotal$.epred_peel <- m_type_pred3b$.epred
m_type_predtotal$.epred_misstrike <- m_type_pred4$.epred
m_type_predtotal$.epred_itemflies <- m_type_pred4b$.epred

## pounds, repositions and peels
m_type_predtotal1 <- melt(m_type_predtotal[,c("item", "Age", "anviltype", "subjectID",
                                         ".epred_pound", ".epred_reposition", ".epred_peel")], 
                          measure.vars = c(".epred_pound", ".epred_reposition", ".epred_peel"))

# colors for pound, reposit, peel
cols_a <- viridis(4, option = "viridis", end = 0.8)

# png("detailedtools/RDS/poundrepositpeel.png", width = 8, height = 7, units = 'in', res = 300)
ggplot(data = m_type_predtotal1, aes(x = Age)) + 
  geom_violin(aes(x = Age, y = value, col = variable, fill = variable), alpha = 0.4) +
  stat_summary(detseq_oi, inherit.aes = FALSE, 
               mapping=aes(x = as.numeric(Age) -0.3, y = n_pounds), col = cols_a[1], 
               geom = "point", fun = "mean", size = 3)+
  stat_summary(detseq_oi, inherit.aes = FALSE, mapping=aes(x = Age, y = n_itemreposit), 
               col = cols_a[2], geom = "point", fun = "mean",size = 3) +   
  stat_summary(detseq_oi, inherit.aes = FALSE, mapping=aes(x = as.numeric(Age) +0.3, y = n_peel), 
               col =  cols_a[3], geom = "point", fun = "mean", size = 3) +
  scale_fill_viridis_d(option = "viridis", end = 0.8, labels = c("pounds", "repositions", "peels")) +
  scale_color_viridis_d(option = "viridis", end = 0.8) +
  guides(fill = guide_legend("Type"), color = "none") +
  labs(x = "Age", y = "Average number per sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 
#dev.off()  

## mistrikes and item flies
m_type_predtotal2 <- melt(m_type_predtotal[,c("item", "Age", "anviltype", "subjectID",
                                              ".epred_misstrike", ".epred_itemflies")], 
                          measure.vars = c(".epred_misstrike", ".epred_itemflies"))

cols2 <- viridis(2, option = "inferno", end = 0.8)

# png("detailedtools/RDS/mistakes.png", width = 8, height = 7, units = 'in', res = 300)
ggplot(data = m_type_predtotal2, aes(x = Age)) + 
  geom_violin(aes(x = Age, y = value, col = variable, fill = variable), alpha = 0.4) +
  stat_summary(detseq_oi, inherit.aes = FALSE, 
               mapping=aes(x = as.numeric(Age) -0.25, y = n_miss), col = cols2[1], 
               geom = "point", fun = "mean", size = 3)+
  stat_summary(detseq_oi, inherit.aes = FALSE, mapping=aes(x = as.numeric(Age) +0.25, y = n_flies), 
               col = cols2[2], geom = "point", fun = "mean",size = 3) +   
  scale_fill_viridis_d(option = "inferno", end = 0.8, labels = c("misstrikes", "item flies")) +
  scale_color_viridis_d(option = "inferno", end = 0.8) +
  guides(fill = guide_legend("Type"), color = "none") +
  labs(x = "Age", y = "Average number per sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) + ylim(0,2)
#dev.off()

##
### INDIVIDUAL VARIATION & DEVELOPMENT ####
##

##### Individual variation ######
# filter to individuals with enough data
ftable(detseq_oi$subjectID)
knownids <- data.frame(ID = c("TER", "SPT", "BAL", "TOM", "PEA", "SMG", "ZIM", "MIC", "LAR"))
for (i in 1:nrow(knownids)) {
  knownids$nrow[i] <- nrow(detseq_o[which(detseq_o$subjectID == knownids$ID[i]),])
}

detseq_o2 <- left_join(detseq_o[which(detseq_o$subjectID %in% knownids$ID),], knownids, by = c("subjectID" = "ID"))

# make subjectID factor that is ordered based on age
detseq_o2$subjectID <- factor(detseq_o2$subjectID, levels = c("ZIM", "PEA", "BAL", "TER", "MIC", "LAR", "SPT", "TOM", "SMG"))

# number of repositions
ggplot(detseq_o2, aes(x=subjectID, y=n_reposit, color = Age, fill = Age)) + 
  geom_violin(alpha = 0.4) + geom_text(aes(y = 9.5, x = subjectID, label = nrow)) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  labs(x = "Age", y = "Average number of repositions per tool use sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

# individual variation in how many mistakes are made
ggplot(detseq_o2, aes(x=subjectID, y=n_flies, color = Age, fill = Age)) + 
  geom_violin(alpha = 0.4) + geom_text(aes(y = 7.5, x = subjectID, label = nrow)) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  labs(x = "Age", y = "Average number of item flying per tool use sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

ggplot(detseq_o2, aes(x=subjectID, y=n_hloss, color = Age, fill = Age)) + 
  geom_violin(alpha = 0.4) + geom_text(aes(y = 7.5, x = subjectID, label = nrow)) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  labs(x = "Age", y = "Average number of hammer loss per tool use sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

ggplot(detseq_o2, aes(x=subjectID, y=n_miss, color = Age, fill = Age)) + 
  geom_violin(alpha = 0.4) + geom_text(aes(y = 7.5, x = subjectID, label = nrow)) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  labs(x = "Age", y = "Average number of misstrikes per tool use sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

# comparing n_pounds, n_miss, n_reposit for known individuals
melt_detseq <- melt(detseq_o2, measure.vars = c("n_pounds", "n_misstotal", "n_itemreposit"))

ggplot(melt_detseq) + geom_violin(aes(y = value, x = variable, color = variable, fill = variable), alpha = 0.4) +
  stat_summary(melt_detseq, inherit.aes = FALSE, mapping=aes(x = variable, y = value, color = variable), geom = "point", fun = "mean",
               size = 4) +
  facet_wrap(~subjectID) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  labs(x = "Age", y = "Average number of repositions per sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

## What items they process over time
ggplot(detseq_o2, aes(x = month(mediadate), fill = item)) + geom_histogram() + facet_wrap(~ subjectID) + theme_bw()

##### Development ######

# First just visualizing how number of pounds/repositions/mistakes change over time

# how n_pounds changes over time, only for brown almendras
ggplot(detseq_o2[detseq_o2$item == c("almendrabrown"),]) + geom_point(aes(x = mediadate, y = n_pounds, shape = location), alpha = 0.4, color = "brown", size = 3) + geom_smooth(aes(x = mediadate, y = n_pounds), color = "brown") + 
  facet_wrap(~subjectID) + theme_bw() + ggtitle("Brown Almendra") + labs(x = "Date", y = "Number of pounds to open item") +theme(axis.text = element_text(size = 12),
                                                                                  axis.title = element_text(size = 14)) 
# n_pounds for green almendras
ggplot(detseq_o2[detseq_o2$item == c("almendragreen"),]) + geom_point(aes(x = mediadate, y = n_pounds, shape = location), alpha = 0.4, color = "darkgreen", size = 3) + geom_smooth(aes(x = mediadate, y = n_pounds), color = "darkgreen") + 
  facet_wrap(~subjectID) + theme_bw() + ggtitle("Green Almendra") + labs(x = "Date", y = "Number of pounds to open item") +theme(axis.text = element_text(size = 12),
                                                                                                                                 axis.title = element_text(size = 14)) 
# how misstrikes (real misses) change over time, brown almendras
ggplot(detseq_o2[detseq_o2$item == c("almendrabrown"),]) + geom_point(aes(x = mediadate, y = n_miss, shape = location), alpha = 0.4, color = "brown", size = 3) + geom_smooth(aes(x = mediadate, y = n_miss), color = "brown") + 
  facet_wrap(~subjectID) + theme_bw() + ggtitle("Brown Almendra") + labs(x = "Date", y = "Number of mistakes") +theme(axis.text = element_text(size = 12),
                                                                                                                                axis.title = element_text(size = 14)) + ylim(c(0,10))
# how itemrepositions change over time
ggplot(detseq_o2[detseq_o2$item == c("almendrabrown"),]) + geom_point(aes(x = mediadate, y = n_itemreposit, shape = location), alpha = 0.4, color = "brown", size = 3) + geom_smooth(aes(x = mediadate, y = n_itemreposit), color = "brown") + 
  facet_wrap(~subjectID) + theme_bw() + ggtitle("Brown Almendra") + labs(x = "Date", y = "Number of repositions per sequence") +theme(axis.text = element_text(size = 12),                                                                                                                  axis.title = element_text(size = 14)) + ylim(c(0,10))

# combining n_pounds, n_miss and n_repositions for brown almendras
# png("detailedtools/RDS/eff_dev.png", width = 9, height = 7, units = 'in', res = 300)
ggplot(detseq_o2[detseq_o2$item == c("almendrabrown"),]) + 
  geom_point(aes(x = mediadate, y = n_pounds, color = "Pounds", shape = location), alpha = 0.2, size = 2) + 
  geom_point(aes(x = mediadate, y = n_itemreposit, color = "Repositions", shape = location), alpha = 0.2, size = 2) + 
  geom_smooth(aes(x = mediadate, y = n_pounds, color = "Pounds")) +  
  geom_smooth(aes(x = mediadate, y = n_itemreposit, color = "Repositions")) + 
  geom_point(aes(x = mediadate, y = n_miss, color = "Misses", shape = location), alpha = 0.2, size = 3) + 
  geom_smooth(aes(x = mediadate, y = n_pounds, color = "Pounds")) +  
  geom_smooth(aes(x = mediadate, y = n_miss, color = "Misses")) + theme_bw() +
  theme(legend.position = "top", legend.title = element_blank()) +
  facet_wrap2(.~subjectID, nrow = 3, ncol = 3,
              strip = strip_themed(
                background_x = list(element_rect(fill = cols[1]),
                                    element_rect(fill = cols[1]),
                                    element_rect(fill = cols[1]),
                                    element_rect(fill = cols[1]),
                                    element_rect(fill = cols[2]),
                                    element_rect(fill = cols[2]),
                                    element_rect(fill = cols[2]),
                                    element_rect(fill = cols[3]),
                                    element_rect(fill = cols[3])))) + 
  theme(strip.text = element_text(colour = 'white')) +
  scale_color_manual("", breaks = c("Pounds", "Misses", "Repositions"), values = c(cols_a[1], cols_a[4], cols_a[2])) +
  scale_x_datetime(labels = date_format("%b")) +
  labs(x = "Date", y = "Number per sequence") +theme(axis.text = element_text(size = 12),  axis.title = element_text(size = 14)) 
# dev.off()

###### GAMs of tool use development ######
detseq_o2$time <- as.numeric(difftime(detseq_o2$mediadate, min(detseq_o2$mediadate), unit = "days"))
head(detseq_o2$time)
detseq_o2$subjectID_F <- as.factor(detseq_o2$subjectID)
detseq_o2$Age_f <- as.factor(detseq_o2$Age)

## Model dev_gam1 ##
# outcome: n_pounds
# fixed effects: time (days since start of deployment), with a curve for each age class
# random effect of subjectID

# MGCV
dev_gam1 <- gam(n_pounds ~ s(time, by = Age_f) + Age_f + s(subjectID_F, bs = "re"), 
                data = detseq_o2[detseq_o2$item == "almendrabrown",], family = "poisson", method = "REML")
summary(dev_gam1)
draw(dev_gam1)
plot(dev_gam1)

# BRMS
dev_gam1b <- brm(n_pounds ~ s(time, by = Age_f) + Age_f + s(subjectID_F, bs = "re"), 
                 data=detseq_o2[detseq_o2$item == "almendrabrown",], family="poisson", 
               chains=2, cores = 2, backend = "cmdstanr", save_pars = save_pars(all = TRUE),
               iter = 2000)
# dev_gam1b <- add_criterion(dev_gam1b, c("loo", "loo_R2", "bayes_R2"), reloo = TRUE, backend = "cmdstanr", ndraws = 2000) 

# saving and loading model
# saveRDS(dev_gam1b, "detailedtools/RDS/dev_gam1b.rds")
# dev_gam1b <- readRDS("detailedtools/RDS/dev_gam1b.rds")
plot(conditional_smooths(dev_gam1b))
plot(conditional_effects(dev_gam1b))

## BRMS linear
dev_m1 <- brm(n_pounds ~ time*subjectID_F, data=detseq_o2[detseq_o2$item == "almendrabrown",], family="poisson", 
              chains=2, cores = 2, backend = "cmdstanr", save_pars = save_pars(all = TRUE),
              iter = 2000)
# dev_m1 <- add_criterion(dev_m1, c("loo", "loo_R2", "bayes_R2"), reloo = TRUE, backend = "cmdstanr", ndraws = 2000) 

# saving and loading model
# saveRDS(dev_m1, "detailedtools/RDS/dev_m1.rds")
# dev_m1 <- readRDS("detailedtools/RDS/dev_m1.rds")
summary(dev_m1)

## Model dev_gam2 ##
# outcome: n_pounds
# fixed effects: time (days since start of deployment), with a curve for each individual

# MGCV
dev_gam2 <- gam(n_pounds ~ s(time, by = subjectID_F) + subjectID_F, 
                data = detseq_o2[detseq_o2$item == "almendrabrown",], family = "poisson", method = "REML")
summary(dev_gam2)
draw(dev_gam2)
plot(dev_gam2)

# BRMS
dev_gam2b <- brm(n_pounds ~ s(time, by = subjectID_F) + subjectID_F, data=detseq_o2[detseq_o2$item == "almendrabrown",], family="poisson", 
                 chains=2, cores = 2, backend = "cmdstanr", save_pars = save_pars(all = TRUE),
                 iter = 2000)
# dev_gam2b <- add_criterion(dev_gam2b, c("loo", "loo_R2", "bayes_R2"), reloo = TRUE, backend = "cmdstanr", ndraws = 2000) 

# saving and loading model
# saveRDS(dev_gam2b, "detailedtools/RDS/dev_gam2b.rds")
# dev_gam2b <- readRDS("detailedtools/RDS/dev_gam2b.rds")
plot(conditional_smooths(dev_gam2b))
plot(conditional_effects(dev_gam2b))

###### Old age cut-off? ######
# do we see ABE use tools at all in this dataset?
ftable(detseq$subjectID)
# we don't see ABE use tools in this dataset at all

# does ABE use tools in old data? 
ABE_only <- subset(agouticlean, agouticlean$name == "ABE (Abraham)")
ftable(ABE_only$tooluse) # 73 instances of tool use
plot(ABE_only$seq_start, ABE_only$tooluse)

# do we see ABE displace/scrounge? 
ABEcomment <- dettools_r2[str_detect(dettools_r2$comment, "ABE|Abraham"),]
# 13 recorded instances of him displacing and scrounging
# how many of displacements is ABE?
displacers <- append(socatt_seq$disp_ID1, socatt_seq$disp_ID2)
displacers <- displacers[!is.na(displacers)]
ftable(displacers)
ftable(socatt_final$n_disp) # 165 total displacements
ftable(socatt_final$n_scr) # 400 total scrounging events

###
### SOCIAL ATTENTION ####
###

# in how many sequences do we see social attention?
ftable(socatt_final$socialattention[!duplicated(socatt_final$sequenceID)])

# what items are being processed
ftable(socatt_final$item[!duplicated(socatt_final$sequenceID)])
# maybe pool together items than almendras? (because others are so rare) 
# or pool together
socatt_final$item2 <- ifelse(str_detect(socatt_final$item, "almendra") == FALSE, "other", socatt_final$item)
ftable(socatt_final$item2[!duplicated(socatt_final$sequenceID)])

ftable(socatt_final$tooluserID)
ftable(socatt_final$tooluser_age)
ftable(socatt_final$observerID)
ftable(socatt_final$observer_agesex, socatt_final$socatt)

# check for NAs
which(is.na(socatt_final), arr.ind = TRUE)
str(socatt_final)
# make factors
socatt_final$tooluser_age <- factor(socatt_final$tooluser_age, levels = c("Juvenile", "Subadult", "Adult"))
socatt_final$observer_agesex <- factor(socatt_final$observer_agesex, levels = c("juvenile", "subadultmale", "adultmale", "adultfemale"))
socatt_final$location <- as.factor(socatt_final$location)
socatt_final$tooluserID <- as.factor(socatt_final$tooluserID)
socatt_final$observerID <- as.factor(socatt_final$observerID)
socatt_final$item2 <- as.factor(socatt_final$item2)
socatt_final$sequenceID <- as.factor(socatt_final$sequenceID)

# get IDs of who pays social attention
socialattentioners <- append(socatt_seq$socatt_ID1, socatt_seq$socatt_ID2)
socialattentioners <- append(socialattentioners, socatt_seq$socatt_ID3)
socialattentioners <- socialattentioners[!is.na(socialattentioners)]
ftable(socialattentioners)

##### What predicts social attention? #####

### Model socatt_bm1 ###
# every line in the dataframe is a capuchin present during a tool use sequence
# Outcome: social attention yes/no 
# Fixed effects: age of the tool user, age-sex of the observer, total number of capuchins present, location, total number of scroungers
# Random effect: sequenceID (to account for repeated rows per sequence)
# Offset of log(sequence duration) reflecting opportunity (longer sequences mean more chances for social attention)
socatt_bm1 <- brm(socatt ~ tooluser_age + observer_agesex + location + p_total + n_scr + (1|sequenceID) + offset(log(seqduration)),
                  data = socatt_final, family = bernoulli(),
                  iter = 2000, chains=2, cores = 2, backend = "cmdstanr", save_pars = save_pars(all = TRUE))
# socatt_bm1 <- add_criterion(socatt_bm1, c("loo", "loo_R2", "bayes_R2"), reloo = TRUE, backend = "cmdstanr", ndraws = 2000) 

# saving and loading model
# saveRDS(socatt_bm1, "detailedtools/RDS/socatt_bm1.rds")
# socatt_bm1 <- readRDS("detailedtools/RDS/socatt_bm1.rds")

# Diagnostics
summary(socatt_bm1)
mcmc_plot(socatt_bm1)
plot(conditional_effects(socatt_bm1))
plot(socatt_bm1)

loo(socatt_bm1) # all cases good
loo_R2(socatt_bm1) # 0.16
round(bayes_R2(socatt_bm1),2) # 0.33

# Interpretation
ggpredict(socatt_bm1, term = c("tooluser_age"))
ggpredict(socatt_bm1, term = c("location", "observer_agesex[juvenile]"))

hypothesis(socatt_bm1, "Intercept + itemalmendrared > Intercept", alpha = 0.05)
hypothesis(socatt_bm1, "Intercept + observer_agesexsubadultmale < Intercept + observer_agesexjuvenile", alpha = 0.05)
hypothesis(socatt_bm1, "Intercept > Intercept + anviltypewood", alpha = 0.05)

# Visualization
socatt_bm1_pred <- socatt_bm1 %>% 
  epred_draws(newdata = tibble(tooluser_age = socatt_final$tooluser_age,
                               observer_agesex = socatt_final$observer_agesex,
                               location = socatt_final$location,
                               p_total = socatt_final$p_total,
                               n_scr = socatt_final$n_scr,
                               sequenceID = socatt_final$sequenceID,
                               seqduration = socatt_final$seqduration))

# who receives and pays social attention
cols_obs <- c("#D35171FF", "#FCA636FF", "#0D0887FF", "#B12A90FF" )
toolusertitles <- c(`Juvenile` = "Juvenile tool user",
                    `Subadult` = "Subadult tool user",
                    `Adult` = "Adult tool user")

# png("detailedtools/RDS/socatt_ages.png", width = 8, height = 6, units = 'in', res = 300)
ggplot(data = socatt_bm1_pred, aes(x = observer_agesex, y = .epred)) + geom_boxplot(aes(color = observer_agesex, fill = observer_agesex), alpha = 0.4,  outlier.color = NA) +
  stat_summary(socatt_final, inherit.aes = FALSE, mapping=aes(x = observer_agesex, y = socatt, color = observer_agesex), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_manual(values = cols_obs) +
  scale_color_manual(values = cols_obs) +
  guides(color = "none", fill = "none") +
  labs(x = "Age-sex class of observer", y = "Likelihood to pay social attention") + facet_wrap(~tooluser_age, labeller = as_labeller(toolusertitles)) +
  theme_bw() + theme(axis.text = element_text(size = 12), strip.text.x = element_text(size = 12),
                     axis.title = element_text(size = 14)) + theme(axis.text.x = element_text(angle =45, hjust = 1))
# dev.off()

# scrounging and total presence
soc_present <- plot(conditional_effects(socatt_bm1), plot = FALSE)[[4]]
soc_scr <- plot(conditional_effects(socatt_bm1), plot = FALSE)[[5]]

presentplot <- data.frame(p_total = c(1,2,3,4,5,6,7))
scrplot <- data.frame(n_scr = c(1,2,3,4,5))

scrplot$socatt <- aggregate(socatt_final[, c("socatt")], list(socatt_final$n_scr), mean)[2]$x
scrplot$sizescrounge <- aggregate(socatt_final[, c("socatt")], list(socatt_final$n_scr), length)[2]$x
presentplot$socatt <- aggregate(socatt_final[, c("socatt")], list(socatt_final$p_total), mean)[2]$x
presentplot$sizepresent <- aggregate(socatt_final[, c("socatt")], list(socatt_final$p_total), length)[2]$x

# png("detailedtools/RDS/socatt_scrpre.png", width = 8, height = 6, units = 'in', res = 300)
soc_present  + theme_bw() +
  stat_summary(data = presentplot, inherit.aes = FALSE, aes(x = p_total, y = socatt, fill = sizepresent), 
  geom = "point", fun = "mean", size = 3, shape = 24, alpha = 0.5) + scale_fill_viridis(limits = c(0,1000)) + 
  labs(x = "Individuals present", y = "Likelihood of social attention") + guides(fill = "none") +
  theme(axis.title = element_text(size = 14),  axis.text = element_text(size = 12)) +
  soc_scr + theme_bw() +   scale_fill_viridis(limits = c(0, 1000)) +
  stat_summary(data = scrplot, inherit.aes = FALSE, aes(x = n_scr, y = socatt, fill = sizescrounge), 
  geom = "point", fun = "mean", size = 3, shape = 24, alpha = 0.5) + labs(y = NULL, x = "Individuals scrounging", fill = "Sample size") +
  theme(legend.text =  element_text(size = 12), 
        legend.title = element_text(size =14),axis.title = element_text(size = 14), axis.text = element_text(size = 12))
#dev.off()  

##### Social attention to efficient tool users #####

# need to subset to known tool users, and where they successfully opened the item
# because if they didn't open it, we don't know how efficient the sequence was
socatt_finali <- droplevels.data.frame(socatt_final[!socatt_final$tooluserID %in% c("adultmale", "juvenileunknown", "subadultmale") &
                                socatt_final$outcome == "opened",])
nrow(socatt_finali)

### Model socatt_bm2 ###
# Outcome: social attention 1/0
# Fixed effects: identity of tool user, age of tool user, interaction with n_pounds, agesex of observer, n_mistakes(total)
# Offset of sequence duration 
socatt_bm1b <- brm(socatt ~ tooluser_age + n_pounds + observer_agesex + n_misstotal + offset(log(seqduration)) + (1|tooluserID), 
                   data = socatt_finali, family = bernoulli(),
                   iter = 2000, chains=2, cores = 4, backend = "cmdstanr", save_pars = save_pars(all = TRUE))
# socatt_bm1b <- add_criterion(socatt_bm1b, c("loo", "loo_R2", "bayes_R2"), reloo = TRUE, backend = "cmdstanr", ndraws = 2000) 

# saving and loading model
# saveRDS(socatt_bm1b, "detailedtools/RDS/socatt_bm1b.rds")
# socatt_bm1b <- readRDS("detailedtools/RDS/socatt_bm1b.rds")

# Diagnostics
summary(socatt_bm1b)
mcmc_plot(socatt_bm1b)
plot(conditional_effects(socatt_bm1b, re_formula = NULL))

loo(socatt_bm1b) # all cases good
loo_R2(socatt_bm1b) # 0.04
round(bayes_R2(socatt_bm1b),2) # 0.09

hypothesis(socatt_bm1b, "Intercept  > Intercept + n_pounds", alpha = 0.05)
ggpredict(socatt_bm1b, term = c("n_pounds"), bias_correction = TRUE)

# Visualization
socatt_bm1b_pred <- socatt_bm1b %>% 
  epred_draws(newdata = tibble(tooluser_age = socatt_finali$tooluser_age,
                               observer_agesex = socatt_finali$observer_agesex,
                               n_pounds = socatt_finali$n_pounds,
                               n_misstotal = socatt_finali$n_misstotal,
                               sequenceID = socatt_finali$sequenceID,
                               tooluserID = socatt_finali$tooluserID,
                               seqduration = socatt_finali$seqduration))

socatt_idgraph <- data.frame(tooluserID = sort(unique(socatt_finali$tooluserID)))
graphvalues_socatt <- aggregate(socatt_finali$socatt[socatt_finali$socatt == 1], list(socatt_finali$tooluserID[socatt_finali$socatt == 1]), length)
graphvalues_socatt <- rbind(graphvalues_socatt,c("ZIM",0), c("JOE", 0))
graphvalues_socatt <- graphvalues_socatt[order(graphvalues_socatt$Group.1),]
graphopp_socatt <- aggregate(socatt_finali$socatt[socatt_finali$socatt == 0], list(socatt_finali$tooluserID[socatt_finali$socatt == 0]), length)
graphopp_socatt <- graphopp_socatt[sort(graphopp_socatt$Group.1),]
socatt_idgraph$value <- graphvalues_socatt$x
socatt_idgraph$opportunity <- graphopp_socatt$x

# png("detailedtools/RDS/socatt_ids.png", width = 8, height = 7, units = 'in', res = 300)
ggplot(data = socatt_bm1b_pred, aes(x = tooluserID, y = .epred)) + geom_violin(aes(color = tooluser_age, fill = tooluser_age), alpha = 0.4) +
  stat_summary(socatt_finali, inherit.aes = FALSE, mapping=aes(x = tooluserID, y = socatt, color = tooluser_age), geom = "point", fun = "mean",
               size = 4) + scale_fill_manual(values = cols_obs) + scale_color_manual(values = cols_obs) +
  guides(color = "none") +   
    geom_text(socatt_idgraph, inherit.aes = FALSE, mapping=aes(x = tooluserID, y = 0.97, label=value), position=position_dodge(width=0.9), vjust=-0.25) +
  geom_text(socatt_idgraph, inherit.aes = FALSE, mapping=aes(x = tooluserID, y = 0.9, label = opportunity), position = position_dodge(width = 0), vjust = -0.5, fontface = "italic" ) +
  labs(x = "Identity of tool user", y = "Likelihood to receive social attention", fill = "Age") +
  theme_bw() + theme(axis.text = element_text(size = 12), strip.text.x = element_text(size = 12),
                     axis.title = element_text(size = 14)) + theme(axis.text.x = element_text(angle =45, hjust = 1))
# dev.off()

# consider doing a GAM with social attention depending on hour of the day?
# might be non-linear relationship?

