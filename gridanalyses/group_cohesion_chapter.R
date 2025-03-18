## Investigating group cohesion in grid data of fixed anvil tool-using vs non-tool-using groups 
## MPI-AB; Z Goldsborough

## Packages required (check if still required after cleaning)
library(stringr)
library(ggplot2)
library(mgcv)
library(gratia)
library(reshape2)
library(asnipe)
library(igraph)
library(activity)
library(lme4)
library(brms)
library(dplyr)
library(geodist)
library(sna)
library(assortnet)
library(ggnewscale)
library(bisonR)
library(hms)
library(tidyverse)
library(oneinfl)
library(tidybayes)
library(emmeans)
library(broom)
library(broom.mixed)
library(mapview)
library(sf)
library(fitdistrplus)


# STILL DO
# PRIOR SIMULATIONS AND SET MY PRIORS
# RUN FINAL MODELS WITH 3 CHAINS 3000 ITERATIONS AND ADD_CRITERION ETC



##### DIAGNOSTICS & FILTERING ####
# Load in grid data
# observation level dataframe
gridclean <- read.csv("gridanalyses/RDS/gridclean.csv")
ftable(gridclean$locationName)
# sequence level dataframe
gridsequence <- read.csv("gridanalyses/RDS/gridsequence.csv")
ftable(gridsequence$locationName)

## exclude grid cameras that are blank due to malfunctions
# NTU-151 and TU-168 have only blanks, TU-152 was pointed at the ground and therefore had mostly blanks from the start. Exclude these three.
gridclean_c <- gridclean[! gridclean$locationName %in% c("NTU-151", "TU-168", "TU-152"),]
# differentiate tool user and non tool user grid
gridclean_c$gridtype <- ifelse(str_detect(gridclean_c$locationName, "NTU") == TRUE, "NTU", "TU")
gridclean_c <- droplevels.data.frame(gridclean_c)

gridsequence_c <-  gridsequence[! gridsequence$locationName %in% c("NTU-151", "TU-168", "TU-152"),]
gridsequence_c$gridtype <- ifelse(str_detect(gridsequence_c$locationName, "NTU") == TRUE, "NTU", "TU")
gridsequence_c <- droplevels.data.frame(gridsequence_c)

## Did we have at least one capuchin detection at all remaining cameras?
ftable(gridsequence_c[which(gridsequence_c$capuchin == 1),]$locationName)
## any weird time issues?
hist(gridsequence_c[which(gridsequence_c$capuchin == 1),]$hour)

# filter down to only capuchin detections and minor cleaning
gridseq_oc <- gridsequence_c[gridsequence_c$capuchin == 1,]
gridseq_oc$gridtype <- as.factor(gridseq_oc$gridtype)
gridseq_oc$deplengthhours <- gridseq_oc$dep_length_hours
gridseq_oc <- droplevels.data.frame(gridseq_oc)

## Filter out detections of unfamiliar individuals
unfamiliars <-gridseq_oc$sequenceID[which(str_detect(gridseq_oc$observationComments, "unfamiliar|Unfamiliar") == TRUE)]
gridseq_ocf <- gridseq_oc[! gridseq_oc$sequenceID %in% unfamiliars,]
# how successful were we at assigning IDs and age sex
ftable(gridseq_ocf$agesex, gridseq_ocf$toolusers)
ftable(gridseq_ocf$individualID, gridseq_ocf$toolusers)

# visualize activity at the different cameras
gridcamerasmap <- as.data.frame(ftable(gridseq_oc$locationName))
colnames(gridcamerasmap) <- c("locationName", "ncapseq")
gridcamerasmap <- left_join(gridcamerasmap, gridseq_oc[!duplicated(gridseq_oc$locationName),c("locationName", "longitude", "latitude")], by = "locationName")
# locations with unfamiliar sightings
unfamiliarsloc <- gridseq_oc$locationName[which(str_detect(gridseq_oc$observationComments, "unfamiliar|Unfamiliar") == TRUE)]
gridcamerasmap$unfamiliar <- ifelse(gridcamerasmap$locationName %in% unfamiliarsloc, 1.2, 1)
gridcammap <- st_as_sf(gridcamerasmap, coords = c("longitude", "latitude"), crs = 4326)

mapview(gridcammap, zcol = "ncapseq", at = c(0, 50, 100, 200, 300, 400), legend = TRUE)

###### Exposure ####
## How many camera trapping days at TU vs NTU
griddays <- gridsequence_c
griddays$dayloc <- paste(griddays$locationfactor, griddays$seqday, sep = " ")
griddays2 <- griddays[!duplicated(griddays$dayloc),]

# make overview of deployments we have and their start and end days
gridlocations_t <- data.frame(uniqueloctag = unique(gridsequence_c$uniqueloctag)) 
gridlocations_t <- left_join(gridlocations_t, gridsequence_c[,c("uniqueloctag", "dep_start", "dep_end", "locationfactor", "gridtype")], by = "uniqueloctag")
gridlocations_t <- gridlocations_t[!duplicated(gridlocations_t$uniqueloctag),]
# take time off and keep just date variable
gridlocations_t$dep_startday <- as.Date(gridlocations_t$dep_start, tz = "America/Panama", "%Y-%m-%d")
gridlocations_t$dep_endday <- as.Date(gridlocations_t$dep_end, tz = "America/Panama", "%Y-%m-%d")
# calculate days in each deployment (round up)
gridlocations_t$dep_days <- ceiling(difftime(gridlocations_t$dep_end, gridlocations_t$dep_start, units = c("days")))
# number of rows in the griddays2 dataframe (so how many days we have)
for (i in 1:nrow(gridlocations_t)) {
  gridlocations_t$nrow[i] <- nrow(griddays2[griddays2$uniqueloctag == gridlocations_t$uniqueloctag[i],])
}

gridlocations_t2 <- aggregate(gridlocations_t$dep_days, list(locationfactor  = gridlocations_t$locationfactor, gridtype = gridlocations_t$gridtype), FUN = sum)

sum(gridlocations_t$dep_days[gridlocations_t$gridtype == "NTU"])
sum(gridlocations_t$dep_days[gridlocations_t$gridtype == "TU"])
## so comparable number of trapping days, but more in NTU than TU grid

# How many locations
ftable(gridlocations_t2$gridtype) # one more location in NTU grid
# Average number of trapping days per location
summary(as.numeric(gridlocations_t2$x[gridlocations_t2$gridtype == "NTU"]))
summary(as.numeric(gridlocations_t2$x[gridlocations_t2$gridtype == "TU"]))
# slightly longer deployments in NTU grid than TU grid, but no dramatic differences

####### Distance from placed cameras to planned GPS coordinates ####
# load in files with original GPS coordinates
head(gridsequence_c)

# real locations
TUgridcams <- gridseq_oc[!duplicated(gridseq_oc$locationfactor) & gridseq_oc$gridtype == "TU", c("locationfactor", "latitude", "longitude")]
TUgridcams <- TUgridcams[order(TUgridcams$locationfactor),]
NTUgridcams <- gridseq_oc[!duplicated(gridseq_oc$locationfactor) & gridseq_oc$gridtype == "NTU", c("locationfactor", "latitude", "longitude")]
NTUgridcams <- NTUgridcams[order(NTUgridcams$locationfactor),]

# planned locations
plannedgridcams <- read.csv("gridanalyses/plannedgridlocations.csv")
TUgridcams_planreal <- left_join(TUgridcams, plannedgridcams, by = "locationfactor")
TUgridcams_planreal <- TUgridcams_planreal %>%
  mutate(
    dist = geosphere::distHaversine(cbind(longitude, latitude), cbind(X, Y))
  )

TUgridcams_planreal

NTUgridcams_planreal <- left_join(NTUgridcams, plannedgridcams, by = "locationfactor")
NTUgridcams_planreal <- NTUgridcams_planreal %>%
  mutate(
    dist = geosphere::distHaversine(cbind(longitude, latitude), cbind(X, Y))
  )
NTUgridcams_planreal

summary(c(TUgridcams_planreal$dist, NTUgridcams_planreal$dist))

###### One NTU group? ####
# checking if we are capturing a single NTU group or multiple?
# max group size seen
NTUgridseq <- gridseq_ocf[gridseq_ocf$gridtype == "NTU",]
NTUgridseq[which(NTUgridseq$n == max(NTUgridseq$n)),]
# in max group size, see 4 adult females, 4 adult males, 5 juveniles (of which one infant)

# look at supposed group composition (max number of adult males and adult females seen in one sequence and how many we have IDed)
max(NTUgridseq$nAF) # max of 4 adult females (have identified 5)
max(NTUgridseq$nAM) # max of 4 adult males (have identified 5, potentially 6)
max(NTUgridseq$nJU) # max of 6 juveniles
max(NTUgridseq$nSM) # max of 2 subadult males ( have identified 2, maybe 3)

# co-occurrence of identifiable individuals (SNA network)
# make dataframe with individual information
gridagesex <- gridclean_c[,c("individualID","lifeStage", "sex", "gridtype")]
gridagesex <- gridagesex[! is.na(gridagesex$individualID) == TRUE & ! duplicated(gridagesex$individualID),]
NTUgridagesex <- gridagesex[gridagesex$gridtype == "NTU",]
# for now just add real names in manually, later use key file
NTUgridagesex$col <- ifelse(NTUgridagesex$sex == "male", "lightblue", "pink")
NTUgridagesex$col <- ifelse(NTUgridagesex$lifeStage == "adult", NTUgridagesex$col, "lightgreen")
NTUgridagesex <- NTUgridagesex[order(NTUgridagesex$individualID),]

# I think data format needs to be sequenceID/individualID
# go to only NTU grid data and only sequence ID and individual ID (when individual ID was known)
NTUassoc <- gridclean_c[gridclean_c$gridtype == "NTU" & ! is.na(gridclean_c$individualID) == TRUE, c("sequenceID", "individualID")]

# then go from long to wide?
NTUassoc_w <- dcast(NTUassoc, sequenceID ~ individualID, fun.aggregate = length )
NTUassoc_w2 <- NTUassoc_w
NTUassoc_w2[,2:15] <- as.numeric(unlist(NTUassoc_w2[,2:15]))
rownames(NTUassoc_w2) <- NTUassoc_w2$sequenceID
NTUassoc_w2 <- NTUassoc_w2[,-c(1,2)]
## now we have a dataframe with all associations (whenever individuals were seen together in the same sequence) in GBI (group by individual) format
# use this with asnipe package to get a network 
adj.m <- get_network(NTUassoc_w2, association_index = "SRI")
assoc.g <- graph_from_adjacency_matrix(adj.m, "undirected", weighted = T)
plot(assoc.g, vertex.color =NTUgridagesex$col, dge.width = E(assoc.g)$weight*100)

net_NTU <- graph.adjacency(adj.m, mode = "undirected", weighted = TRUE, diag = FALSE)
# png("gridanalyses/RDS/SNA_NTU.png", width = 8, height = 6, units = 'in', res = 300)
plot(net_NTU, vertex.color = NTUgridagesex$col, edge.width = E(assoc.g)$weight*100)
# dev.off()
coms_NTU <- fastgreedy.community(net_NTU) #identify communities
NTUgridagesex$COM <- membership(coms_NTU) #assign membership of communities
plot(net_NTU, vertex.color =NTUgridagesex$col, edge.with = 20*E(net_NTU)$weight^2, mark.groups = coms_NTU)

# largely appears to be one group, only Drop and Kai are unsure, but they were also only seen very rarely 

##### DAILY ACTIVITY PATTERN #####
# Visually, what time of day do we see activity of capuchins? 
# colors for two histograms in one
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")

### Tool users vs non tool users
histTU <- hist(gridseq_ocf$hour[gridseq_ocf$gridtype == "TU"], breaks = seq(from = 0, to = 24, by = 1), xlim = c(0, 24), freq = FALSE)
histNTU <- hist(gridseq_ocf$hour[gridseq_ocf$gridtype == "NTU"], breaks = seq(from = 0, to = 24, by = 1), xlim = c(0, 24), freq = FALSE)
plot(histNTU, col = c2, freq = FALSE, main = "Tool users (blue) vs non-tool users (red)", xlab = "Time of Day", ylab = "Proportion of sequences with capuchins", ylim = c(0, 0.12))
plot(histTU, col = c1, freq = FALSE, add = TRUE)

# in general TU group appears to be active more later in the afternoon
# all early morning and late evening activity is from the TU group
nightowls <- gridseq_ocf[gridseq_ocf$hour < 5 | gridseq_ocf$hour > 19,]
# all of these cameras seem to have the correct time set. So this means we truly have a capuchin detection at midnight and one at 4 AM!

## Activity analysis
# set to solar time (activity classified as being during day or during night)
gridseq_ocf$timestamp <- ymd_hms(gridseq_ocf$seq_start, tz = "America/Panama")

tmp <- solartime(gridseq_ocf$timestamp,
                 gridseq_ocf$latitude,
                 gridseq_ocf$longitude,
                 tz = -5,
                 format = "%Y-%m-%d %H:%M:%S")

gridseq_ocf$solar <- tmp$solar
gridseq_ocf$clock <- tmp$clock

plot(gridseq_ocf$solar, gridseq_ocf$clock)

# compare TU to NTU
# TU
act_m1 <- fitact(gridseq_ocf$solar[gridseq_ocf$gridtype == "TU"], sample = "model", reps = 1000) 
#saveRDS(act_m1, "gridanalyses/RDS/act_m1.RDS")
#act_m1 <- readRDS("gridanalyses/RDS/act_m1.RDS")
plot(act_m1)
act_m1@act[1] * 24
# this means they spend 0.38 * 24 = 9 hours per day active

# NTU
act_m2 <- fitact(gridseq_ocf$solar[gridseq_ocf$gridtype == "NTU"], sample = "model", reps = 1000) 
#saveRDS(act_m2, "gridanalyses/RDS/act_m2.RDS")
#act_m2 <- readRDS("gridanalyses/RDS/act_m2.RDS")
plot(act_m2)
act_m2@act[1] * 24
# this means they spend 0.38 * 24 = 10 hours per day active

# plot both together on same axis
# png("gridanalyses/RDS/dailyactivity.png", width = 8, height = 6, units = 'in', res = 300)
plot(act_m1, yunit="density", data="none", las=1, lwd=2,
     tline=list(lwd=3, col = "#C8800F"), # Thick line 
     cline=list(lty=3, col = "#C8800F")) # Supress confidence intervals

plot(act_m2, yunit="density", data="none", add=TRUE, 
     tline=list(lty = 5, col="#81A956", lwd=3),
     cline=list(lty=3, col="#81A956"),
     points(y = rep(0,6), x = nightowls$hour, pch = 19, col = "#C8800F", cex = 3))

legend("topright", c("TU", "NTU"), col=c("#C8800F","#81A956"), lty=c(1,5), lwd=2)
#dev.off()

# overlap between the two
compareCkern(act_m1, act_m2, reps = 100)
# 0.896448990 lot of overlap

##### PARTY SIZE ####

###### 1: Mean party size ####
mean(gridseq_ocf$n[gridseq_ocf$gridtype == "NTU"])
mean(gridseq_ocf$n[gridseq_ocf$gridtype == "TU"])
hist(gridseq_ocf$n[gridseq_ocf$gridtype == "NTU"])
hist(gridseq_ocf$n[gridseq_ocf$gridtype == "TU"])

ftable(gridseq_ocf$n) # 68.94 percent of all data is 1's 
max(gridseq_ocf$n[gridseq_ocf$gridtype == "NTU"])
max(gridseq_ocf$n[gridseq_ocf$gridtype == "TU"])

## Social party size rather than total party size
gridseq_ocf$partysize <- gridseq_ocf$n - 1
hist(gridseq_ocf$partysize)

sps_bm1a <-  brm(bf(partysize ~ gridtype + (1|locationfactor) + offset(log(deplengthhours)), hu ~ gridtype + (1|locationfactor)), 
                 data = gridseq_ocf, family = hurdle_poisson(), iter = 2000, chain = 2, core = 2, backend = "cmdstanr")
#saveRDS(sps_bm1a, "gridanalyses/RDS/sps_bm1a.rds")
#sps_bm1a <- readRDS("gridanalyses/RDS/sps_bm1a.rds")

summary(sps_bm1a)
plot(sps_bm1a)
plot(conditional_effects(sps_bm1a))
pp_check(sps_bm1a)
hypothesis(sps_bm1a, "Intercept > Intercept + gridtypeTU")

# Interpretation
# https://www.andrewheiss.com/blog/2022/05/09/hurdle-lognormal-gaussian-brms/
tidy(sps_bm1a)
# probability of a 0 (party of 1) in NTU
hurdle_intercept <- tidy(sps_bm1a) |> 
  filter(term == "hu_(Intercept)") |> 
  pull(estimate)
plogis(hurdle_intercept)
# probability of a 0 (party of 1) in TU
hurdle_TU <- tidy(sps_bm1a) |>
  filter(term == "hu_gridtypeTU") |>
  pull(estimate)

(plogis(hurdle_intercept + hurdle_TU) - plogis(hurdle_intercept)) * 100
# tool using group increases probability of a 0 (party of 1) by 2.03 percentage points, on average

# conditional effects of 0/not 0 
plot(conditional_effects(sps_bm1a, dpar = "hu"))
plot(conditional_effects(sps_bm1a, dpar = "mu"))
# okay so see here clearly that at TU we are more likely to see a party of 1, and at NTU we see larger parties if it is not party of 1
# and all together also slightly more likely to have larger party at NTU

# using emmeans to extract estimates on real scale
# only non-zero mu part
emmeans(sps_bm1a, "gridtype", dpar = "hu", regrid = "response")
emmeans(sps_bm1a, "gridtype", dpar = "mu", regrid = "response")
emmeans(sps_bm1a, "gridtype", regrid = "response")


###### 2: Variability in party size ####
ggplot(data = gridseq_ocf, aes( x= hour, y = n, col = gridtype, shape = gridtype)) + geom_point(alpha = 0.5, aes(shape = gridtype)) + geom_smooth()
ggplot(data = gridseq_ocf, aes( x= seqday, y = n, col = gridtype)) + geom_point(alpha = 0.5,) + geom_smooth() + facet_wrap(~gridtype)

ggplot(data = gridseq_ocf[gridseq_ocf$gridtype == "TU",], aes( x= hour, y = n, col = day(seqday))) + geom_point(alpha = 0.5) + facet_wrap(~month) + scale_colour_viridis_c() 
ggplot(data = gridseq_ocf[gridseq_ocf$gridtype == "NTU",], aes( x= hour, y = n, col = day(seqday))) + geom_point(alpha = 0.5) + facet_wrap(~month) + scale_colour_viridis_c() 

head(gridseq_ocf)

# per day per camera have sd deviation of party size
gridseq_daysd <- do.call(data.frame, aggregate(gridseq_ocf$n, by = list(gridseq_ocf$gridtype, gridseq_ocf$locationName, gridseq_ocf$seqday), 
                                               FUN = function(x) c(sd = sd(x), n = length(x), mn = mean(x))))
colnames(gridseq_daysd) <- c("gridtype", "locationfactor", "seqday", "party_sd", "party_n", "party_mean")
head(gridseq_daysd)
# exclude all the NAs (which is when one party was observed)
gridseq_daysd <- gridseq_daysd[is.na(gridseq_daysd$party_sd) == FALSE,]
hist(gridseq_daysd$party_sd)
ftable(gridseq_daysd$party_sd[gridseq_daysd$party_n == 2])
## determine distribution
descdist(gridseq_daysd$party_sd)

testdist1.1 <- fitdist(gridseq_daysd$party_sd, "lognormal")
plot(testdist1.1)

# compare sd in party size between NTU and TU. locationfactor as random effect. number of parties as offset (?) or fixed effect
ps_bm1b <- brm(party_sd ~ gridtype + (1|locationfactor) + offset(log(party_n)), data = gridseq_daysd, 
               family = hurdle_gamma(), iter= 2000, chain =2, core = 2, backend = "cmdstanr")
#saveRDS(ps_bm1b), "gridanalyses/RDS/ps_bm1b.rds")
#ps_bm1b <- readRDS("gridanalyses/RDS/ps_bm1b.rds")
summary(ps_bm1b)
plot(conditional_effects(ps_bm1b))
hypothesis(ps_bm1b, "Intercept  > Intercept + gridtypeTU", alpha = 0.05)

##### PARTY COMPOSITION ####

###### 1: Adult males and adult females ####

# try zero-inflated poisson number of adult  females in a sequence (in dataframe only capuchin detections. TU vs NTU)
# same for adult males
# then combine the likelihood functions into one model? 
pc_bm1 <- brm(nAF ~ gridtype + offset(log(deplengthhours)) + (1|locationfactor), data = gridseq_ocf, family = zero_inflated_poisson(link = "log", link_zi = "logit"), iter = 2000, chain = 2, core = 2, backend = "cmdstanr")
#saveRDS(pc_bm1, "gridanalyses/RDS/pc_bm1.rds")
#pc_bm1 <- readRDS("gridanalyses/RDS/pc_bm1.rds")
summary(pc_bm1)
plot(conditional_effects(pc_bm1))
pp_check(pc_bm1)
hypothesis(pc_bm1, "Intercept > Intercept + gridtypeTU")

# significantly higher number of females per sequence in NTU than TU grid

#males
pc_bm2 <- brm(nAM ~ gridtype + offset(log(deplengthhours)) + (1|locationfactor), data = gridseq_ocf, family = zero_inflated_poisson(link = "log", link_zi = "logit"), iter = 2000, chain = 2, core = 2, backend = "cmdstanr")
#saveRDS(pc_bm2, "gridanalyses/RDS/pc_bm2.rds")
#pc_bm2 <- readRDS("gridanalyses/RDS/pc_bm2.rds")
plot(conditional_effects(pc_bm2))
pp_check(pc_bm2)

# also higher number of males per sequence in NTU than TU grid

#combining?
pc_bm3 <- brm(nAF ~ gridtype*nAM + offset(log(deplengthhours)) + (1|locationfactor), data = gridseq_ocf, family = zero_inflated_poisson(link = "log", link_zi = "logit"), iter = 2000, chain = 2, core = 2, backend = "cmdstanr")
#saveRDS(pc_bm3, "gridanalyses/RDS/pc_bm3.rds")
#pc_bm3 <- readRDS("gridanalyses/RDS/pc_bm3.rds")
plot(conditional_effects(pc_bm3))
pp_check(pc_bm3)
mcmc_plot(pc_bm3)
summary(pc_bm3)
hypothesis(pc_bm3, "Intercept > Intercept + gridtypeTU")

# plot with real data plotted over it
partycomp <- plot(conditional_effects(pc_bm3), plot = FALSE)[[3]]

pcompplot <- gridseq_ocf %>% 
  group_by(nAM, gridtype) %>%
  summarize_at(vars("nAF"), list(mean = mean, sd = sd, nsample = length))

pcompplot$se <- pcompplot$sd/sqrt(pcompplot$nsample)

# with real means on
ggplot() +  scale_color_manual(values = c("#81A956", "#C8800F")) + scale_fill_manual(values = c("#81A956", "#C8800F")) +
  geom_point(data = pcompplot, aes(x = nAM, y = mean, color = gridtype, group = gridtype, alpha = nsample), size = 3, inherit.aes = FALSE) + 
  scale_alpha_continuous(range = c(0.3,1)) +
  geom_line(data = partycomp$data, aes(x = nAM, y = estimate__, color = gridtype, group = gridtype), size = 1.5) +
  geom_ribbon(data = partycomp$data, aes(x = nAM, ymin = lower__, ymax = upper__, group = gridtype, fill = gridtype), alpha = 0.2) +
  labs(x = "Number of adult males", y = "Number of adult females") +  theme_bw() + 
  theme(strip.text.x = element_text(size = 16), axis.title = element_text(size = 16), legend.text =  element_text(size = 14), legend.title = element_text(size =14),
        axis.text = element_text(size = 12)) 

# with means and error bars
ggplot() +  scale_color_manual(values = c("#81A956", "#C8800F")) + scale_fill_manual(values = c("#81A956", "#C8800F")) +
  geom_point(data = pcompplot, aes(x = nAM, y = mean, color = gridtype),  size = 3, inherit.aes = FALSE, alpha = 0.5) +  
  geom_errorbar(data = pcompplot, aes(x = nAM, ymin = mean - se, 
                                      ymax =  mean + se, color= gridtype),
                width=.2, linewidth = 1, alpha = 0.5) +
  geom_line(data = partycomp$data, aes(x = nAM, y = estimate__, color = gridtype, group = gridtype), size = 1.5) +
  geom_ribbon(data = partycomp$data, aes(x = nAM, ymin = lower__, ymax = upper__, group = gridtype, fill = gridtype), alpha = 0.2) +
  labs(x = "Number of adult males", y = "Number of adult females") +  theme_bw() + 
  theme(strip.text.x = element_text(size = 16), axis.title = element_text(size = 16), legend.text =  element_text(size = 14), legend.title = element_text(size =14),
        axis.text = element_text(size = 12)) 

###### 2: Adults and juveniles ####
pc_bm4 <- brm(nJuvenile ~ gridtype*nAdults + offset(log(deplengthhours)) + (1|locationfactor), data = gridseq_ocf, family = zero_inflated_poisson(link = "log", link_zi = "logit"), iter = 2000, chain = 2, core = 2, backend = "cmdstanr")


###### SPATIAL COHESION #####

##### CO-OCCURRENCES ####

## Step 1: Generate distance matrix showing distance between each camera per grid
TUdistmat <- geodist::geodist(TUgridcams)
rownames(TUdistmat) <- TUgridcams$locationfactor
colnames(TUdistmat) <- TUgridcams$locationfactor
TUdistmat

NTUdistmat <- geodist::geodist(NTUgridcams)
rownames(NTUdistmat) <- NTUgridcams$locationfactor
colnames(NTUdistmat) <- NTUgridcams$locationfactor
NTUdistmat

### Using all TU data (not just grid), looking for co-occurrences within 60 seconds >150 m away
# exclude sequences with unfamiliar individuals (from grid) and CEBUS-03 (duplicate)
agoutiseq_jt <- agoutisequence_c[agoutisequence_c$capuchin == 1 & agoutisequence_c$island == "Jicaron" & agoutisequence_c$tool_site == 1 & !agoutisequence_c$locationfactor == "CEBUS-03",]
unfamiliars <- agouticlean$sequenceID[which(str_detect(agouticlean$observationComments, "unfamiliar") == TRUE & agouticlean$tool_site == 1 & agouticlean$island == "Jicaron")]
agoutiseq_jt <- agoutiseq_jt[! agoutiseq_jt$sequenceID %in% unfamiliars,]

# make distance matrix for all cameras
TUcams <- agoutiseq_jt[!duplicated(agoutiseq_jt$locationfactor), c("locationfactor", "latitude", "longitude")]
TUcams <- TUcams[order(TUcams$locationfactor),]
TUdistmat_all <- geodist::geodist(TUcams)
rownames(TUdistmat_all) <- TUcams$locationfactor
colnames(TUdistmat_all) <- TUcams$locationfactor

# make blank co-occurrence info data frame? 
agoutiseq_jt$cooccurrence <- 0
agoutiseq_jt$cooc_ID <- NA

cooccurrences <- data.frame(cooc_ID = "seqid", seqstart = NA, seqday = NA, cam1 = NA, cam2 = NA, cam3 =NA, distcam12 = 0, distcam13 = 0, nrseq = 0, nrcap_1 = 0, nrcap_2 = 0, nrcap_3 = 0,
                            nAdult_1 = 0, nAdult_2 = 0, nAdult_3 = 0, nSubadult_1 = 0, nSubadult_2 = 0, nSubadult_3 = 0, nJuvenile_1 = 0, nJuvenile_2 = 0,
                            nJuvenile_3 = 0, nUU_1 = 0, nUU_2 = 0, nUU_3 = 0, tooluse_1 = NA, tooluse_2 = NA, tooluse_3 = NA)

## taking 150 m < 60 seconds
for (i in 1:nrow(agoutiseq_jt)) {
  ## at beginning have some kind of check if the sequenceID is already in the co-occurence dataframe, if so can skip everything
  if(sum(str_detect(cooccurrences$cooc_ID, paste(agoutiseq_jt$sequenceID[i]))) == 0) {
    dist <- as.data.frame(subset(TUdistmat_all, rownames(TUdistmat_all) %in% agoutiseq_jt$locationfactor[i])) 
    cand_locs <- colnames(dist[,dist > 150]) # make list of candidate locations for co-occurrence (>150 m away)
    # filter to sequence that are at candidate location and on same day as sequence we're looking at 
    cand_seq <- agoutiseq_jt[agoutiseq_jt$locationfactor %in% cand_locs & agoutiseq_jt$seqday == agoutiseq_jt$seqday[i], c("sequenceID", "locationfactor", "seqday", "seq_start", "seq_end", "n", "nAdult", "nJuvenile","nSubadult", "nUU", "tooluse")]
    dist_m <- melt(dist)
    cand_seq$locationfactor <- as.character(cand_seq$locationfactor)
    dist_m$variable <- as.character(dist_m$variable)
    cand_seq <- left_join(cand_seq, dist_m, by = c("locationfactor" = "variable"))
    # see if there are any co-occurrences
    # if there is anything, then extract information from those sequences, both add to agoutiseq_jt dataframe, and to co-occurrence dataframe?
    if(nrow(cand_seq) > 0) {
      if(min(abs(difftime(agoutiseq_jt$seq_start[i], cand_seq$seq_start, unit = "s"))) < 60) {
        cand_seq$dtime <- difftime(agoutiseq_jt$seq_start[i], cand_seq$seq_start, unit = "s")
        cand_seq_t <- cand_seq[abs(cand_seq$dtime) < 60,]
        cand_seq_t <- cand_seq_t[!duplicated(cand_seq_t$locationfactor),]
        if(nrow(cand_seq_t) > 0) {
          agoutiseq_jt$cooccurrence <- 1
          agoutiseq_jt$cooc_ID[i] <- ifelse(nrow(cand_seq_t) == 1, paste(agoutiseq_jt$sequenceID[i], cand_seq_t$sequenceID[1], sep = ","),
                                            paste(agoutiseq_jt$sequenceID[i], cand_seq_t$sequenceID[1], cand_seq_t$sequenceID[2], sep = ","))
          cooccurrences[nrow(cooccurrences) +1,] <- c(agoutiseq_jt$cooc_ID[i], paste(agoutiseq_jt$seq_start[i]), paste(agoutiseq_jt$seqday[i]), paste(agoutiseq_jt$locationfactor[i]), 
                                                      paste(cand_seq_t$locationfactor[1]), paste(cand_seq_t$locationfactor[2]), cand_seq_t$value[1], cand_seq_t$value[2], nrow(cand_seq_t), agoutiseq_jt$n[i], 
                                                      cand_seq_t$n[1], cand_seq_t$n[2], agoutiseq_jt$nAdult[i], cand_seq_t$nAdult[1], cand_seq_t$nAdult[2],agoutiseq_jt$nSubadult[i], cand_seq_t$nSubadult[1], 
                                                      cand_seq_t$nSubadult[2], agoutiseq_jt$nJuvenile[i], cand_seq_t$nJuvenile[1], cand_seq_t$nJuvenile[2], agoutiseq_jt$nUU[i], cand_seq_t$nUU[1], 
                                                      cand_seq_t$nUU[2], paste(agoutiseq_jt$tooluse[i]), paste(cand_seq_t$tooluse[1]), paste(cand_seq_t$tooluse[2]))
        }
      }
    }
  }
  print(i)
}

# Using lenient criterion (the further, the more time)
# make curve for co-occurrences
# so 150 meters < 60 seconds, but at 300 meters < 120 seconds (or different relationship)
disttime <- data.frame(distance = seq(150, 1020, by = 10), time = 60)
str(disttime)
disttime$time[which(disttime$distance > 150)] <- 60 + ((disttime$distance[which(disttime$distance>150)]-150) * 0.4)

for (i in 1:nrow(agoutiseq_jt)) {
  ## at beginning have some kind of check if the sequenceID is already in the co-occurence dataframe, if so can skip everything
  if(sum(str_detect(cooccurrences$cooc_ID, paste(agoutiseq_jt$sequenceID[i]))) == 0) {
    dist <- as.data.frame(subset(TUdistmat_all, rownames(TUdistmat_all) %in% agoutiseq_jt$locationfactor[i])) 
    cand_locs <- colnames(dist[,dist > 150]) # make list of candidate locations for co-occurrence (>150 m away)
    # filter to sequence that are at candidate location and on same day as sequence we're looking at 
    cand_seq <- agoutiseq_jt[agoutiseq_jt$locationfactor %in% cand_locs & agoutiseq_jt$seqday == agoutiseq_jt$seqday[i], c("sequenceID", "locationfactor", "seqday", "seq_start", "seq_end", "n", "nAdult", "nJuvenile","nSubadult", "nUU", "tooluse")]
    dist_m <- melt(dist)
    cand_seq$locationfactor <- as.character(cand_seq$locationfactor)
    dist_m$variable <- as.character(dist_m$variable)
    cand_seq <- left_join(cand_seq, dist_m, by = c("locationfactor" = "variable"))
    # see if there are any co-occurrences
    # if there is anything, then extract information from those sequences, both add to agoutiseq_jt dataframe, and to co-occurrence dataframe?
    if(nrow(cand_seq) > 0) {
      # identify the difftime cutoff for every row
      cand_seq$cutoff <- disttime$time[findInterval(cand_seq$value, disttime$distance)]
      cand_seq$flag <- ifelse(abs(difftime(agoutiseq_jt$seq_start[i], cand_seq$seq_start, unit = "s")) < cand_seq$cutoff, 1, 0)
      if(sum(cand_seq$flag) > 0)  {
        cand_seq$dtime <- difftime(agoutiseq_jt$seq_start[i], cand_seq$seq_start, unit = "s")
        cand_seq_t <- cand_seq[cand_seq$flag == 1,]
        cand_seq_t <- cand_seq_t[!duplicated(cand_seq_t$locationfactor),]
        if(nrow(cand_seq_t) > 0) {
          agoutiseq_jt$cooccurrence <- 1
          agoutiseq_jt$cooc_ID[i] <- ifelse(nrow(cand_seq_t) == 1, paste(agoutiseq_jt$sequenceID[i], cand_seq_t$sequenceID[1], sep = ","),
                                            paste(agoutiseq_jt$sequenceID[i], cand_seq_t$sequenceID[1], cand_seq_t$sequenceID[2], sep = ","))
          cooccurrences[nrow(cooccurrences) +1,] <- c(agoutiseq_jt$cooc_ID[i], paste(agoutiseq_jt$seq_start[i]), paste(agoutiseq_jt$seqday[i]), paste(agoutiseq_jt$locationfactor[i]), 
                                                      paste(cand_seq_t$locationfactor[1]), paste(cand_seq_t$locationfactor[2]), cand_seq_t$value[1], cand_seq_t$value[2], nrow(cand_seq_t), agoutiseq_jt$n[i], 
                                                      cand_seq_t$n[1], cand_seq_t$n[2], agoutiseq_jt$nAdult[i], cand_seq_t$nAdult[1], cand_seq_t$nAdult[2],agoutiseq_jt$nSubadult[i], cand_seq_t$nSubadult[1], 
                                                      cand_seq_t$nSubadult[2], agoutiseq_jt$nJuvenile[i], cand_seq_t$nJuvenile[1], cand_seq_t$nJuvenile[2], agoutiseq_jt$nUU[i], cand_seq_t$nUU[1], 
                                                      cand_seq_t$nUU[2], paste(agoutiseq_jt$tooluse[i]), paste(cand_seq_t$tooluse[1]), paste(cand_seq_t$tooluse[2]))
        }
      }
    }
  }
  print(i)
}


cooccurrences <- cooccurrences[-1,]
cooccurrences[,7:22] <- as.numeric(unlist(cooccurrences[,7:22]))

# > 150 m < 60 seconds criterion
#saveRDS(cooccurrences, "gridanalyses/RDS/cooccurrences_fixedgps150.RDS")
#cooccurrences <- readRDS("gridanalyses/RDS/cooccurrences_fixedgps150.RDS")

# flexible criterion depending on distance
#saveRDS(cooccurrences, "gridanalyses/RDS/cooccurrences_fixedgpscurve.RDS")
#cooccurrences <- readRDS("gridanalyses/RDS/cooccurrences_fixedgpscurve.RDS")

cooccurrences_old <- readRDS("gridanalyses/RDS/cooccurrences_wronggps150.RDS")
cooccurrences_new <- readRDS("gridanalyses/RDS/cooccurrences_fixedgps150.RDS")
missed <- cooccurrences[-which(cooccurrences_new$cooc_ID %in% cooccurrences_old$cooc_ID),]
# so the reason we lost so many co occurrences when fixing the long/lat is because CEBUS-04 had a totally wrong location, so many cooccurrences where CEBUs-04 with CEBUS-01 or 02 (which are actually next to each other)
# SO THE NEW CO-OCCURRENCES ARE CORRECT! 

## Analyses

sum(cooccurrences$tooluse_1 == TRUE | cooccurrences$tooluse_2 == TRUE)
summary(cooccurrences$nrcap_1)
summary(cooccurrences$nrcap_2)
allcapnrs <- c(cooccurrences$nrcap_1, cooccurrences$nrcap_2)
sum(cooccurrences$nrcap_1 == 1 | cooccurrences$nrcap_2 == 1)
cooccurrences$ratio <- cooccurrences$nrcap_1/(cooccurrences$nrcap_1 + cooccurrences$nrcap_2)
summary(cooccurrences$ratio)

cooccurrences$tooluse_any <- factor(ifelse(cooccurrences$tooluse_1 == "TRUE" | cooccurrences$tooluse_2 == "TRUE" | cooccurrences$tooluse_3 == "TRUE", "TRUE", "FALSE"))

hist(hour(cooccurrences$seqstart))
ggplot(data = cooccurrences, aes(x = hour(cooccurrences$seqstart), fill = tooluse_any)) + geom_histogram(bins = 12, col = "black") + facet_wrap(~tooluse_any) + scale_x_continuous(breaks = 6:18) +
  theme_bw() + labs(fill = "Tool use", x = "Hour", y = "Number of co-occurrences") + theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), legend.title = element_text(size = 14),
                                                                                           legend.text = element_text(size = 12))

ggplot(data = cooccurrences, aes(x = hour(cooccurrences$seqstart))) + geom_histogram(bins = 12, col = "black", fill = "#C8800F") + scale_x_continuous(breaks = 6:18) +
  theme_bw() + labs(fill = "Tool use", x = "Hour", y = "Number of co-occurrences") + theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), legend.title = element_text(size = 14),
                                                                                           legend.text = element_text(size = 12))

max(unique(hour(cooccurrences$seqstart)))
sum(cooccurrences$nJuvenile_1 > 0 | cooccurrences$nJuvenile_2 > 0)
sum(cooccurrences$nAdult_1 > 0 | cooccurrences$nAdult_2 > 0)
sum(cooccurrences$nSubadult_1 > 0 | cooccurrences$nSubadult_2 > 0)

ggplot(data = cooccurrences, aes(x = month(cooccurrences$seqstart), fill = as.factor(year(cooccurrences$seqstart)))) + geom_histogram(bins = 12, col = "black")  + theme_bw() + scale_x_continuous(breaks = 0:12) +
  labs(fill = "Year") + xlab("Month")

# co-occurrences depending on distance
hist(cooccurrences$distcam12)
blub2 <- aggregate(cooccurrences$distcam12, by = list(cam1 = cooccurrences$cam1, cam2 = cooccurrences$cam2, distcam = cooccurrences$distcam12), FUN = ftable)


summary(cooccurrences$distcam12)
agoutiseq_jt$year <- year(agoutiseq_jt$seq_start)
length(unique(agoutiseq_jt$locationfactor[which(agoutiseq_jt$year == "2022")]))

str(cooccurrences)


cooccurrences$nrcap_lowest <- cooccurrences$nrcap_1
cooccurrences$nrcap_lowest[cooccurrences$nrcap_2 <= cooccurrences$nrcap_1] <- cooccurrences$nrcap_2[cooccurrences$nrcap_2 <= cooccurrences$nrcap_1] 
cooccurrences$nrcap_highest <- cooccurrences$nrcap_2
cooccurrences$nrcap_highest[cooccurrences$nrcap_2 <= cooccurrences$nrcap_1] <- cooccurrences$nrcap_1[cooccurrences$nrcap_2 <= cooccurrences$nrcap_1]

cooccurrences$tooluseoption_1 <- ifelse(str_detect(cooccurrences$cam1, "CEBUS") == TRUE, TRUE, FALSE)
cooccurrences$tooluseoption_2 <- ifelse(str_detect(cooccurrences$cam2, "CEBUS") == TRUE, TRUE, FALSE)

cooccurrences$truetooluse_1 <- cooccurrences$tooluse_1
cooccurrences$truetooluse_1[which(cooccurrences$tooluseoption_1 == FALSE)] <- NA
cooccurrences$truetooluse_2 <- cooccurrences$tooluse_2
cooccurrences$truetooluse_2[which(cooccurrences$tooluseoption_2 == FALSE)] <- NA

cooccurrences$tooluse <- ifelse(cooccurrences$truetooluse_1 == TRUE & cooccurrences$truetooluse_2 == TRUE, "both", 
                                ifelse(cooccurrences$truetooluse_1 == FALSE & cooccurrences$truetooluse_2 == FALSE, "neither", 
                                       ifelse(cooccurrences$truetooluse_1 == TRUE & cooccurrences$truetooluse_2 == FALSE | cooccurrences$truetooluse_1 == FALSE & cooccurrences$truetooluse_2 == TRUE, "one", NA)))

table(cooccurrences$truetooluse_1, cooccurrences$tooluseoption_1)

library(scales)
ggplot(data = cooccurrences, aes(x = nrcap_lowest, y = nrcap_highest)) + geom_bin2d(binwidth = c(1,1)) + theme_bw() +
  stat_bin2d(geom = "text", aes(label = ..count..), binwidth = 1, col = "white") + scale_fill_gradient(limits = c(0, 75), oob = squish) + 
  theme(axis.text.x = element_text(hjust = -10, size = 14), axis.text.y = element_text(vjust = -1, size = 14), 
        axis.ticks = element_blank(), legend.position = "none", axis.title = element_text(size = 18)) + labs(x = "Number of capuchins in smallest party", y = "Number of capuchins in largest party")

### co-occurrences only grid
# widening timeframe, 2 min instead of 1
disttime2 <- data.frame(distance = seq(150, 1020, by = 10), time = 120)
str(disttime2)
disttime2$time[which(disttime2$distance > 150)] <- 120 + ((disttime2$distance[which(disttime2$distance>150)]-150) * 0.8)

## TU 
gridseq_ocTU <- gridseq_ocf[gridseq_ocf$gridtype == "TU",]
# distance matrix is already made
TUdistmat

# make blank co-occurrence info data frame? 
gridseq_ocTU$cooccurrence <- 0
gridseq_ocTU$cooc_ID <- NA

cooccurrences_TU_TU <- data.frame(cooc_ID = "seqid", seqstart = NA, seqday = NA, cam1 = NA, cam2 = NA, cam3 =NA, distcam12 = 0, distcam13 = 0, nrseq = 0, nrcap_1 = 0, nrcap_2 = 0, nrcap_3 = 0,
                                  nAdult_1 = 0, nAdult_2 = 0, nAdult_3 = 0, nSubadult_1 = 0, nSubadult_2 = 0, nSubadult_3 = 0, nJuvenile_1 = 0, nJuvenile_2 = 0,
                                  nJuvenile_3 = 0, nUU_1 = 0, nUU_2 = 0, nUU_3 = 0, tooluse_1 = NA, tooluse_2 = NA, tooluse_3 = NA)

for (i in 1:nrow(gridseq_ocTU)) {
  ## at beginning have some kind of check if the sequenceID is already in the co-occurence dataframe, if so can skip everything
  if(sum(str_detect(cooccurrences_TU$cooc_ID, paste(gridseq_ocTU$sequenceID[i]))) == 0) {
    dist <- as.data.frame(subset(TUdistmat_all, rownames(TUdistmat_all) %in% gridseq_ocTU$locationfactor[i])) 
    cand_locs <- colnames(dist[,dist > 150]) # make list of candidate locations for co-occurrence (>150 m away)
    # filter to sequence that are at candidate location and on same day as sequence we're looking at 
    cand_seq <- gridseq_ocTU[gridseq_ocTU$locationfactor %in% cand_locs & gridseq_ocTU$seqday == gridseq_ocTU$seqday[i], c("sequenceID", "locationfactor", "seqday", "seq_start", "seq_end", "n", "nAdult", "nJuvenile","nSubadult", "nUU", "tooluse")]
    dist_m <- melt(dist)
    cand_seq$locationfactor <- as.character(cand_seq$locationfactor)
    dist_m$variable <- as.character(dist_m$variable)
    cand_seq <- left_join(cand_seq, dist_m, by = c("locationfactor" = "variable"))
    # see if there are any co-occurrences
    # if there is anything, then extract information from those sequences, both add to gridseq_ocTU dataframe, and to co-occurrence dataframe?
    if(nrow(cand_seq) > 0) {
      # identify the difftime cutoff for every row
      cand_seq$cutoff <- disttime2$time[findInterval(cand_seq$value, disttime2$distance)]
      cand_seq$flag <- ifelse(abs(difftime(gridseq_ocTU$seq_start[i], cand_seq$seq_start, unit = "s")) < cand_seq$cutoff, 1, 0)
      if(sum(cand_seq$flag) > 0)  {
        cand_seq$dtime <- difftime(gridseq_ocTU$seq_start[i], cand_seq$seq_start, unit = "s")
        cand_seq_t <- cand_seq[cand_seq$flag == 1,]
        cand_seq_t <- cand_seq_t[!duplicated(cand_seq_t$locationfactor),]
        if(nrow(cand_seq_t) > 0) {
          gridseq_ocTU$cooccurrence <- 1
          gridseq_ocTU$cooc_ID[i] <- ifelse(nrow(cand_seq_t) == 1, paste(gridseq_ocTU$sequenceID[i], cand_seq_t$sequenceID[1], sep = ","),
                                            paste(gridseq_ocTU$sequenceID[i], cand_seq_t$sequenceID[1], cand_seq_t$sequenceID[2], sep = ","))
          cooccurrences_TU[nrow(cooccurrences_TU) +1,] <- c(gridseq_ocTU$cooc_ID[i], paste(gridseq_ocTU$seq_start[i]), paste(gridseq_ocTU$seqday[i]), paste(gridseq_ocTU$locationfactor[i]), 
                                                            paste(cand_seq_t$locationfactor[1]), paste(cand_seq_t$locationfactor[2]), cand_seq_t$value[1], cand_seq_t$value[2], nrow(cand_seq_t), gridseq_ocTU$n[i], 
                                                            cand_seq_t$n[1], cand_seq_t$n[2], gridseq_ocTU$nAdult[i], cand_seq_t$nAdult[1], cand_seq_t$nAdult[2],gridseq_ocTU$nSubadult[i], cand_seq_t$nSubadult[1], 
                                                            cand_seq_t$nSubadult[2], gridseq_ocTU$nJuvenile[i], cand_seq_t$nJuvenile[1], cand_seq_t$nJuvenile[2], gridseq_ocTU$nUU[i], cand_seq_t$nUU[1], 
                                                            cand_seq_t$nUU[2], paste(gridseq_ocTU$tooluse[i]), paste(cand_seq_t$tooluse[1]), paste(cand_seq_t$tooluse[2]))
        }
      }
    }
  }
  print(i)
}

cooccurrences_TU <- cooccurrences_TU[-1,]
cooccurrences_TU[,7:22] <- as.numeric(unlist(cooccurrences_TU[,7:22]))

hist(cooccurrences_TU$distcam12)

## NTU 
gridseq_ocNTU <- gridseq_ocf[gridseq_ocf$gridtype == "NTU",]

# distance matrix is already made
NTUdistmat

# make blank co-occurrence info data frame? 
gridseq_ocNTU$cooccurrence <- 0
gridseq_ocNTU$cooc_ID <- NA

cooccurrences_NTU <- data.frame(cooc_ID = "seqid", seqstart = NA, seqday = NA, cam1 = NA, cam2 = NA, cam3 =NA, distcam12 = 0, distcam13 = 0, nrseq = 0, nrcap_1 = 0, nrcap_2 = 0, nrcap_3 = 0,
                                nAdult_1 = 0, nAdult_2 = 0, nAdult_3 = 0, nSubadult_1 = 0, nSubadult_2 = 0, nSubadult_3 = 0, nJuvenile_1 = 0, nJuvenile_2 = 0,
                                nJuvenile_3 = 0, nUU_1 = 0, nUU_2 = 0, nUU_3 = 0, tooluse_1 = NA, tooluse_2 = NA, tooluse_3 = NA)

for (i in 1:nrow(gridseq_ocNTU)) {
  ## at beginning have some kind of check if the sequenceID is already in the co-occurence dataframe, if so can skip everything
  if(sum(str_detect(cooccurrences_NTU$cooc_ID, paste(gridseq_ocNTU$sequenceID[i]))) == 0) {
    dist <- as.data.frame(subset(TUdistmat_all, rownames(TUdistmat_all) %in% gridseq_ocNTU$locationfactor[i])) 
    cand_locs <- colnames(dist[,dist > 150]) # make list of candidate locations for co-occurrence (>150 m away)
    # filter to sequence that are at candidate location and on same day as sequence we're looking at 
    cand_seq <- gridseq_ocNTU[gridseq_ocNTU$locationfactor %in% cand_locs & gridseq_ocNTU$seqday == gridseq_ocNTU$seqday[i], c("sequenceID", "locationfactor", "seqday", "seq_start", "seq_end", "n", "nAdult", "nJuvenile","nSubadult", "nUU", "tooluse")]
    dist_m <- melt(dist)
    cand_seq$locationfactor <- as.character(cand_seq$locationfactor)
    dist_m$variable <- as.character(dist_m$variable)
    cand_seq <- left_join(cand_seq, dist_m, by = c("locationfactor" = "variable"))
    # see if there are any co-occurrences
    # if there is anything, then extract information from those sequences, both add to gridseq_ocNTU dataframe, and to co-occurrence dataframe?
    if(nrow(cand_seq) > 0) {
      # identify the difftime cutoff for every row
      cand_seq$cutoff <- disttime2$time[findInterval(cand_seq$value, disttime2$distance)]
      cand_seq$flag <- ifelse(abs(difftime(gridseq_ocNTU$seq_start[i], cand_seq$seq_start, unit = "s")) < cand_seq$cutoff, 1, 0)
      if(sum(cand_seq$flag) > 0)  {
        cand_seq$dtime <- difftime(gridseq_ocNTU$seq_start[i], cand_seq$seq_start, unit = "s")
        cand_seq_t <- cand_seq[cand_seq$flag == 1,]
        cand_seq_t <- cand_seq_t[!duplicated(cand_seq_t$locationfactor),]
        if(nrow(cand_seq_t) > 0) {
          gridseq_ocNTU$cooccurrence <- 1
          gridseq_ocNTU$cooc_ID[i] <- ifelse(nrow(cand_seq_t) == 1, paste(gridseq_ocNTU$sequenceID[i], cand_seq_t$sequenceID[1], sep = ","),
                                             paste(gridseq_ocNTU$sequenceID[i], cand_seq_t$sequenceID[1], cand_seq_t$sequenceID[2], sep = ","))
          cooccurrences_NTU[nrow(cooccurrences_NTU) +1,] <- c(gridseq_ocNTU$cooc_ID[i], paste(gridseq_ocNTU$seq_start[i]), paste(gridseq_ocNTU$seqday[i]), paste(gridseq_ocNTU$locationfactor[i]), 
                                                              paste(cand_seq_t$locationfactor[1]), paste(cand_seq_t$locationfactor[2]), cand_seq_t$value[1], cand_seq_t$value[2], nrow(cand_seq_t), gridseq_ocNTU$n[i], 
                                                              cand_seq_t$n[1], cand_seq_t$n[2], gridseq_ocNTU$nAdult[i], cand_seq_t$nAdult[1], cand_seq_t$nAdult[2],gridseq_ocNTU$nSubadult[i], cand_seq_t$nSubadult[1], 
                                                              cand_seq_t$nSubadult[2], gridseq_ocNTU$nJuvenile[i], cand_seq_t$nJuvenile[1], cand_seq_t$nJuvenile[2], gridseq_ocNTU$nUU[i], cand_seq_t$nUU[1], 
                                                              cand_seq_t$nUU[2], paste(gridseq_ocNTU$tooluse[i]), paste(cand_seq_t$tooluse[1]), paste(cand_seq_t$tooluse[2]))
        }
      }
    }
  }
  print(i)
}

cooccurrences_NTU <- cooccurrences_NTU[-1,]
cooccurrences_NTU[,7:22] <- as.numeric(unlist(cooccurrences_NTU[,7:22]))

hist(cooccurrences_NTU$distcam12)



######## GAUSSIAN PROCESSES ############

## For Brendan: can start running from here
# load data in that I send you
#TUdistmat <- readRDS("gridanalyses/RDS/TUdistmat.rds")
#NTUdistmat <- readRDS("gridanalyses/RDS/NTUdistmat.rds")
#gridseq_ocf <- readRDS("gridanalyses/RDS/gridseq_ocf.rds")

dTU <- gridseq_ocf[gridseq_ocf$gridtype == "TU",]
dTU <- droplevels.data.frame(dTU)
# index cameras
dTU$camera <- as.numeric(dTU$locationfactor)

dNTU <- gridseq_ocf[gridseq_ocf$gridtype == "NTU",]
dNTU <- droplevels.data.frame(dNTU)
# index cameras
dNTU$camera <- as.numeric(dNTU$locationfactor)

## visualize the data we have

# TU and NTU separate
# per month (?)
# information captured:
# one row per camera, number of sequences in a day (or a hash per sequence). Color can reflect number of capuchins per sequence. 
# Could do it on time and show all the times of the sequences within a day (limit to daytime hours?)

# subset to things we need for graph
dTU_graph <- dTU[,c("seq_startday", "locationfactor", "n", "seq_length", "seq_start", "camera")]
dNTU_graph <- dNTU[,c("seq_startday", "locationfactor", "n", "seq_length", "seq_start", "camera")]

rbPal <- colorRampPalette(c("blue", "red"))
colors <- data.frame(n = 1:16, color = rbPal(16))

cameraTU_index <- data.frame(index = sort(unique(dTU_graph$camera)), cameralocation = sort(unique(dTU_graph$locationfactor)))
cameraNTU_index <- data.frame(index = sort(unique(dNTU_graph$camera)), cameralocation = sort(unique(dNTU_graph$locationfactor)))

#TU group - first half of deployment (until 18th of September 2022)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(camera ~ seq_start, data = dTU_graph[dTU_graph$seq_start < "2022-09-19",], pch = 21, yaxt = "n",  ylab = NA, xaxs ="i", col = colors$color[n], xlim = c(as.POSIXct("2022-05-12"), as.POSIXct("2022-09-20")))
axis(2, at=cameraTU_index$index , labels=cameraTU_index$cameralocation, las=2 , cex.axis=1 , tick=TRUE , lwd.ticks=0 , hadj=0.8)
legend("right",  inset=c(-0.15,0), pch=21, col=colors$color, as.character(colors$n), horiz=FALSE , cex=1 ,  xjust = 0.5, title = "# capuchins" )
title("Tool-users")

# TU group - second half of deployment (18th of september - end)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(camera ~ seq_start, data = dTU_graph[dTU_graph$seq_start > "2022-09-18",], pch = 21, yaxt = "n",  ylab = NA, xaxs ="i", col = colors$color[n], xlim = c(as.POSIXct("2022-09-17"), as.POSIXct("2023-01-26")))
axis(2, at=cameraTU_index$index , labels=cameraTU_index$cameralocation, las=2 , cex.axis=1 , tick=TRUE , lwd.ticks=0 , hadj=0.8)
legend("right",  inset=c(-0.15,0), pch=21, col=colors$color, as.character(colors$n), horiz=FALSE , cex=1 ,  xjust = 0.5, title = "# capuchins" )
title("Tool-users")

# NTU group - first half of deployment (until 18th of september)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(camera ~ seq_start, data = dNTU_graph[dNTU_graph$seq_start < "2022-09-19",], pch = 21, yaxt = "n",  ylab = NA, xaxs ="i", col = colors$color[n], xlim = c(as.POSIXct("2022-05-09"), as.POSIXct("2022-09-20")))
axis(2, at=cameraNTU_index$index , labels=cameraNTU_index$cameralocation, las=2 , cex.axis=1 , tick=TRUE , lwd.ticks=0 , hadj=0.8)
legend("right",  inset=c(-0.15,0), pch=21, col=colors$color, as.character(colors$n), horiz=FALSE , cex=1 ,  xjust = 0.5, title = "# capuchins" )
title("Non-tool-users")

# NTU group - second half of deployment (18th of september - end)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(camera ~ seq_start, data = dNTU_graph[dNTU_graph$seq_start > "2022-09-18",], pch = 21, yaxt = "n",  ylab = NA, xaxs ="i", col = colors$color[n], xlim = c(as.POSIXct("2022-09-17"), as.POSIXct("2023-01-27")))
axis(2, at=cameraNTU_index$index , labels=cameraNTU_index$cameralocation, las=2 , cex.axis=1 , tick=TRUE , lwd.ticks=0 , hadj=0.8)
legend("right",  inset=c(-0.15,0), pch=21, col=colors$color, as.character(colors$n), horiz=FALSE , cex=1 ,  xjust = 0.5, title = "# capuchins" )
title("Non-tool-users")

# make sure statistical rethinking is installed
require(rethinking)

# put distance matrices on a scale that makes sense
# now doing hundreds of meters
TUdistmat2 <- round(TUdistmat/100,2)
NTUdistmat2 <- round(NTUdistmat/100,2)

## work with number of capuchins per sequence (poisson). 

# try removing the zeros
# try keeping the zeros in


## without 0s
### Tool-users


# island example is collapsed to have one row per island. So we would need one row per camera. 
# he uses poisson, if we want poisson too we could do total number of capuchins seen (very crude but ok)
dTU_total <- aggregate(list(n = dTU$n, seq_length = dTU$seq_length), by = list(camera = dTU$camera, long = dTU$longitude, lat = dTU$latitude,  location = dTU$locationName), FUN = "sum")
## add distance to coast per camera
dist2coast_all <- read.csv("tide_analysis/allcams_gps.csv", header = TRUE)
dist2coast_all <- dist2coast_all[ , c(2,5)]

dTU_total <- left_join(dTU_total, dist2coast_all, by = c("location" = "camera_id"))
#standardize
dTU_total$distcoast_z <- as.numeric(scale(dTU_total$distcoast, center = TRUE, scale = TRUE))

# more accurate would be average number of capuchins/rate. So this would be a gamma. Try that too
# could do average number of capuchins seen? over all sequences?
## IMPORTANT, HAVE ONE SEQUENCE 01b1e2b3-3a65-4165-a1ed-088f903de735 with a seq_length of 0, because it is only one picture. Not sure how this happened
# for now just make it last 1 second
dTU$seq_length[dTU$seq_length == 0] <- 1
dTU_meannumber <- aggregate(dTU$n, by = list(camera = dTU$camera, long = dTU$longitude, lat = dTU$latitude, location = dTU$locationName), FUN = "mean")
# later would probably need to do something like "capuchins per second"
dTU$rate <- dTU$n/dTU$seq_length
dTU_meanrate <- aggregate(dTU$rate, by = list(camera = dTU$camera,  long = dTU$longitude, lat = dTU$latitude,  location = dTU$locationName), FUN = "mean")

## first just model the spatial covariance, null-model without any predictors
# simulate priors
n <- 30
etasq <- rexp(n,2)
rhosq <- rexp(n,0.5)

plot(NULL, xlim = c(0,10), ylim = c(0,2), 
     xlab = "distance (hundred m)",
     ylab = "covariance")

for(i in 1:n){
  curve(etasq[i]*exp(-rhosq[i]*x^2),
        add = TRUE, lwd = 4,
        col = col.alpha(2,0.5))
}


dat_list <- list(
  N = dTU$n,
  C = dTU$camera,
  D = TUdistmat2 )

mTdist_TU <- ulam(
  alist(
    N ~ dpois(lambda),
    log(lambda) <- abar + a[C],
    vector[24]:a ~ multi_normal(0, K),
    matrix[24,24]:K <- cov_GPL2(D, etasq, rhosq, 0.01),
    abar ~ normal(3, 0.5),
    etasq ~ dexp(2),
    rhosq ~ dexp(0.5)
  ), data = dat_list, chains = 4, cores = 4, iter = 4000)

# save and load model
#saveRDS(mTdist_TU, "gridanalyses/RDS/mTdist_TU.rds")
#mTdist_TU <- readRDS("gridanalyses/RDS/mTdist_TU.rds")

precis(mTdist_TU, 2)

# visualize posterior
post <- extract.samples(mTdist_TU)

# plot posterior median covariance function
plot(NULL, xlab = "distance (hundred m)", ylab = "covariance",
     xlim = c(0,10), ylim = c(0,3))

# compute posterior mean covariance
x_seq <- seq(from=0, to = 10, length.out = 100)
pmcov <- sapply(x_seq, function(x) post$etasq*exp(-post$rhosq*x^2))
pmcov_mu <- apply(pmcov, 2, mean)
lines(x_seq, pmcov_mu, lwd = 3)

# prior in black
for(i in 1:n){
  curve(etasq[i]*exp(-rhosq[i]*x^2),
        add = TRUE, lwd = 2,
        col = col.alpha("black",0.25))
}

# plot 60 functions sampled from the posterior
for (i in 1:50) {
  curve( post$etasq[i]*exp(-post$rhosq[i]*x^2) , add=TRUE , lwd = 2,
         col=col.alpha("red",0.5) )
}

# compute posterior median covariance among societies
K <- matrix(0, nrow = 24, ncol = 24)
for (i in 1:24)
  for (j in 1:24)
    K[i,j] <- median(post$etasq) *
  exp( - median(post$rhosq) * TUdistmat2[i,j]^2)
diag(K) <- median(post$etasq) + 0.01

# convert to correlation matrix
Rho <- round(cov2cor(K), 2)
# add row/col names for convenience
colnames(Rho) <- colnames(TUdistmat2)
rownames(Rho) <- colnames(Rho)

# plot raw data and labels
plot(dTU_total$long , dTU_total$lat , xlab="longitude" , ylab="latitude" ,
     col="red"  , pch=16, xlim = c(-81.824, -81.815))
labels <- as.character(dTU_total$location)
text( dTU_total$long , dTU_total$lat , labels=labels , cex=0.7 , pos=c(2,4,3,3,4,1,3,2,4,2) )
# overlay lines shaded by Rho
for( i in 1:24 )
  for ( j in 1:24 )
    if ( i < j )
      lines( c( dTU_total$long[i],dTU_total$long[j] ) , c( dTU_total$lat[i],dTU_total$lat[j] ) ,
             lwd=2 , col=col.alpha("black",Rho[i,j]^2) )

## including "population", in our case, total sequence length
dat_list <- list(
  N = dTU$n,
  C = dTU$camera,
  P = dTU$seq_length,
  D = TUdistmat2 )

mTdist_TU2 <- ulam(
  alist(
    N ~ dpois(lambda),
    lambda <- (abar*P^b/g)*exp(a[C]),
    vector[24]:a ~ multi_normal(0, K),
    transpars>matrix[24,24]:K <- 
      cov_GPL2(D, etasq, rhosq, 0.01),
    c(abar, b,g) ~ dexp(1),
    etasq ~ dexp(2),
    rhosq ~ dexp(0.5)
  ), data = dat_list, chains = 4, cores = 4, iter = 4000)

# save and load model
#saveRDS(mTdist_TU2, "gridanalyses/RDS/mTdist_TU2.rds")
#mTdist_TU2 <- readRDS("gridanalyses/RDS/mTdist_TU2.rds")

precis(mTdist_TU2, 3)

# visualize posterior
post2 <- extract.samples(mTdist_TU2)

# plot posterior median covariance function
plot(NULL, xlab = "distance (hundred m)", ylab = "covariance",
     xlim = c(0,10), ylim = c(0,3))

# prior in black
for(i in 1:n){
  curve(etasq[i]*exp(-rhosq[i]*x^2),
        add = TRUE, lwd = 2,
        col = col.alpha("black",0.25))
}

## null model in blue
for (i in 1:50) {
  curve( post$etasq[i]*exp(-post$rhosq[i]*x^2) , add=TRUE , lwd = 2,
         col=col.alpha("blue",0.5) )
}

# plot 60 functions sampled from the posterior
for (i in 1:50) {
  curve( post2$etasq[i]*exp(-post2$rhosq[i]*x^2) , add=TRUE , lwd = 2,
         col=col.alpha("red",0.5) )
}

# compute posterior median covariance among societies
K <- matrix(0, nrow = 24, ncol = 24)
for (i in 1:24)
  for (j in 1:24)
    K[i,j] <- median(post2$etasq) *
  exp( - median(post2$rhosq) * TUdistmat2[i,j]^2)
diag(K) <- median(post2$etasq) + 0.01

# convert to correlation matrix
Rho <- round(cov2cor(K), 2)
# add row/col names for convenience
colnames(Rho) <- colnames(TUdistmat2)
rownames(Rho) <- colnames(Rho)

# plot raw data and labels
plot(dTU_total$long , dTU_total$lat , xlab="longitude" , ylab="latitude" ,
     col="red"  , pch=16, xlim = c(-81.824, -81.815))
labels <- as.character(dTU_total$location)
text( dTU_total$long , dTU_total$lat , labels=labels , cex=0.7 , pos=c(2,4,3,3,4,1,3,2,4,2) )
# overlay lines shaded by Rho
for( i in 1:24 )
  for ( j in 1:24 )
    if ( i < j )
      lines( c( dTU_total$long[i],dTU_total$long[j] ) , c( dTU_total$lat[i],dTU_total$lat[j] ) ,
             lwd=2 , col=col.alpha("black",Rho[i,j]^2) )

# now on log seqlength scale , still want to make this plot
# youtube video https://www.youtube.com/watch?v=PIuqxOBJqLU around 35:00
# statistical rethinking around page 475 (figure 14.12)
# compute posterior median relationship, ignoring distance 
max(log(dTU_total$seq_length))

# something goes wrong with code below
logpop.seq <- seq( from=4 , to=9 , length.out=30 )
lambda <- sapply( logpop.seq , function(lp) exp( post2$a + post2$bp*lp ) )
lambda.median <- apply( lambda , 2 , median )
lambda.PI80 <- apply( lambda , 2 , PI , prob=0.8 )
# plot raw data and labels
plot( d$logpop , d$total_tools , col=rangi2 , cex=psize , pch=16 ,
      xlab="log population" , ylab="total tools" )
text( d$logpop , d$total_tools , labels=labels , cex=0.7 ,
      pos=c(4,3,4,2,2,1,4,4,4,2) )
# display posterior predictions
lines( logpop.seq , lambda.median , lty=2 )
lines( logpop.seq , lambda.PI80[1,] , lty=2 )
lines( logpop.seq , lambda.PI80[2,] , lty=2 )
# overlay correlations
for( i in 1:10 )
  for ( j in 1:10 )
    if ( i < j )
      lines( c( d$logpop[i],d$logpop[j] ) ,
             c( d$total_tools[i],d$total_tools[j] ) ,
             lwd=2 , col=col.alpha("black",Rho[i,j]^2) )

## including distance to coast as predictor
dat_list <- list(
  N = dTU_total$n,
  C = 1:24,
  P = dTU_total$seq_length,
  D = TUdistmat2,
  X = dTU_total$distcoast_z)


## discuss with Brendan how to incorporate distcoast in the ulam
# now just added bX and X to lambda but that doesn't seem to be right?
mTdist_TU3 <- ulam(
  alist(
    N ~ dpois(lambda),
    lambda <- (abar*P^b/g)*exp(a[C]) + bX*X,
    vector[24]:a ~ multi_normal(0, K),
    transpars>matrix[24,24]:K <- 
      cov_GPL2(D, etasq, rhosq, 0.01),
    c(abar, b,g) ~ dexp(1),
    bX ~ normal(0, 0.5),
    etasq ~ dexp(2),
    rhosq ~ dexp(0.5)
  ), data = dat_list, chains = 4, cores = 4, iter = 4000)

precis(mTdist_TU3, 3)

# visualize posterior
post3 <- extract.samples(mTdist_TU3)

# plot posterior median covariance function
plot(NULL, xlab = "distance (hundred m)", ylab = "covariance",
     xlim = c(0,10), ylim = c(0,2))

# prior in black
for(i in 1:n){
  curve(etasq[i]*exp(-rhosq[i]*x^2),
        add = TRUE, lwd = 2,
        col = col.alpha("black",0.5))
}

## null model in blue
for (i in 1:50) {
  curve( post$etasq[i]*exp(-post$rhosq[i]*x^2) , add=TRUE , lwd = 2,
         col=col.alpha("blue",0.5) )
}

# population added
for (i in 1:50) {
  curve( post2$etasq[i]*exp(-post2$rhosq[i]*x^2) , add=TRUE , lwd = 2,
         col=col.alpha("red",0.5) )
}

# population and distcoast
for (i in 1:50) {
  curve( post3$etasq[i]*exp(-post3$rhosq[i]*x^2) , add=TRUE , lwd = 2,
         col=col.alpha("green",0.5) )
}

# compute posterior median covariance among societies
K <- matrix(0, nrow = 24, ncol = 24)
for (i in 1:24)
  for (j in 1:24)
    K[i,j] <- median(post3$etasq) *
  exp( - median(post3$rhosq) * TUdistmat2[i,j]^2)
diag(K) <- median(post3$etasq) + 0.01

# convert to correlation matrix
Rho <- round(cov2cor(K), 2)
# add row/col names for convenience
colnames(Rho) <- colnames(TUdistmat2)
rownames(Rho) <- colnames(Rho)

# plot raw data and labels
plot(dTU_total$long , dTU_total$lat , xlab="longitude" , ylab="latitude" ,
     col="red"  , pch=16, xlim = c(-81.824, -81.815))
labels <- as.character(dTU_total$location)
text( dTU_total$long , dTU_total$lat , labels=labels , cex=0.7 , pos=c(2,4,3,3,4,1,3,2,4,2) )
# overlay lines shaded by Rho
for( i in 1:24 )
  for ( j in 1:24 )
    if ( i < j )
      lines( c( dTU_total$long[i],dTU_total$long[j] ) , c( dTU_total$lat[i],dTU_total$lat[j] ) ,
             lwd=2 , col=col.alpha("black",Rho[i,j]^2) )

# not sure how to measure if distcoast affects covariance?


### Non-tool-users

# island example is collapsed to have one row per island. So we would need one row per camera. 
# he uses poisson, if we want poisson too we could do total number of capuchins seen (very crude but ok)
dNTU_total <- aggregate(list(n = dNTU$n, seq_length = dNTU$seq_length), by = list(camera = dNTU$camera, long = dNTU$longitude, lat = dNTU$latitude,  location = dNTU$locationName), FUN = "sum")

dNTU_total <- left_join(dNTU_total, dist2coast_all, by = c("location" = "camera_id"))
#standardize
dNTU_total$distcoast_z <- as.numeric(scale(dNTU_total$distcoast, center = TRUE, scale = TRUE))


# more accurate would be average number of capuchins/rate. So this would be a gamma. Try that too
# could do average number of capuchins seen? over all sequences?
dNTU_meannumber <- aggregate(dNTU$n, by = list(camera = dNTU$camera, long = dNTU$longitude, lat = dNTU$latitude, location = dNTU$locationName), FUN = "mean")
# later would probably need to do something like "capuchins per second"
dNTU$rate <- dNTU$n/dNTU$seq_length
dNTU_meanrate <- aggregate(dNTU$rate, by = list(camera = dNTU$camera,  long = dNTU$longitude, lat = dNTU$latitude,  location = dNTU$locationName), FUN = "mean")

## first just model the spatial covariance, null-model without any predictors
# simulate priors
dat_list <- list(
  N = dNTU$n,
  C = dNTU$camera,
  D = NTUdistmat2 )

mTdist_NTU <- ulam(
  alist(
    N ~ dpois(lambda),
    log(lambda) <- abar + a[C],
    vector[25]:a ~ multi_normal(0, K),
    matrix[25,25]:K <- cov_GPL2(D, etasq, rhosq, 0.01),
    abar ~ normal(3, 0.5),
    etasq ~ dexp(2),
    rhosq ~ dexp(0.5)
  ), data = dat_list, chains = 4, cores = 4, iter = 4000)

# save and load model
#saveRDS(mTdist_NTU, "gridanalyses/RDS/mTdist_NTU.rds")
#mTdist_NTU <- readRDS("gridanalyses/RDS/mTdist_NTU.rds")

precis(mTdist_NTU, 2)

# visualize posterior
post <- extract.samples(mTdist_NTU)

# plot posterior median covariance function
plot(NULL, xlab = "distance (hundred m)", ylab = "covariance",
     xlim = c(0,10), ylim = c(0,2))

# compute posterior mean covariance
x_seq <- seq(from=0, to = 10, length.out = 100)
pmcov <- sapply(x_seq, function(x) post$etasq*exp(-post$rhosq*x^2))
pmcov_mu <- apply(pmcov, 2, mean)
lines(x_seq, pmcov_mu, lwd = 3)

# prior in black
for(i in 1:n){
  curve(etasq[i]*exp(-rhosq[i]*x^2),
        add = TRUE, lwd = 2,
        col = col.alpha("black",0.25))
}

# plot 60 functions sampled from the posterior
for (i in 1:50) {
  curve( post$etasq[i]*exp(-post$rhosq[i]*x^2) , add=TRUE , lwd = 2,
         col=col.alpha("red",0.5) )
}

# compute posterior median covariance among societies
K <- matrix(0, nrow = 25, ncol = 25)
for (i in 1:25)
  for (j in 1:25)
    K[i,j] <- median(post$etasq) *
  exp( - median(post$rhosq) * NTUdistmat2[i,j]^2)
diag(K) <- median(post$etasq) + 0.01

# convert to correlation matrix
Rho <- round(cov2cor(K), 2)
# add row/col names for convenience
colnames(Rho) <- colnames(NTUdistmat2)
rownames(Rho) <- colnames(Rho)

# plot raw data and labels
plot(dNTU_total$long , dNTU_total$lat , xlab="longitude" , ylab="latitude" ,
     col="red"  , pch=16, xlim = c(-81.798, -81.790))
labels <- as.character(dNTU_total$location)
text( dNTU_total$long , dNTU_total$lat , labels=labels , cex=0.7 , pos=c(2,4,3,3,4,1,3,2,4,2) )
# overlay lines shaded by Rho
for( i in 1:25 )
  for ( j in 1:25 )
    if ( i < j )
      lines( c( dNTU_total$long[i],dNTU_total$long[j] ) , c( dNTU_total$lat[i],dNTU_total$lat[j] ) ,
             lwd=2 , col=col.alpha("black",Rho[i,j]^2) )


## including "population", in our case, total sequence length
dat_list <- list(
  N = dNTU$n,
  C = dNTU$camera,
  P = dNTU$seq_length,
  D = NTUdistmat2 )

mTdist_NTU2 <- ulam(
  alist(
    N ~ dpois(lambda),
    lambda <- (abar*P^b/g)*exp(a[C]),
    vector[25]:a ~ multi_normal(0, K),
    transpars>matrix[25,25]:K <- 
      cov_GPL2(D, etasq, rhosq, 0.01),
    c(abar, b,g) ~ dexp(1),
    etasq ~ dexp(2),
    rhosq ~ dexp(0.5)
  ), data = dat_list, chains = 4, cores = 4, iter = 4000)

#saveRDS(mTdist_NTU2, "gridanalyses/RDS/mTdist_NTU2.rds")
#mTdist_NTU2 <- readRDS("gridanalyses/RDS/mTdist_NTU2.rds")

precis(mTdist_NTU2, 2)

# visualize posterior
post2 <- extract.samples(mTdist_NTU2)

# plot posterior median covariance function
plot(NULL, xlab = "distance (hundred m)", ylab = "covariance",
     xlim = c(0,10), ylim = c(0,2))

# compute posterior mean covariance
x_seq <- seq(from=0, to = 10, length.out = 100)
pmcov <- sapply(x_seq, function(x) post$etasq*exp(-post$rhosq*x^2))
pmcov_mu <- apply(pmcov, 2, mean)
lines(x_seq, pmcov_mu, lwd = 3)

# prior in black
for(i in 1:n){
  curve(etasq[i]*exp(-rhosq[i]*x^2),
        add = TRUE, lwd = 2,
        col = col.alpha("black",0.25))
}

## null model in blue
for (i in 1:50) {
  curve( post$etasq[i]*exp(-post$rhosq[i]*x^2) , add=TRUE , lwd = 2,
         col=col.alpha("blue",0.5) )
}

# plot 60 functions sampled from the posterior
for (i in 1:50) {
  curve( post2$etasq[i]*exp(-post2$rhosq[i]*x^2) , add=TRUE , lwd = 2,
         col=col.alpha("red",0.5) )
}

# compute posterior median covariance among societies
K <- matrix(0, nrow = 25, ncol = 25)
for (i in 1:25)
  for (j in 1:25)
    K[i,j] <- median(post2$etasq) *
  exp( - median(post2$rhosq) * NTUdistmat2[i,j]^2)
diag(K) <- median(post2$etasq) + 0.01

# convert to correlation matrix
Rho <- round(cov2cor(K), 2)
# add row/col names for convenience
colnames(Rho) <- colnames(NTUdistmat2)
rownames(Rho) <- colnames(Rho)

# plot raw data and labels
plot(dNTU_total$long , dNTU_total$lat , xlab="longitude" , ylab="latitude" ,
     col="red"  , pch=16, xlim = c(-81.798, -81.790))
labels <- as.character(dNTU_total$location)
text( dNTU_total$long , dNTU_total$lat , labels=labels , cex=0.7 , pos=c(2,4,3,3,4,1,3,2,4,2) )
# overlay lines shaded by Rho
for( i in 1:24 )
  for ( j in 1:24 )
    if ( i < j )
      lines( c( dNTU_total$long[i],dNTU_total$long[j] ) , c( dNTU_total$lat[i],dNTU_total$lat[j] ) ,
             lwd=2 , col=col.alpha("black",Rho[i,j]^2) )

# so my first interpretation of this is that including distance only already doesn't show much covariance between cameras for the NTU but does for TU
# hmmmm
# is then interpretation correct that for NTU close cameras are similar and far apart are not similar and for TU there is more covariance between far cameras too?
# look at this in more detail later


### TRIALS #####

### looking at how many 1 individual sequences we get 
hist(gridsequence_c$n[gridsequence_c$gridtype == "TU"])
hist(gridsequence_c$n[gridsequence_c$gridtype == "NTU"])

## this is on all data tool site
# level of observation
str(agouticlean)
agouticlean_jt <- agouticlean[agouticlean$island == "Jicaron" & agouticlean$tool_site == 1,]
v <- as.data.frame(ftable(agouticlean_jt$individualID))
# IDs where we have many observations
ids <- as.character(v$Var1[which(v$Freq > 20)])

agouticlean_jti <- agouticlean_jt[agouticlean_jt$individualID %in% ids,]
# number of 1 sequences
ggplot(data = agouticlean_jti, aes(x = n)) + geom_histogram(stat = "count") + facet_wrap(~name, scales = "free")

# hist of sequences depending on age-sex
# filter to age-sex that we have enough data for
agesex <- as.data.frame(ftable(agouticlean_jt$agesexF)) 
agouticlean_jtas <- agouticlean_jt[agouticlean_jt$agesexF %in% agesex$Var1[which(agesex$Freq > 50)],]
ggplot(data = agouticlean_jtas, aes(x = n)) + geom_histogram(stat = "count") + facet_wrap(~agesexF, scales = "free")

# only grid data
ggplot(data = gridclean_c[gridclean_c$gridtype == "TU",], aes(x = n)) + geom_histogram(stat = "count") + facet_wrap(~agesexF, scales = "free")
ggplot(data = gridclean_c[gridclean_c$gridtype == "NTU",], aes(x = n)) + geom_histogram(stat = "count") + facet_wrap(~agesexF, scales = "free")

# also just out of interest, when we have the ID data, do we see individuals more at specific locations?
ggplot(data = agouticlean_jti, aes(x = locationName)) + geom_histogram(stat = "count") + facet_wrap(~name)

# trying out looking at explicitly spatial model

agoutiseq_jt_order <- agoutiseq_jt[order(agoutiseq_jt$seq_start),]
test <- agoutiseq_jt_order[15000:16000, c("longitude", "latitude", "uniqueloctag", "seqday", "seq_start")]

ftable(test$seqday)
unique(test$uniqueloctag)
plot(test$longitude, test$latitude, type = "l")

#one day
test2 <- agoutiseq_jt_order[agoutiseq_jt_order$seqday == "2022-02-20", c("longitude", "latitude", "locationfactor")]

devtools::install_version("DiagTest3Grp",version="1.6")

P.move = df2move(test2, proj = "+init=EPSG:32616 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", x = "longitude", y = "latitude", time = "seq_start")
m <- align_move(P.move, res = 3, unit = "secs")
frames <- frames_spatial(m, 
                         map_service = "mapbox", map_type = "satellite", alpha = 0.5, path_size = 1.5, trace_show = TRUE, map_token = "pk.eyJ1IjoiZ29vbmVyNjE5NSIsImEiOiJja2cxMHB0ZG4wY3M1Mnhxd2QxMzFuZHdnIn0.k8zIUKu6KcV5ZXne3a04kg") %>%
  add_labels(title = "Herd directionality in Oryx exposed to simulated predators", x = "Longitude", y = "Latitude") %>% # add some customizations, such as axis labels
  add_northarrow(colour = "white", position = "bottomleft") %>% 
  add_scalebar(colour = "black", position = "upperright", distance = 0.15) %>% 
  add_timestamps(m, type = "label") %>%
  add_progress()
animate_frames(frames, out_file = "Oryx.movements.mp4")

####### CHRONOLOGICAL SPACE #####
# looking how far way each subsequent sequence detection is from previous one

# need
TUdistmat
NTUdistmat
head(gridseq_ocTU)
head(gridseq_ocNTU)

# this is all capuchin sightings of TU and NTU grid with unfamiliars removed (important!)

# step one: order chronologically
gridseq_ocTUc <- gridseq_ocTU[order(gridseq_ocTU$seq_start),]  
gridseq_ocNTUc <- gridseq_ocNTU[order(gridseq_ocNTU$seq_start),]

# make variable of distance to previous line
# make blank variable
gridseq_ocTUc$space <- NA

for(i in 2:nrow(gridseq_ocTUc)) {
  gridseq_ocTUc$space[i] <- TUdistmat[row.names(TUdistmat) == gridseq_ocTUc$locationName[i], colnames(TUdistmat) == gridseq_ocTUc$locationName[i-1]]
}

ftable(gridseq_ocTUc$space)

test <- gridseq_ocTUc[c("locationName", "seq_start", "space"),]
gridseq_ocTUc$hour
gridseq_ocTUc$hms <- as_hms(gridseq_ocTUc$seq_start)
head(gridseq_ocTUc$hms)


# make variable of distance to previous line
# make blank variable
gridseq_ocNTUc$space <- NA

for(i in 2:nrow(gridseq_ocNTUc)) {
  gridseq_ocNTUc$space[i] <- NTUdistmat[row.names(NTUdistmat) == gridseq_ocNTUc$locationName[i], colnames(NTUdistmat) == gridseq_ocNTUc$locationName[i-1]]
}

ftable(gridseq_ocNTUc$space)

test <- gridseq_ocNTUc[c("locationName", "seq_start", "space"),]
gridseq_ocNTUc$hms <- as_hms(gridseq_ocNTUc$seq_start)
head(gridseq_ocNTUc$hms)
ggplot(data = gridseq_ocNTUc, aes(x = hms, y = space, col = seqday)) + geom_point() 

ggplot(data= gridseq_ocTUc, aes(x = hms, y = space, col = seqday)) + geom_point() + geom_smooth() 
ggplot(data= gridseq_ocTUc, aes(x = hour, y = space, col = seqday)) +  stat_summary(data = gridseq_ocTUc, aes(x = hour, y = space), fun = mean, geom = "point", inherit.aes = FALSE) 
ggplot(data= gridseq_ocNTUc, aes(x = hour, y = space, col = seqday)) +  stat_summary(data = gridseq_ocNTUc, aes(x = hour, y = space), fun = mean, geom = "point", inherit.aes = FALSE) 

## THIS NEEDS SOME WORK. i THINK YOU ALWAYS WANT TO EXCLUDE THE FIRST ONE OF A DAY? 
# OR SET THE FIRST ONE OF A DAY TO NA AND ONLY LOOK AT DISTANCES WITHIN A DAY?
## need to somehow do this within a day and plot than the average of all days
# think about this. Already interesting how it looks different though, but have to wonder why that is

# maybe number observations per day and then have plotted average distance?
gridseq_ocNTUc$obsnumber <- as.numeric(ave(gridseq_ocNTUc$seqday, gridseq_ocNTUc$seqday, FUN = seq_along))
gridseq_ocTUc$obsnumber <- as.numeric(ave(gridseq_ocTUc$seqday, gridseq_ocTUc$seqday, FUN = seq_along))

ggplot(data= gridseq_ocTUc[!gridseq_ocTUc$obsnumber == 1,], aes(x = obsnumber, y = space)) + geom_point() + geom_smooth() 
ggplot(data= gridseq_ocNTUc[!gridseq_ocNTUc$obsnumber == 1,], aes(x = obsnumber, y = space)) + geom_point() + geom_smooth() 

ggplot(data= gridseq_ocTUc[!gridseq_ocTUc$obsnumber == 1,], aes(x = obsnumber, y = space)) + geom_point() + geom_point(aes(col = seqday), alpha = 0.5) + 
  stat_summary(data = gridseq_ocTUc[!gridseq_ocTUc$obsnumber == 1,], aes(x = obsnumber, y = space), fun = mean, geom = "point", inherit.aes = FALSE, size = 3, shape = 15) +
  ggtitle("Tool-Users")

ggplot(data= gridseq_ocNTUc[!gridseq_ocNTUc$obsnumber == 1,], aes(x = obsnumber, y = space)) + geom_point(aes(col = seqday), alpha = 0.5) +
  stat_summary(data = gridseq_ocNTUc[!gridseq_ocNTUc$obsnumber == 1,], aes(x = obsnumber, y = space), fun = mean, geom = "point", inherit.aes = FALSE, size = 3, shape = 15) +
  ggtitle("Non-tool-users")

## add in TIME to previous sighting
gridseq_ocTUc$time <- NA

for(i in 2:nrow(gridseq_ocTUc)) {
  gridseq_ocTUc$time[i] <- difftime(gridseq_ocTUc$seq_start[i],gridseq_ocTUc$seq_start[i-1], units = "m")
}

hist(gridseq_ocTUc$time)

gridseq_ocNTUc$time <- NA

for(i in 2:nrow(gridseq_ocNTUc)) {
  gridseq_ocNTUc$time[i] <- difftime(gridseq_ocNTUc$seq_start[i],gridseq_ocNTUc$seq_start[i-1], units = "m")
}

hist(gridseq_ocNTUc$time)


ggplot(data= gridseq_ocTUc[!gridseq_ocTUc$obsnumber == 1,], aes(x = obsnumber, y = space)) + geom_point() + geom_point(aes(col = time), alpha = 0.3, size = 3) + 
  stat_summary(data = gridseq_ocTUc[!gridseq_ocTUc$obsnumber == 1,], aes(x = obsnumber, y = space), fun = mean, geom = "point", inherit.aes = FALSE, size = 3, color = "red", shape = 15) +
  ggtitle("Tool-Users") + scale_colour_viridis_c() + 
  labs(x = "Number of consecutive observation per day", y = "Distance (meters) to previous sighting", color = "Time (seconds) to previous sighting") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) + ylim(0,850)

ggplot(data= gridseq_ocNTUc[!gridseq_ocNTUc$obsnumber == 1,], aes(x = obsnumber, y = space)) + geom_point(aes(col = time), alpha = 0.3, size = 3) +
  stat_summary(data = gridseq_ocNTUc[!gridseq_ocNTUc$obsnumber == 1,], aes(x = obsnumber, y = space), fun = mean, geom = "point", inherit.aes = FALSE, size = 3, color = "red", shape = 15) +
  ggtitle("Non-tool-users") +scale_colour_viridis_c() + 
  labs(x = "Number of consecutive observation per day", y = "Distance (meters) to previous sighting", color = "Time (seconds) to previous sighting") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) +ylim(0,850)

test <- gridseq_ocNTUc[c("seq_start","space","time")]

require(viridis)
inferncol <- viridis_pal(option = "B")(10)
mybreaks <- seq(0, 0.0001, length.out = 100)
mybreaks <- mybreaks[c(0:10,100)]
breaklabel <- function(x){
  labels<- paste0(mybreaks[1:10], "-", mybreaks[2:11])
  labels[1:x]
}



head(gridseq_ocNTUc)
ggplot(data = gridseq_ocNTUc[!gridseq_ocNTUc$obsnumber == 1 & gridseq_ocNTUc$space < 201 & gridseq_ocNTUc$time < 201,], aes(x = space, y = time)) +
  geom_density2d()
ggplot(data = gridseq_ocTUc[!gridseq_ocTUc$obsnumber == 1 & gridseq_ocTUc$space < 201 & gridseq_ocTUc$time < 201,], aes(x = space, y = time)) +
  geom_density2d()

ggplot(data = gridseq_ocTUc[!gridseq_ocTUc$obsnumber == 1,], aes(x = space, y = time)) + 
  geom_density2d_filled(breaks = mybreaks, show.legend = TRUE) + scale_fill_manual(values = inferncol, name = "Change nr of capuchins") 

ggplot(data = gridseq_ocNTUc[!gridseq_ocNTUc$obsnumber == 1,], aes(x = space, y = time)) + 
  geom_density2d_filled(breaks = mybreaks, show.legend = TRUE) + scale_fill_manual(values = inferncol, name = "Change nr of capuchins") 

# change binning so it has above certain number all yellow and other colors for the rest

ggplot(data = gridseq_ocTUc[!gridseq_ocTUc$obsnumber == 1,], aes(x = space, y = time)) + geom_point(alpha = 0.05)
ggplot(data = gridseq_ocNTUc[!gridseq_ocNTUc$obsnumber == 1,], aes(x = space, y = time)) + geom_point(alpha = 0.05)

# color reflect frequency/intensity? 

?stat_summary


# plot real data as points to get idea of where data is



######## MAP ############
library(mapview)
library(rgdal)
library(sf)

camcoords <- gridsequence[!duplicated(gridsequence$locationName), c("locationName", "longitude", "latitude")]
all_cams <- st_as_sf(camcoords , coords = c("longitude", "latitude"), crs = 4326) #do it again if rading csv
all_cams_map <- mapview(all_cams , col.regions="black"  )
all_cams_map
