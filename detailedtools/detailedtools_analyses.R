## Detailed tool use analyses -- Analyses
## MPI-AB; Z Goldsborough

## Analysing efficiency, technique and selectivity in tool use

## packages needed
library(stringr)
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

## NOTES ###
# If we have enough data, ideally exclude sequences that are split across multiple videos (split == TRUE)
# Similarly, if we have enough data with identified individuals, only work on that dataset so individual can be a random effect

### Load datasets, cleaned in detailedtools.R script ####
#detseq <- readRDS("detailedtools/RDS/detseq.rds")
#dettools_r2 <- readRDS("detailedtools/RDS/dettools_r2.rds")

# for most analyses, will only look at sequences that were successful (so outcome = opened)
# also filter the items to only the most abundant ones (for now)
ftable(detseq$item)
detseq_o <- detseq[detseq$outcome == "opened" & detseq$item %in% c("almendrabrown", "almendragreen", "almendraunknown", "almendrared"),]

### Diagnostics #### 

# first let's examine some general things, and see whether we have good inter-rater agreement and how much data we have
ggplot(detseq_o, aes(x=coder, y=n_pounds)) + 
  geom_violin()

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

### Efficiency ####
## Comparing efficiency between age classes on opened sequences, multiple measures
# 1. Sequence duration (seconds) 
# 2. Number of pounds to open item 
# 3. Number of misstrikes
# 4. Number of repositions

## Other things that likely affect these variables are:
# The specific item being processed
# The type of anvil (wood or stone)
# The hammerstone (when we have reliable IDs, include hammerstone ID as random effect)
# Identity of tool user (when we only use subset with identified individuals, random effect)
# Whether sequence is split across two videos or not (with pounds miss some, durations are longer)

### IF SEQUENCE IS SPLIT then we can "reliably" still look at the duration of the sequence in seconds (though likely overestimation?)
# but to be safe probably shouldnt include these when looking at number of pounds as we miss at least 1 (maybe more) pounds in between
# maybe exclude them in general if we have enough data, think about it
# definitely include whether it was split yes/no as a random effect/fixed effect 

nrow(detseq_o)
# so at the moment have 2184 sequences
# filter to the best of the best (where individuals are known), of course should still check as many unidentified ones as possible
# also exclude splits
detseq_oi <- detseq_o[!detseq_o$subjectID %in% c("adultmale", "subadultmale", "juvenileunknown") & detseq_o$split == FALSE,]
# exclude ones were we missed all pounds (so n_pounds = NA)
detseq_oi <- detseq_oi[!is.na(detseq_oi$n_pounds) == TRUE,]
# then end up with 2839 sequences

#### 1. Sequence duration #####

## determine distribution, what family is best?
descdist(detseq_oi$seqduration)
hist(detseq_oi$seqduration)

testdist1.1 <- fitdist(detseq_oi$seqduration, "norm")
plot(testdist1.1)

testdist1.2 <- fitdist(detseq_oi$seqduration, "gamma")
plot(testdist1.2)

# Gamma appears best (for now)

## Model 1: Duration depending on age, including item, anviltype, and individual ID as random effect
m_e1 <- brm(seqduration ~ Age + item*anviltype + (1|subjectID), data = detseq_oi, iter = 1000, chain = 2, core = 2, backend = "cmdstanr", family = "gamma")
# saving model if you dont want to have to run it again
#saveRDS(m_e1, "detailedtools/RDS/m_e1.rds")
#m_e1 <- readRDS("detailedtools/RDS/m_e1.rds")

# diagnostics
summary(m_e1)
mcmc_plot(m_e1)
pp_check(m_e1)
plot(conditional_effects(m_e1))

# visualize
# make violin plot
m_type_pred <- m_e1 %>% 
  epred_draws(newdata = tibble(item = detseq_oi$item,
                               Age = detseq_oi$Age,
                               anviltype = detseq_oi$anviltype,
                               subjectID = detseq_oi$subjectID))

# age difference in duration to open item
ggplot(data = m_type_pred, aes(x = Age, y = .epred)) + geom_violin(aes(color = Age, fill = Age), alpha = 0.4) + ylim(0,50) +
  stat_summary(detseq_o, inherit.aes = FALSE, mapping=aes(x = Age, y = seqduration, color = Age), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Age", y = "Seconds required to open item") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14))

# item difference in duration to open item
ggplot(data = m_type_pred, aes(x = item, y = .epred)) + geom_violin(aes(color = item, fill = item), alpha = 0.4) + ylim(0,100) +
  stat_summary(detseq_o, inherit.aes = FALSE, mapping=aes(x = item, y = seqduration, color = item), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Item type", y = "Seconds required to open item") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

#### 2. Number of pounds #####
descdist(detseq_oi$n_pounds)

testdist2.1 <- fitdist(detseq_oi$n_pounds, "pois")
plot(testdist2.1)

# Model 2: Number of pounds depending on age, including item, anviltype, and individual ID as random effect
m_e2 <- brm(n_pounds ~ Age + item*anviltype + (1|subjectID), data = detseq_oi, family = "poisson", iter = 1000, chain = 2, core = 2, backend = "cmdstanr")
# save and load model
# saveRDS(m_e2, "detailedtools/RDS/m_e2.rds")
# m_e2 <- readRDS("detailedtools/RDS/m_e2.rds")

# diagnostics
summary(m_e2)
pp_check(m_e2)
mcmc_plot(m_e2)
plot(conditional_effects(m_e2))

# to test hypotheses
hypothesis(m_e2, "Intercept  > Intercept + anviltypewood", alpha = 0.05)
hypothesis(m_e2, "Intercept  < Intercept + AgeJuvenile", alpha = 0.05)

# visualizing
# make violin plot
m_type_pred2 <- m_e2 %>% 
  epred_draws(newdata = tibble(item = detseq_oi$item,
                               Age = detseq_oi$Age,
                               anviltype = detseq_oi$anviltype,
                               subjectID = detseq_oi$subjectID))

# age difference in number of pounds to open item
ggplot(data = m_type_pred2, aes(x = Age, y = .epred)) + geom_violin(aes(color = Age, fill = Age), alpha = 0.4) +
  stat_summary(detseq_o, inherit.aes = FALSE, mapping=aes(x = Age, y = n_pounds, color = Age), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Age", y = "Number of pounds required to open item") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14))

# item difference in number of pounds to open item
ggplot(data = m_type_pred2, aes(x = item, y = .epred)) + geom_violin(aes(color = item, fill = item), alpha = 0.4) +
  stat_summary(detseq_o[which(detseq_o$item %in% c("almendrabrown", "almendragreen")),], inherit.aes = FALSE, mapping=aes(x = item, y = n_pounds, color = item), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Item type", y = "Number of pounds required to open item") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

#### 2a. Number of pounds per second (with offset of sequence duration) #####
# relationship between nr of pounds and how long the sequence lasts
ggplot(detseq_oi, aes(y = seqduration, x = n_pounds, color = Age, shape = Age)) + geom_point(size = 3) + geom_smooth() + 
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  labs(y = "Seconds needed to open item", x = "Number of pounds needed to open item") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 


# Model 2a: Number of pounds depending on age, including item, anviltype, and individual ID as random effect, and sequence duration as offset
m_e2a <- brm(n_pounds ~ Age + item*anviltype + (1|subjectID) + offset(log(seqduration)), data = detseq_oi, family = "poisson", iter = 1000, chain = 2, core = 2, backend = "cmdstanr")
# save and load model
# saveRDS(m_e2a, "detailedtools/RDS/m_e2a.rds")
# m_e2a <- readRDS("detailedtools/RDS/m_e2a.rds")

# diagnostics
summary(m_e2a)
pp_check(m_e2a)
mcmc_plot(m_e2a)
plot(conditional_effects(m_e2a))

# to test hypotheses
hypothesis(m_e2a, "Intercept  > Intercept + anviltypewood", alpha = 0.05)
hypothesis(m_e2a, "Intercept  < Intercept + AgeJuvenile", alpha = 0.05)

# visualizing
# make violin plot
m_type_pred2 <- m_e2a %>% 
  epred_draws(newdata = tibble(item = detseq_oi$item,
                               Age = detseq_oi$Age,
                               anviltype = detseq_oi$anviltype,
                               subjectID = detseq_oi$subjectID,
                               seqduration = detseq_oi$seqduration))

# age difference in number of pounds to open item
ggplot(data = m_type_pred2, aes(x = Age, y = .epred)) + geom_violin(aes(color = Age, fill = Age), alpha = 0.4) +
  stat_summary(detseq_o, inherit.aes = FALSE, mapping=aes(x = Age, y = n_pounds, color = Age), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Age", y = "Number of pounds required to open item") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14))

#### 3. Number of mistakes #####

# look at item flying, hammer loss and missing separately
# are items more likely to fly on stone or wood?
t.test(detseq_oi$n_flies ~ as.factor(detseq_oi$anviltype))
# seems like flying happens more at stone than wood ( but also coded more there now)

## 3a: True misses (n_miss)
descdist(detseq_oi$n_miss)
testdist3.1 <- fitdist(detseq_oi$n_miss, "pois")
plot(testdist3.1)

# Model 3a: Number of misstrikes (true misses) depending on age, item, anviltype and subjectID as random effect
## ZERO-INFLATED POISSON
m_e3a <- brm(n_miss ~ Age + item*anviltype + (1|subjectID), data = detseq_o, family = zero_inflated_poisson, iter = 1000, chain = 2, core = 2, backend = "cmdstanr", control = list(adapt_delta = 0.99))
# saving and loading model
# saveRDS(m_e3a, "detailedtools/RDS/m_e3a.rds")
# m_e3a <- readRDS("detailedtools/RDS/m_e3a.rds")

# diagnostics
summary(m_e3a)
mcmc_plot(m_e3a)
pp_check(m_e3a)
plot(conditional_effects(m_e3a))

# make violin plot
m_type_pred3 <- m_e3a %>% 
  epred_draws(newdata = tibble(Age = detseq_oi$Age,
                               item = detseq_oi$item,
                               anviltype = detseq_oi$anviltype,
                               subjectID = detseq_oi$subjectID))

# age difference in misstrikes
ggplot(data = m_type_pred3, aes(x = Age, y = .epred)) + geom_boxplot(aes(color = Age, fill = Age), alpha = 0.4) +
  stat_summary(detseq_o, inherit.aes = FALSE, mapping=aes(x = Age, y = n_miss, color = Age), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Age", y = "Average number of mistakes per tool use sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

## 3b: Item flying (itemflies)

# Model 3b: Number of item flies depending on age, item, anviltype and subjectID as random effect
## ZERO-INFLATED POISSON
m_e3b <- brm(n_flies ~ Age + item*anviltype + (1|subjectID), data = detseq_o, family = zero_inflated_poisson, iter = 1000, chain = 2, core = 2, backend = "cmdstanr", control = list(adapt_delta = 0.99))
# saving and loading model
# saveRDS(m_e3b, "detailedtools/RDS/m_e3b.rds")
# m_e3b <- readRDS("detailedtools/RDS/m_e3b.rds")

# diagnostics
summary(m_e3b)
mcmc_plot(m_e3b)
pp_check(m_e3b)
plot(conditional_effects(m_e3b))

# make violin plot
m_type_pred3b <- m_e3b %>% 
  epred_draws(newdata = tibble(Age = detseq_oi$Age,
                               item = detseq_oi$item,
                               anviltype = detseq_oi$anviltype,
                               subjectID = detseq_oi$subjectID))

# age difference in items flying
ggplot(data = m_type_pred3b, aes(x = Age, y = .epred)) + geom_boxplot(aes(color = Age, fill = Age), alpha = 0.4) +
  stat_summary(detseq_o, inherit.aes = FALSE, mapping=aes(x = Age, y = n_flies, color = Age), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Age", y = "Average number of items flying per tool use sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

## 3c: Losing hammer (hammerlost)

# Model 3c: Number of hammer losses depending on age, item, anviltype and subjectID as random effect
## ZERO-INFLATED POISSON
m_e3c <- brm(n_hloss ~ Age + item*anviltype + (1|subjectID), data = detseq_o, family = zero_inflated_poisson, iter = 1000, chain = 2, core = 2, backend = "cmdstanr", control = list(adapt_delta = 0.99))
# saving and loading model
# saveRDS(m_e3c, "detailedtools/RDS/m_e3c.rds")
# m_e3c <- readRDS("detailedtools/RDS/m_e3c.rds")

# diagnostics
summary(m_e3c)
mcmc_plot(m_e3c)
pp_check(m_e3c)
plot(conditional_effects(m_e3c))

# make violin plot
m_type_pred3c <- m_e3c %>% 
  epred_draws(newdata = tibble(Age = detseq_oi$Age,
                               item = detseq_oi$item,
                               anviltype = detseq_oi$anviltype,
                               subjectID = detseq_oi$subjectID))

# age difference in hammer losses
ggplot(data = m_type_pred3c, aes(x = Age, y = .epred)) + geom_boxplot(aes(color = Age, fill = Age), alpha = 0.4) +
  stat_summary(detseq_o, inherit.aes = FALSE, mapping=aes(x = Age, y = n_hloss, color = Age), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Age", y = "Average number of hammer losses per tool use sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

# individual variation in how many mistakes are made
# only take individuals with enough sequences observed
ftable(detseq_oi$subjectID)
knownids <- data.frame(ID = c("TER", "SPT", "BAL", "TOM", "PEA", "SMG", "ZIM", "MIC", "LAR"))
for (i in 1:nrow(knownids)) {
  knownids$nrow[i] <- nrow(detseq_o[which(detseq_o$subjectID == knownids$ID[i]),])
}

detseq_o2 <- left_join(detseq_o[which(detseq_o$subjectID %in% knownids$ID),], knownids, by = c("subjectID" = "ID"))

ggplot(detseq_o2, aes(x=subjectID, y=n_misstotal, color = Age, fill = Age)) + 
  geom_violin(alpha = 0.4) + geom_text(aes(y = 7.5, x = subjectID, label = nrow)) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  labs(x = "Age", y = "Average number of mistakes per tool use sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

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


#### 4. Number of repositions #####

## 4a: repositions of item
# Model 4a: Number of repositions depending on age, item, anviltype and subject ID as random effect
m_e4a <- brm(n_itemreposit ~ Age + item*anviltype + (1|subjectID), data = detseq_oi, family = "poisson", iter = 1000, chain = 2, core = 2, backend = "cmdstanr")
# saving and loading model
# saveRDS(m_e4a, "detailedtools/RDS/m_e4a.rds")
# m_e4a <- readRDS("detailedtools/RDS/m_e4a.rds")

# diagnostics
summary(m_e4a)
mcmc_plot(m_e4a)
pp_check(m_e4a)
plot(conditional_effects(m_e4a))

# make violin plot
m_type_pred4 <- m_e4a %>% 
  epred_draws(newdata = tibble(Age = detseq_oi$Age,
                               item = detseq_oi$item,
                               anviltype = detseq_oi$anviltype,
                               subjectID = detseq_oi$subjectID))

# age difference in number of item repositions
ggplot(data = m_type_pred4, aes(x = Age, y = .epred)) + geom_violin(aes(color = Age, fill = Age), alpha = 0.4) +
  stat_summary(detseq_o, inherit.aes = FALSE, mapping=aes(x = Age, y = n_itemreposit, color = Age), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Age", y = "Average number of repositions per sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

## 4b: repositions of hammers
# Model 4b: Number of hammer repositions depending on age, item, anviltype and subject ID as random effect
m_e4b <- brm(n_hamreposit ~ Age + item*anviltype + (1|subjectID), data = detseq_oi, family = "poisson", iter = 1000, chain = 2, core = 2, backend = "cmdstanr")
# saving and loading model
# saveRDS(m_e4b, "detailedtools/RDS/m_e4b.rds")
# readRDS("detailedtools/RDS/m_e4b.rds")

# diagnostics
summary(m_e4b)
mcmc_plot(m_e4b)
pp_check(m_e4b)
plot(conditional_effects(m_e4b))

# make violin plot
m_type_pred4b <- m_e4b %>% 
  epred_draws(newdata = tibble(Age = detseq_oi$Age,
                               item = detseq_oi$item,
                               anviltype = detseq_oi$anviltype,
                               subjectID = detseq_oi$subjectID))

# age difference in number of hammer repositions
ggplot(data = m_type_pred4b, aes(x = Age, y = .epred)) + geom_violin(aes(color = Age, fill = Age), alpha = 0.4) +
  stat_summary(detseq_o, inherit.aes = FALSE, mapping=aes(x = Age, y = n_hamreposit, color = Age), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Age", y = "Average number of repositions per sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

## 4c: peeling
# Model 4c: Number of peels depending on age, item, anviltype and subject ID as random effect
m_e4c <- brm(n_peel ~ Age + item*anviltype + (1|subjectID), data = detseq_oi, family = "poisson", iter = 1000, chain = 2, core = 2, backend = "cmdstanr")
# saving and loading model
# saveRDS(m_e4c, "detailedtools/RDS/m_e4c.rds")
# m_e4c <- readRDS("detailedtools/RDS/m_e4c.rds")

# diagnostics
summary(m_e4c)
mcmc_plot(m_e4c)
pp_check(m_e4c)
plot(conditional_effects(m_e4c))

# make violin plot
m_type_pred4c <- m_e4c %>% 
  epred_draws(newdata = tibble(Age = detseq_oi$Age,
                               item = detseq_oi$item,
                               anviltype = detseq_oi$anviltype,
                               subjectID = detseq_oi$subjectID))

# age difference in number of peels
ggplot(data = m_type_pred4c, aes(x = Age, y = .epred)) + geom_violin(aes(color = Age, fill = Age), alpha = 0.4) +
  stat_summary(detseq_o, inherit.aes = FALSE, mapping=aes(x = Age, y = n_itemreposit, color = Age), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Age", y = "Average number of peels per sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 


# individual variation

ggplot(detseq_o2, aes(x=subjectID, y=n_reposit, color = Age, fill = Age)) + 
  geom_violin(alpha = 0.4) + geom_text(aes(y = 9.5, x = subjectID, label = nrow)) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  labs(x = "Age", y = "Average number of repositions per tool use sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

### Exploring individual variation and development ####

# focus only on identifiable individuals and data from R11 (which is complete)
detseq_o2c <- detseq_o2[detseq_o2$deployment == "R11",]
ftable(detseq_o2c$item)
ftable(detseq_o2c$subjectID)
# make subjectID factor that is ordered based on age
detseq_o2c$subjectID <- factor(detseq_o2c$subjectID, levels = c("ZIM", "PEA", "BAL", "TER", "MIC", "LAR", "SPT", "TOM", "SMG"))

## What items they process over time
ggplot(detseq_o2c, aes(x = month(mediadate), fill = item)) + geom_histogram() + facet_wrap(~ subjectID) + theme_bw()

# old plot showing change in n_pound, n_miss and n_reposit (but not accounting for itemtype)
ggplot(detseq_o2c[detseq_o2c$location == "EXP-ANV-01",]) + geom_smooth(aes(x = mediadate, y = n_miss, color = "n_miss")) + geom_smooth(aes(x = mediadate, y = n_pounds, color = "n_pounds")) + ylim(0,15) +
  geom_smooth(aes(x = mediadate, y = n_itemreposit,  color = "n_repositions"))  + facet_wrap(~subjectID)  + theme_bw() + scale_color_manual("", breaks = c("n_miss", "n_pounds", "n_repositions"),
                                                                                                                                                       values = c("red", "blue", "green"))

# how n_pounds changes over time, but ideally staying within each itemtype.
ggplot(detseq_o2c[detseq_o2c$item == c("almendrabrown"),]) + geom_point(aes(x = mediadate, y = n_pounds, shape = location), alpha = 0.4, color = "brown", size = 3) + geom_smooth(aes(x = mediadate, y = n_pounds), color = "brown") + 
  facet_wrap(~subjectID) + theme_bw() + ggtitle("Brown Almendra") + labs(x = "Date", y = "Number of pounds to open item") +theme(axis.text = element_text(size = 12),
                                                                                  axis.title = element_text(size = 14)) 

ggplot(detseq_o2c[detseq_o2c$item == c("almendragreen"),]) + geom_point(aes(x = mediadate, y = n_pounds, shape = location), alpha = 0.4, color = "darkgreen", size = 3) + geom_smooth(aes(x = mediadate, y = n_pounds), color = "darkgreen") + 
  facet_wrap(~subjectID) + theme_bw() + ggtitle("Green Almendra") + labs(x = "Date", y = "Number of pounds to open item") +theme(axis.text = element_text(size = 12),
                                                                                                                                 axis.title = element_text(size = 14)) 

# how misstrikes (real misses) change over time, within each itemtype
ggplot(detseq_o2c[detseq_o2c$item == c("almendrabrown"),]) + geom_point(aes(x = mediadate, y = n_miss, shape = location), alpha = 0.4, color = "brown", size = 3) + geom_smooth(aes(x = mediadate, y = n_miss), color = "brown") + 
  facet_wrap(~subjectID) + theme_bw() + ggtitle("Brown Almendra") + labs(x = "Date", y = "Number of mistakes") +theme(axis.text = element_text(size = 12),
                                                                                                                                 axis.title = element_text(size = 14)) + ylim(c(0,10))

# how itemrepositions change over time
ggplot(detseq_o2c[detseq_o2c$item == c("almendrabrown"),]) + geom_point(aes(x = mediadate, y = n_itemreposit, shape = location), alpha = 0.4, color = "brown", size = 3) + geom_smooth(aes(x = mediadate, y = n_itemreposit), color = "brown") + 
  facet_wrap(~subjectID) + theme_bw() + ggtitle("Brown Almendra") + labs(x = "Date", y = "Number of repositions per sequence") +theme(axis.text = element_text(size = 12),
                                                                                                                      axis.title = element_text(size = 14)) + ylim(c(0,10))

# combining n_pounds and n_repositions
ggplot(detseq_o2c[detseq_o2c$item == c("almendrabrown"),]) + geom_point(aes(x = mediadate, y = n_pounds, color = "Pounds"), shape = 16, alpha = 0.4, size = 3) + 
  geom_point(aes(x = mediadate, y = n_itemreposit, color = "Item repositions"), shape = 17, alpha = 0.4, size = 3) + geom_smooth(aes(x = mediadate, y = n_pounds, color = "Pounds")) +  
  geom_smooth(aes(x = mediadate, y = n_itemreposit, color = "Item repositions")) + 
  geom_point(aes(x = mediadate, y = n_miss, color = "Misses"), shape = 18, alpha = 0.4, size = 3) + geom_smooth(aes(x = mediadate, y = n_pounds, color = "Pounds")) +  
  geom_smooth(aes(x = mediadate, y = n_miss, color = "Misses")) +
  facet_wrap(~subjectID, scales = "free_y") + theme_bw() + ggtitle("Brown Almendra") + scale_color_manual("", breaks = c("Pounds", "Misses", "Item repositions"),
                                                                                                          values = c("blue", "darkred", "darkgreen")) +
  labs(x = "Date", y = "Number per sequence") +theme(axis.text = element_text(size = 12),  axis.title = element_text(size = 14)) 

# comparing n_pounds, n_miss, n_reposit for known individuals
melt_detseq <- melt(detseq_o2c, measure.vars = c("n_pounds", "n_misstotal", "n_itemreposit"))

ggplot(melt_detseq) + geom_violin(aes(y = value, x = variable, color = variable, fill = variable), alpha = 0.4) +
  stat_summary(melt_detseq, inherit.aes = FALSE, mapping=aes(x = variable, y = value, color = variable), geom = "point", fun = "mean",
               size = 4) +
  facet_wrap(~subjectID) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  labs(x = "Age", y = "Average number of repositions per sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 


# would probably have to be some kind of GAM, that also includes other things that affects these things (item, anviltype)
# either work with actual number or do maybe a julian day or something? and then split by ID? 
detseq_o2$time <- as.numeric(difftime(detseq_o2$mediadate, min(detseq_o2$mediadate), unit = "days"))
head(detseq_o2c$time)
detseq_o2$subjectID_F <- as.factor(detseq_o2$subjectID)
detseq_o2$Age_f <- as.factor(detseq_o2$Age)

str(detseq_o2c)
dev_gam1 <- gam(n_pounds ~ s(time, by = Age_f) + Age_f + s(subjectID_F, bs = "re"), data = detseq_o2[detseq_o2$item == "almendrabrown",], family = "poisson", method = "REML")
summary(dev_gam1)
draw(dev_gam1)
plot(dev_gam1)

# in brms (trial)
dev_gam1b <- brm(n_pounds ~ s(time, by = Age_f) + Age_f + s(subjectID_F, bs = "re"), data=detseq_o2[detseq_o2$item == "almendrabrown",], family="poisson", 
               chains=2, cores = 4, backend = "cmdstanr", save_pars = save_pars(all = TRUE),
               iter = 1000)
# saving and loading model
# saveRDS(dev_gam1b, "detailedtools/RDS/dev_gam1b.rds")
# dev_gam1b <- readRDS("detailedtools/RDS/dev_gam1b.rds")
plot(conditional_smooths(dev_gam1b))
plot(conditional_effects(dev_gam1b))


### What is ABE doing? ###
# do we see ABE use tools at all in this dataset?
head(dettools_r2)
head(detseq)
ftable(detseq$subjectID)
# we don't see ABE use tools in this dataset at all (INK also does it much less...)

# sidenote: does ABE use tools in agouti data? 
head(agoutisequence_c)
str(agoutisequence_c)
ftable(agoutisequence_c$name)
# make abe subset
ABE_only <- subset(agouticlean, agouticlean$name == "ABE (Abraham)")
ftable(ABE_only$tooluse)
plot(ABE_only$seq_start, ABE_only$tooluse)

# do we see ABE displace/scrounge? 
ABEcomment <- dettools_r2[str_detect(dettools_r2$comment, "ABE|Abraham"),]
# 5 recorded instances of him displacing and scrounging
# all displace
displacements <- dettools_r2[dettools_r2$displacement == "anvildisp" | dettools_r2$displacement == "fulldisp" |
dettools_r2$displacement == "hammerdisp",]

ftable(dettools_r2$displacement)
# how many of displacements is ABE?
ftable(detseq$displacement)
# 96 total displacements
ftable(displacements[!duplicated(displacements$sequenceID),]$scrounging)
ftable(detseq$scrounging)

### Social Attention ####
# how many sequences is someone paying social attention
soc_att <- detseq[detseq$socatt != "None",]
nas <- soc_att[!complete.cases(soc_att),]
# so the NAs are from when they did not open the item so the n_pounds is not there
# filter out when they did not open the item?
ftable(soc_att$outcome[soc_att$socatt == "socialattention"])

# maybe filter out other items than almendras? (because others are so rare) 
# or pool together
soc_att$item2 <- ifelse(str_detect(soc_att$item, "almendra") == FALSE, "other", soc_att$item)
soc_att$attention <- ifelse(soc_att$socatt == "socialattention", 1, 0)

table(soc_att$socatt, soc_att$displacement)

#who is the tool user in these cases
ggplot(data = soc_att[!is.na(soc_att$age_of),], aes(x = age_of, fill = attention)) + geom_histogram(stat = "count")
# what are they processing
ggplot(data = soc_att[str_detect(soc_att$item, "almendra") == TRUE,], aes(x = item, fill = attention)) + geom_histogram(stat = "count")
str(soc_att)
# more likely to pay attention with more efficient ones?
ggplot(data = soc_att, aes(x = n_pounds, fill = attention)) + geom_histogram(stat = "count")
str(soc_att)
ggplot(data = soc_att, aes(x = n_misstotal, fill = attention)) + geom_histogram(stat = "count")
soc_att$age_f <- factor(soc_att$Age, levels = c("Juvenile", "Subadult", "Adult"))
soc_att$hour <- hour(soc_att$videostart)

# make list of videoIDs with social attention
# will need to code all videos with capuchin present what their age/sex are
# also if they don't pay social attention
ftable(soc_att$socatt)
# probably easiest is to code all these videos in a separate BORIS with separate ethogram
# but with same video ID, and then left_join the information with the full dataset here
# cause otherwise I'd need to go into the BORIS files of Meredith & Leonie which is very hard
# ones I coded
socatt_vidnames <- soc_att[,c("videoID", "coder", "subjectID", "attention", "scrounging", "displacement")]
# filter to ones not coded yet
# load in coding
socatt_c <- read.csv("detailedtools/socialattentioncoding.csv")
tocode <- socatt_vidnames[!socatt_vidnames$videoID %in% socatt_c$Observation.id,]
tocode[0:nrow(tocode),]
# still need to do meredith cebus-02 and last bit of leonie coding 

## NOTE
# if someone was displaced then there was someone else present (the displacer)
# so check if we had social attention 0 but displacing yes, then need to fix that
# make sure to override if someone coded social attention 0 but I coded 1 in the detailed coding (or have them fix it in their BORIS file)

## binomial model predicting whether or not they pay social attention
# first just generally, are they more likely to pay more attention to:
# certain age classes, certain items being opened, depending on outcome, depending on location, offset of duration and random effect of who tool user was
# with item in
soc_att_bm1 <- brm(attention ~ age_f + item2 + outcome + offset(log(seqduration)) + (1|subjectID) + location, data = soc_att, family = bernoulli(),
                   iter = 1000, chains=2, cores = 4, backend = "cmdstanr", save_pars = save_pars(all = TRUE))
# saving and loading model
# saveRDS(soc_att_bm1, "detailedtools/RDS/soc_att_bm1.rds")
# soc_att_bm1 <- readRDS("detailedtools/RDS/soc_att_bm1.rds")

summary(soc_att_bm1)
mcmc_plot(soc_att_bm1)
plot(conditional_effects(soc_att_bm1))
plot(soc_att_bm1)

## NOTE: if we want to include subject ID as random effect, need to exclude the unknowns (or ID them!)

# then only for when items were opened successfully and we know how efficient the tool user was (n_pounds, n_miss). 
# more likely to pay attention to more efficient tool users? 
soc_att_bm1b <- brm(attention ~ age_f + n_pounds + item2 +  n_misstotal + offset(log(seqduration)) + (1|subjectID) + location, data = soc_att[soc_att$outcome == "opened",], family = bernoulli(),
                   iter = 1000, chains=2, cores = 4, backend = "cmdstanr", save_pars = save_pars(all = TRUE))
# saving and loading model
# saveRDS(soc_att_bm1b, "detailedtools/RDS/soc_att_bm1b.rds")
# soc_att_bm1b <- readRDS("detailedtools/RDS/soc_att_bm1b.rds")

summary(soc_att_bm1b)
mcmc_plot(soc_att_bm1b)
plot(conditional_effects(soc_att_bm1b))

hypothesis(soc_att_bm1b, "Intercept  > Intercept + n_pounds", alpha = 0.05)

# consider if interactions are interesting 
# e.g. number of pounds * age of tool user (so more atttention to efficient subadult than adult?)
soc_att_bm1c <- brm(attention ~ age_f* n_pounds + item2 +  n_misstotal + offset(log(seqduration)) + (1|subjectID) + location, data = soc_att[soc_att$outcome == "opened",], family = bernoulli(),
                    iter = 1000, chains=2, cores = 4, backend = "cmdstanr", save_pars = save_pars(all = TRUE))
# saving and loading model
# saveRDS(soc_att_bm1c, "detailedtools/RDS/soc_att_bm1c.rds")
# soc_att_bm1c <- readRDS("detailedtools/RDS/soc_att_bm1c.rds")

summary(soc_att_bm1c)
mcmc_plot(soc_att_bm1c)
plot(conditional_effects(soc_att_bm1c))

hypothesis(soc_att_bm1b, "Intercept  > Intercept + n_pounds", alpha = 0.05)

# consider doing a GAM with social attention depending on hour of the day?
# might be non-linear relationship?

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
ggplot(detseq_oh[detseq_oh$deployment == "R11" & detseq_oh$location == "EXP-ANV-01",], 
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
