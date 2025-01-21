## Detailed tool use -- Analyses
## MPI-AB; Z Goldsborough

## Analysing efficiency of tool use and social attention paid to tool use

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
library(easystats)

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

###
### EFFICIENCY ####
###

#### Filter and general diagnostics ####
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

ggplot(detseq, 
       aes(x = mediadate, fill = location)) + geom_histogram() + theme_bw()
ftable(detseq$socatt)

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

#### Analyses ####
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
hypothesis(m_e1, "Intercept < Intercept + itemalmendragreen", alpha = 0.05)

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
            iter = 2000, chain = 3, core = 3, save_pars = save_pars(all = TRUE), 
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
hypothesis(m_e1, "Intercept  < Intercept + itemalmendragreen", alpha = 0.05)
hypothesis(m_e2, "Intercept > Intercept + AgeAdult", alpha = 0.05)

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
# png("detailedtools/RDS/duration_pound.png", width = 8, height = 7, units = 'in', res = 300)
ggplot(detseq_oi, aes(y = seqduration, x = n_pounds, color = Age, shape = Age)) + geom_point(size = 3, alpha = 0.3) + geom_smooth() + 
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  labs(y = "Seconds needed to open item", x = "Number of pounds needed to open item") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 
# dev.off()

##### 3. Number of mistakes ######
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

## 3a: True misses (n_miss)
descdist(detseq_oi$n_miss)
testdist3.1 <- fitdist(detseq_oi$n_miss, "pois")
plot(testdist3.1)
# use a zero-inflated poisson

### Model_e3a ### 
# Outcome: number of misstrikes (true misses)
# Fixed effects: age, item, and anviltype
# Random effects: subjectID
m_e3a <- brm(n_miss ~ Age + item + anviltype + (1|subjectID), data = detseq_oi, family = zero_inflated_poisson, 
             iter = 2000, chain = 2, core = 2, save_pars = save_pars(all = TRUE), backend = "cmdstanr", control = list(adapt_delta = 0.99))
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
loo_R2(m_e3a) # 0.17
round(bayes_R2(m_e3a),2) # 0.16

# Interpretation
round(exp(-0.77-3.63),2) * 0.48
hypothesis(m_e3a, "Intercept > Intercept + AgeAdult", alpha = 0.05)

# make violin plot
m_type_pred3 <- m_e3a %>% 
  epred_draws(newdata = tibble(Age = detseq_oi$Age,
                               item = detseq_oi$item,
                               anviltype = detseq_oi$anviltype,
                               subjectID = detseq_oi$subjectID))

# age difference in misstrikes
ggplot(data = m_type_pred3, aes(x = Age, y = .epred)) + geom_boxplot(aes(color = Age, fill = Age), alpha = 0.4) +
  stat_summary(detseq_oi, inherit.aes = FALSE, mapping=aes(x = Age, y = n_miss, color = Age), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Age", y = "Average number of mistakes per tool use sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

## 3b: Item flying (itemflies)

### Model_e3b ### 
# Outcome: number of itemflies (item flying off)
# Fixed effects: age, item, and anviltype
# Random effects: subjectID
m_e3b <- brm(n_flies ~ Age + item + anviltype + (1|subjectID), data = detseq_oi, 
             family = zero_inflated_poisson, iter = 2000, chain = 2, core = 2, 
             save_pars = save_pars(all = TRUE),
             backend = "cmdstanr", control = list(adapt_delta = 0.99))
# m_e3b <- add_criterion(m_e3b, c("loo", "loo_R2", "bayes_R2"), reloo = TRUE, backend = "cmdstanr", ndraws = 2000) 

# saving and loading model
# saveRDS(m_e3b, "detailedtools/RDS/m_e3b.rds")
# m_e3b <- readRDS("detailedtools/RDS/m_e3b.rds")

# diagnostics
summary(m_e3b)
mcmc_plot(m_e3b)
pp_check(m_e3b)
plot(conditional_effects(m_e3b))

loo(m_e3b) # all cases good
loo_R2(m_e3b) # 0.10
round(bayes_R2(m_e3b),2) # 0.13

# Interpretation
round(exp(-1.08),2) * 0.40
hypothesis(m_e3b, "Intercept > Intercept + AgeAdult", alpha = 0.05)

# make violin plot
m_type_pred3b <- m_e3b %>% 
  epred_draws(newdata = tibble(Age = detseq_oi$Age,
                               item = detseq_oi$item,
                               anviltype = detseq_oi$anviltype,
                               subjectID = detseq_oi$subjectID))

# age difference in items flying
ggplot(data = m_type_pred3b, aes(x = Age, y = .epred)) + geom_boxplot(aes(color = Age, fill = Age), alpha = 0.4) +
  stat_summary(detseq_oi, inherit.aes = FALSE, mapping=aes(x = Age, y = n_flies, color = Age), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Age", y = "Average number of items flying per tool use sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

## 3c: Losing hammer (hammerlost)
ftable(detseq_oi$n_hloss, detseq_oi$item)
# Model 3c: Number of hammer losses depending on age
# subjectID as random effect
## ZERO-INFLATED POISSON
m_e3c <- brm(n_hloss ~ Age + (1|subjectID), data = detseq_oi, 
             family = zero_inflated_poisson, iter = 2000, chain = 2, core = 2, 
             save_pars = save_pars(all = TRUE), backend = "cmdstanr", 
             control = list(adapt_delta = 0.99))
 m_e3c <- add_criterion(m_e3c, c("loo", "loo_R2", "bayes_R2"), reloo = TRUE, backend = "cmdstanr", ndraws = 2000) 

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
  stat_summary(detseq_oi, inherit.aes = FALSE, mapping=aes(x = Age, y = n_hloss, color = Age), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Age", y = "Average number of hammer losses per tool use sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

## 3d: all mistakes together

# Model 3d: Number of hammer losses depending on age, item, anviltype and subjectID as random effect
## ZERO-INFLATED POISSON
m_e3d <- brm(n_misstotal ~ Age + item + anviltype + (1|subjectID), data = detseq_oi, 
             family = zero_inflated_poisson, iter = 2000, chain = 2, core = 2, 
             save_pars = save_pars(all = TRUE), backend = "cmdstanr", 
             control = list(adapt_delta = 0.99))
# m_e3d <- add_criterion(m_e3d, c("loo", "loo_R2", "bayes_R2"), reloo = TRUE, backend = "cmdstanr", ndraws = 2000) 

# saving and loading model
# saveRDS(m_e3d, "detailedtools/RDS/m_e3d.rds")
# m_e3d <- readRDS("detailedtools/RDS/m_e3d.rds")




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
m_e4a <- brm(n_itemreposit ~ Age + item*anviltype + (1|subjectID), data = detseq_oi, family = "poisson", 
             iter = 1000, chain = 2, core = 2, save_pars = save_pars(all = TRUE), backend = "cmdstanr")
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
  stat_summary(detseq_oi, inherit.aes = FALSE, mapping=aes(x = Age, y = n_hamreposit, color = Age), geom = "point", fun = "mean",
               size = 4) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  guides(color = "none", fill = "none") +
  labs(x = "Age", y = "Average number of repositions per sequence") +
  theme_bw() + theme(axis.text = element_text(size = 12),
                     axis.title = element_text(size = 14)) 

## 4c: peeling
# Model 4c: Number of peels depending on age, item, anviltype and subject ID as random effect
m_e4c <- brm(n_peel ~ Age + item*anviltype + (1|subjectID), data = detseq_oi,
             save_pars = save_pars(all = TRUE), family = "poisson", iter = 1000, chain = 2, core = 2, backend = "cmdstanr")
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

# focus only on identifiable individuals who had more than 10 observations (so no Joe and Ink)
detseq_o2c <- detseq_o2
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
detseq_o2c$time <- as.numeric(difftime(detseq_o2c$mediadate, min(detseq_o2c$mediadate), unit = "days"))
head(detseq_o2c$time)
detseq_o2c$subjectID_F <- as.factor(detseq_o2c$subjectID)
detseq_o2c$Age_f <- as.factor(detseq_o2c$Age)

str(detseq_o2c)
dev_gam1 <- gam(n_pounds ~ s(time, by = Age_f) + Age_f + s(subjectID_F, bs = "re"), data = detseq_o2c[detseq_o2c$item == "almendrabrown",], family = "poisson", method = "REML")
summary(dev_gam1)
draw(dev_gam1)
plot(dev_gam1)

# in brms (trial)
dev_gam1b <- brm(n_pounds ~ s(time, by = Age_f) + Age_f + s(subjectID_F, bs = "re"), data=detseq_o2c[detseq_o2c$item == "almendrabrown",], family="poisson", 
               chains=2, cores = 2, backend = "cmdstanr", save_pars = save_pars(all = TRUE),
               iter = 1000)
# saving and loading model
# saveRDS(dev_gam1b, "detailedtools/RDS/dev_gam1b.rds")
# dev_gam1b <- readRDS("detailedtools/RDS/dev_gam1b.rds")
plot(conditional_smooths(dev_gam1b))
plot(conditional_effects(dev_gam1b))

dev_gam2b <- brm(n_pounds ~ s(time, by = subjectID_F) + subjectID_F + Age_f, data=detseq_o2c[detseq_o2c$item == "almendrabrown",], family="poisson", 
                 chains=2, cores = 2, backend = "cmdstanr", save_pars = save_pars(all = TRUE),
                 iter = 1000)
# saving and loading model
# saveRDS(dev_gam2b, "detailedtools/RDS/dev_gam2b.rds")
# dev_gam2b <- readRDS("detailedtools/RDS/dev_gam2b.rds")
plot(conditional_smooths(dev_gam2b))
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
## for coding
# ones I coded
socatt_vidnames <- soc_att[,c("videoID", "coder", "subjectID", "attention", "scrounging", "displacement")]
# filter to ones not coded yet
# load in coding
socatt_c <- read.csv("detailedtools/socialattentioncoding.csv")
tocode <- socatt_vidnames[!socatt_vidnames$videoID %in% socatt_c$Observation.id,]
tocode[0:nrow(tocode),]

## first look at the detailed social attention coding
head(socatt_seq)
head(socatt_ct)

# my first hunch says to use the socatt_ct, filtered just to presence
head(socatt_ct)
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
# 990 sequences with capuchins present
# now attach the relevant information from the main dataframe (information on the tool user etc)
head(detseq)

socatt_final <- left_join(socatt_final, detseq[,c("sequenceID", "subjectID", "coder", "location", "item",
                                                    "outcome", "displacement", "scrounging",
                                                   "anviltype", "seqduration", "n_pounds", "n_misstotal",
                                                   "Age", "deployment", "split", "hammerswitches", "anvilswitches")],
                          by = "sequenceID")

# also add more detailed information on number of capuchins present, number scrounging etc, from the social coding
head(socatt_seq)
socatt_final <- left_join(socatt_final, socatt_seq[,c("sequenceID", "n_socatt", "n_disp", "n_scr", "p_total")],
                          by = "sequenceID")

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
# then have 837 sequences
length(unique(socatt_final$sequenceID))

# in how many sequences do we see social attention?
ftable(socatt_final$socialattention[!duplicated(socatt_final$sequenceID)])

# what items are being processed
ftable(socatt_final$item[!duplicated(socatt_final$sequenceID)])
# maybe pool together items than almendras? (because others are so rare) 
# or pool together
socatt_final$item2 <- ifelse(str_detect(socatt_final$item, "almendra") == FALSE, "other", socatt_final$item)
ftable(socatt_final$item2[!duplicated(socatt_final$sequenceID)])

ftable(socatt_final$tooluserID)
ftable(socatt_final$observerID)

# check for NAs
which(is.na(socatt_final), arr.ind = TRUE)

# so we usually have the ID of the tool user, not always of the observer
# could envision doing a model with ID in subsetted to a more thorough dataset
# but I think this will be more for descriptives and for the model we do more general on age
# probably need random effect of sequenceID?


# want to end up with dataframe that's like
# sequence ID, tool user ID, tool user age sex, present age sex, social attention 1/0, total_n capuchins (excl tool user, incl observer),
# item processed, outcome of tool use event, sequence duration, location, n_scrounging (?), n_displacement (?)

#

## binomial model predicting whether or not they pay social attention
# first just generally, is there more likely to be social attention for
# certain age classes (tool user), certain age classes (observer), at certain locations, when more capuchins are present
# in combination with scrounging? Offset of sequence duration (opportunity) and inclusion of random effect (sequenceID)
socatt_bm1 <- brm(socatt ~ tooluser_age + observer_agesex + offset(log(seqduration)) + location + p_total + n_scr + (1|sequenceID),
                   data = socatt_final, family = bernoulli(),
                   iter = 2000, chains=2, cores = 2, backend = "cmdstanr", save_pars = save_pars(all = TRUE))
# saving and loading model
# saveRDS(socatt_bm1, "detailedtools/RDS/socatt_bm1.rds")
# socatt_bm1 <- readRDS("detailedtools/RDS/socatt_bm1.rds")

summary(socatt_bm1)
mcmc_plot(socatt_bm1)
plot(conditional_effects(socatt_bm1))
plot(socatt_bm1)

## alternatively, differentiate known and unknown juveniles (so known tool users, and likely not tool users)
socatt_final$observer_agesex2 <- socatt_final$observer_agesex
socatt_final$observer_agesex2[which(socatt_final$observer_agesex == "juvenile" & 
                                      socatt_final$observerID %in% knownids$ID)] <- "juveniletool"
table(socatt_final$observer_agesex2, socatt_final$socatt)
## I THINK THIS DOESN'T WORK IN A FORMAL ANALYSIS
# because we do reliably recognize juveniles when paying social attention, but not when they are just present
# so I still think this exploration should be more descriptive (just like who receives more social attention)

socatt_bm1j <- brm(socatt ~ tooluser_age + observer_agesex2 + offset(log(seqduration)) + location + p_total + n_scr + (1|sequenceID),
                  data = socatt_final, family = bernoulli(),
                  iter = 2000, chains=2, cores = 2, backend = "cmdstanr", save_pars = save_pars(all = TRUE))
# saving and loading model
# saveRDS(socatt_bm1j, "detailedtools/RDS/socatt_bm1j.rds")
# socatt_bm1j <- readRDS("detailedtools/RDS/socatt_bm1j.rds")

summary(socatt_bm1j)
mcmc_plot(socatt_bm1j)
plot(conditional_effects(socatt_bm1j))
plot(socatt_bm1j)


## NOTE: if we want to include subject ID as random effect, need to exclude the unknowns (or ID them!)
## NOTE: if we'd want to look at outcome or item, we don't have that much variation in the categories. Are those worth it?
# I don't think they are super interesting per se... Is there harm to including them? maybe as random effects?
# Interaction of agesex tool user and observer sounds interesting, but a lot of combinations do not occur so it causes convergence issues
# maybe explore visually rather than in a model? 

## MODEL 2 
# depending on the efficiency of the tool user?
# need to subset to when they were opened successfully

# then only for when items were opened successfully and we know how efficient the tool user was (n_pounds, n_miss). 
# more likely to pay attention to more efficient tool users? 
socatt_bm1b <- brm(socatt ~ tooluser_age + observer_agesex  + n_pounds + n_misstotal + offset(log(seqduration)), 
                   data = socatt_final[socatt_final$outcome == "opened",], family = bernoulli(),
                   iter = 2000, chains=2, cores = 4, backend = "cmdstanr", save_pars = save_pars(all = TRUE))
# saving and loading model
# saveRDS(socatt_bm1b, "detailedtools/RDS/socatt_bm1b.rds")
# socatt_bm1b <- readRDS("detailedtools/RDS/socatt_bm1b.rds")

summary(socatt_bm1b)
mcmc_plot(socatt_bm1b)
plot(conditional_effects(socatt_bm1b))

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

### MODEL 3
# maybe looking at, when the ID of the tool user is known, whether there are specific individuals that tolerate a lot of social attention?
# filter to known tool user IDs
# then have subjectID in (not as random effect but as fixed effect?)
socatt_bm2 <- brm(socatt ~ tooluser_age + observer_agesex  + tooluserID + offset(log(seqduration)), 
                   data = socatt_final[which(socatt_final$tooluserID %in% knownids$ID),], family = bernoulli(),
                   iter = 2000, chains=2, cores = 2, backend = "cmdstanr", save_pars = save_pars(all = TRUE))
# saving and loading model
# saveRDS(socatt_bm2, "detailedtools/RDS/socatt_bm2.rds")
# socatt_bm2 <- readRDS("detailedtools/RDS/socatt_bm2.rds")

summary(socatt_bm1b)
mcmc_plot(socatt_bm1b)
plot(conditional_effects(socatt_bm2))

table(socatt_final$tooluserID[!duplicated(socatt_final$sequenceID)],
                              socatt_final$socialattention[!duplicated(socatt_final$sequenceID)])
str(ind_socatt)
ind_socatt$percentage <- (ind_socatt$Var2/(ind_socatt$Var1+ind_socatt$Var2)) * 100

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

