head(agoutisequence_c)

### number of capuchins inspecting long format

wide_inspect <- agoutisequence_c[,c("sequenceID", "n_inspect_adult", "n_inspect_nonadult", "n_inspect_unknown")]
long_inspect <- melt(wide_inspect, measure.vars = c("n_inspect_adult", "n_inspect_nonadult", "n_inspect_unknown"))
colnames(long_inspect) <- c("sequenceID", "Age", "n_inspect_melted")
long_inspect$Age <- ifelse(long_inspect$Age == "n_inspect_adult", "adult", 
                           ifelse(long_inspect$Age == "n_inspect_nonadult", "nonadult", "unknown"))

## number of capuchins present per age class
agoutisequence_c$nNonadult <- agoutisequence_c$nJuvenile + agoutisequence_c$nSubadult

agoutisequence_cm <- melt(agoutisequence_c, measure.vars = c("nAdult", "nNonadult", "nr_unknown"))
agoutisequence_cm$nr_capuchins <- agoutisequence_cm$value 
agoutisequence_cm$Age <- agoutisequence_cm$variable
agoutisequence_cm$Age <- ifelse(agoutisequence_cm$Age == "n_inspect_adult", "adult", 
                           ifelse(agoutisequence_cm$Age == "n_inspect_nonadult", "nonadult", "unknown"))
head(agoutisequence_cm)

head(long_inspect)
long_inspect$identifier <- paste(long_inspect$sequenceID, long_inspect$Age)
agoutisequence_cm$identifier <- paste(agoutisequence_cm$sequenceID, agoutisequence_cm$Age)

agoutisequence_cm <- left_join(agoutisequence_cm, long_inspect[,c("identifier", "n_inspect_melted")], by = "identifier")
ftable(long_inspect$Age)

cap_inspect <- agoutisequence_cm[agoutisequence_cm$n > 0,c("deploymentID", "sequenceID", "locationName", "tags", "seq_start", "seq_end",
                                                           "seq_length", "mediatype", "uniqueloctag", "dep_start", "dep_end", "dep_length_hours", 
                                                           "month", "season", "island", "tool_anvil", "tool_site", "streambed", "n", "nAdult", "nNonadult",
                                                           "nr_unknown", "n_inspect", "toolusers", "depdays", "depnr", "nr_capuchins", "Age", "n_inspect_melted")]

# on level of the sequence
# inspections: n_inspect is total number of capuchins inspecting in that sequence (repeated 3 times, no idea of age) --> linked to unique sequence ID
# n_inspect_melted is inspecting, by age (not repeated, each line unique) --> linked to the unique combination of sequence ID and age
# nr of capuchins present: n is total number of capuchins irrespective of age (repeated 3 times) --> linked to the unique sequence ID
# nr_capuchins is number of capuchins present by age (so number of adults, number of nonadults) --> linked to unique combination of sequence ID and age

# in model saying n_inspect_melted with offset of nr_capuchins means you are saying:"what proportion of each age class inspects". So if there are 3 adults, how many are inspecting. 
# if you want to see how many inspect depending on how many are present OVERALL then you would n_inspect_melted but with n as an offset. Then question is, if there are 5 capuchins, how many adults/nonadults are inspecting
# before running model would subset on only adults and nonadults (so get rid of the unknowns)
