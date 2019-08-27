# Taking a first look at the Fungal Competition Isotope Data.

setwd("~/Documents/2018-2019/Fungal competition/")

require(tidyverse)

R13C = 0.0112372 # Slater et al. 2011
R15N = 0.0036765 # Slater et al. 2011

# fcdata = read_csv("FC_isotope_data_from_UCSC.csv") # to revert to old analysis, uncomment
fcdata = read_csv("FC_isotope_data_from_UCSC_July.csv")
metadata = read_csv("Fungal_competition_plant_tracking.csv")

head(metadata)
colnames(fcdata)

fcdata = select(fcdata, Plant = plant, Side = side, everything())
metadata$Side = tolower(metadata$Side)
fcdata$Side = tolower(fcdata$Side)

# metadata[is.na(metadata$Batch),]

isowithmeta = left_join(fcdata, metadata)
isowithmeta$enriched = numeric(nrow(isowithmeta))
isowithmeta$enriched[grep("FC", isowithmeta$Batch)] = 1

summary(as.factor(isowithmeta$tissue))
isowithmeta$tissue = tolower(isowithmeta$tissue)
summary(as.factor(isowithmeta$tissue))
isowithmeta$tissue[isowithmeta$tissue == "hyphae/rhizomorphs"] = "hyphae"
isowithmeta$tissue[grep("mycos", isowithmeta$tissue)] = "mycorrhizas"
isowithmeta$tissue[grep("nm|tfr", isowithmeta$tissue)] = "uncolonized_roots"

mixed_compartments_need_attention = c("6002b", "6057a", "6067b", "6098a")
#6002b: Identifier 2 = 3_F10, it's SUIPU, based on mass and notes
#6002b: Identifier 2 = 4_C3, it's THETE, based on mass and notes
#6057a: ID 2 = 3_B12, it's SUIPU
#6057a: ID2 = 2_C6, it's THETE
# 6067b: ID2 = 4_E9, it's SUIPU
# 6067b: ID2 = 4_H4, it's THETE
# 6098a: Apparently only ran one mycorrhiza sample, and it failed. Boo.

SUIPU_ids = c("3_F10", "3_B12", "4_E9")
THETE_ids = c("4_C3", "2_C6", "4_H4")

isowithmeta$mycofungus = isowithmeta$Actual_fungus_by_compartment

# Let's straighten out those mixed compartments

isowithmeta$mycofungus[isowithmeta$`Identifier 2` %in% THETE_ids] = "THETE"

isowithmeta$mycofungus[isowithmeta$`Identifier 2` %in% SUIPU_ids] = "SUIPU"

isowithmeta$uniqueIds = paste(isowithmeta$Plant, isowithmeta$Side, isowithmeta$tissue, isowithmeta$mycofungus)

isowithmeta$replicated = duplicated(isowithmeta$uniqueIds)

isowithmeta[isowithmeta$replicated == TRUE,] # 15 of these.
# Some were TFRs (terminal fine roots) that I recoded as 
# uncolonized roots; others appear to just have been samples
# that I genuinely replicated for some reason.

isowithmeta = subset(isowithmeta, is.na(isowithmeta$d13C) == FALSE)

# Let's average out the replicates.

for (i in 1:nrow(isowithmeta)) {
  for (j in 1:nrow(isowithmeta)) {
    if (isowithmeta$uniqueIds[i] == isowithmeta$uniqueIds[j]) {
      isowithmeta$d13C[i] = mean(isowithmeta$d13C[i], isowithmeta$d13C[j])
      isowithmeta$pctC[i] = mean(isowithmeta$pctC[i], isowithmeta$pctC[j])
      isowithmeta$mcgC[i] = mean(isowithmeta$mcgC[i], isowithmeta$mcgC[j])
      isowithmeta$d15N[i] = mean(isowithmeta$d15N[i], isowithmeta$d15N[j])
      isowithmeta$pctN[i] = mean(isowithmeta$pctN[i], isowithmeta$pctN[j])
      isowithmeta$mcgN[i] = mean(isowithmeta$mcgN[i], isowithmeta$mcgN[j])
    }
  }
}

isowithmeta = subset(isowithmeta, replicated == FALSE)

atompct = function(delta, R) {
  output = 100/(1+(1/(((delta/1000)+1)*R)))
  return(output)
}

isowithmeta$atmpct13C = numeric(nrow(isowithmeta))
isowithmeta$atmpct15N = numeric(nrow(isowithmeta))

for (i in 1:nrow(isowithmeta)) {
  isowithmeta$atmpct13C[i] = atompct(isowithmeta$d13C[i], R13C)
  isowithmeta$atmpct15N[i] = atompct(isowithmeta$d15N[i], R15N)
}

baseline = isowithmeta[isowithmeta$enriched == 0 & isowithmeta$tissue != "hyphae",]


meanbaseline = baseline %>% summarize(meand13C = mean(d13C), 
                                      meand15N = mean(d15N),
                                      meanatmpct13C = mean(atmpct13C),
                                      meanatmpct15N = mean(atmpct15N))

# isowithmeta$mmol13Cexcess = numeric(nrow(isowithmeta))
# isowithmeta$mmol15Nexcess = numeric(nrow(isowithmeta))
# 
# for (i in 1:nrow(isowithmeta)) {
#   deltadiff = (isowithmeta$d13C[i] - meanbaseline$meand13C)
#   isowithmeta$mmol13Cexcess[i] = R13C*deltadiff
#   ndeltadiff = (isowithmeta$d15N[i] - meanbaseline$meand15N)
#   isowithmeta$mmol15Nexcess[i] = R15N*ndeltadiff
# }

# Well, at least I've got a couple of representatives of each host
# in this data set.

isowithmeta$mmol13Cexcess = numeric(nrow(isowithmeta))
isowithmeta$mmol15Nexcess = numeric(nrow(isowithmeta))

isowithmeta$APE13C = numeric(nrow(isowithmeta))
isowithmeta$APE15N = numeric(nrow(isowithmeta))

for (i in 1:nrow(isowithmeta)) {
  deltadiff = (isowithmeta$d13C[i] - meanbaseline$meand13C)
  isowithmeta$mmol13Cexcess[i] = R13C*deltadiff
  
  ndeltadiff = (isowithmeta$d15N[i] - meanbaseline$meand15N)
  isowithmeta$mmol15Nexcess[i] = R15N*ndeltadiff
  
  isowithmeta$APE13C[i] = isowithmeta$atmpct13C[i] - meanbaseline$meanatmpct13C
  isowithmeta$APE15N[i] = isowithmeta$atmpct15N[i] - meanbaseline$meanatmpct15N
  
}

# write_csv(isowithmeta, "Cleaned_processed_FC_isotope_data.csv")
write_csv(isowithmeta, "Cleaned_processed_FC_isotope_data_July.csv")

#### EVERYTHING FROM HERE DOWN IS NOT USED IN DOWNSTREAM PIPELINE ####

isowithmeta$enriched = as.factor(isowithmeta$enriched)

# ggplot(data = isowithmeta) +
#   geom_boxplot(aes(x = enriched, y = mmol13Cexcess, color = N_level))
# ggplot(data = subset(isowithmeta, Actual_fungi_at_harvest != "FAILED")) +
#   geom_boxplot(aes(x = Actual_fungi_at_harvest, y = APE15N, color = N_level))
# 
# forplot = subset(isowithmeta, Actual_fungi_at_harvest != "FAILED" & enriched == 1)
# forplot = subset(isowithmeta, is.na(APE15N) == FALSE)
# minAP15N = min(forplot$APE15N) + 0.00001
# forplot$APE15N = forplot$APE15N - minAP15N
# ggplot(data = forplot) +
#   # geom_boxplot(aes(x = Actual_fungi_at_harvest, y = (APE15N), color = `Receives 15N label?`)) +
#   geom_jitter(width = 0.25, aes(x = Actual_fungi_at_harvest, y = log(APE15N), 
#                   color = `Receives 15N label?`, shape = tissue)) +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1))
# Okay, well, this is disappointing, but I think
# 1) I need to run more of my unenriched control
# samples, and 2) It is possible my panels leaked,
# because NM/NM plants aren't too different from the others.
# ggplot(data = subset(forplot, tissue == "mycorrhizas" & 
#                        (Actual_fungus_by_compartment == "THETE" | Actual_fungus_by_compartment == "SUIPU"))) +
#   # geom_boxplot(aes(x = Actual_fungi_at_harvest, y = (APE15N), color = `Receives 15N label?`)) +
#   geom_jitter(width = 0.25, aes(x = Actual_fungi_at_harvest, y = (APE15N), 
#                                 color = `Receives 15N label?`, 
#                                 shape = Actual_fungus_by_compartment)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
#


mynms = isowithmeta[isowithmeta$Actual_fungi_at_harvest == "NM/NM",]

# How well does rhizomorph 15N & 13C? correlate to mycorrhizas?
rhizosvshyphae_carbon = select(isowithmeta, Plant, Side, tissue, APE13C)
rhizosvshyphae_nitrogen = select(isowithmeta, Plant, Side, tissue, APE15N)

# not working
rhizosvshyphae_carbon %>% spread(APE13C, tissue)

ggplot(data = isowithmeta) +
  geom_boxplot(aes(x = enriched, y = mmol15Nexcess, color = N_level))


carbon = select(isowithmeta, Plant, Side, `Receives 15N label?`, Actual_fungus_by_compartment, Actual_fungi_at_harvest, tissue, mycofungus, mmol13Cexcess)
carbon = carbon[!is.na(carbon$`Receives 15N label?`),]
# small$uniqueId = c(1:nrow(small))
carbonspread = carbon %>% spread(tissue, mmol13Cexcess)

carbonspread = select(carbonspread, everything(), hyphaeC = hyphae, mycosC = mycorrhizas, uncolC = uncolonized_roots)

nitrogen = select(isowithmeta, Plant, Side, `Receives 15N label?`, Actual_fungus_by_compartment, Actual_fungi_at_harvest, tissue, mycofungus, mmol15Nexcess)
nitrogen = nitrogen[!is.na(nitrogen$`Receives 15N label?`),]
nitrogenspread = nitrogen %>% spread(tissue, mmol15Nexcess)
nitrogenspread = select(nitrogenspread, everything(), hyphaeN = hyphae, mycosN = mycorrhizas, uncolN = uncolonized_roots)

pctNtable = select(isowithmeta, Plant, Side, `Receives 15N label?`, Actual_fungus_by_compartment, Actual_fungi_at_harvest, tissue, mycofungus, pctN)
pctNtable = pctNtable[!is.na(pctNtable$`Receives 15N label?`),]
pctNspread = pctNtable %>% spread(tissue, pctN)

everything = full_join(carbonspread, nitrogenspread)
everything = left_join(everything, pctNspread)

everything$hyphaeCratio = numeric(nrow(everything))
everything$mycosCratio = numeric(nrow(everything))
everything$nmCratio = numeric(nrow(everything))
everything$hyphaeNratio = numeric(nrow(everything))
everything$mycosNratio = numeric(nrow(everything))
everything$nmNratio = numeric(nrow(everything))
everything$absorpCratio = numeric(nrow(everything))
everything$pctNratiomycos = numeric(nrow(everything))
everything$pctNratiouncol = numeric(nrow(everything))


for (i in 1:nrow(everything)) {
  for (j in 1:nrow(everything)) {
    if (everything$Plant[i] == everything$Plant[j] &
        everything$Side[i] != everything$Side[j]) {
      if (everything$`Receives 15N label?`[i] == "Y") {
        everything$hyphaeCratio[i] = everything$hyphaeC[i]/everything$hyphaeC[j]
        everything$mycosCratio[i] = everything$mycosC[i]/everything$mycosC[j]
        everything$nmCratio[i] = everything$uncolC[i]/everything$uncolC[j]
        everything$hyphaeNratio[i] = everything$hyphaeN[i]/everything$hyphaeN[j]
        everything$mycosNratio[i] = everything$mycosN[i]/everything$mycosN[j]
        everything$nmNratio[i] = everything$uncolN[i]/everything$uncolN[j]
        everything$pctNratiomycos[i] = everything$mycorrhizas[i]/everything$mycorrhizas[j]
        everything$pctNratiouncol[i] = everything$uncolonized_roots[i]/everything$uncolonized_roots[j]
        
        if (everything$Actual_fungus_by_compartment[i] == "THETE" | everything$Actual_fungus_by_compartment[i] == "SUIPU" &
            everything$Actual_fungus_by_compartment[j] == "THETE" | everything$Actual_fungus_by_compartment[j] == "SUIPU") {
          everything$absorpCratio[i] =  everything$mycosC[i]/everything$mycosC[j]
        } else if (everything$Actual_fungus_by_compartment[i] == "THETE" | everything$Actual_fungus_by_compartment[i] == "SUIPU" &
                   everything$Actual_fungus_by_compartment[j] == "NM") {
          everything$absorpCratio[i] =  everything$mycosC[i]/everything$uncolC[j]
          
        } else if (everything$Actual_fungus_by_compartment[i] == "NM" &
                   everything$Actual_fungus_by_compartment[j] == "THETE" | everything$Actual_fungus_by_compartment[j] == "SUIPU") {
          everything$absorpCratio[i] =  everything$uncolC[i]/everything$mycosC[j]
        } else if (everything$Actual_fungus_by_compartment[i] == "NM" &
                  everything$Actual_fungus_by_compartment[j] == "NM") {
          everything$absorpCratio[i] =  everything$uncolC[i]/everything$uncolC[j]
          
        }
        
      }
    }
  }
  
}

everything = subset(everything, Actual_fungi_at_harvest != "FAILED" & Actual_fungi_at_harvest != "THETE")

# I would like to calculate the log ratio
# of 13C in mycorrhizas vs 15N in coarse roots. 

ggplot(data = subset(everything, `Receives 15N label?` == "Y") )+
  geom_point(aes(y = log(nmCratio), x = log(pctNratiouncol), color = Actual_fungi_at_harvest))

ggplot(data = subset(everything, `Receives 15N label?` == "Y" & Plant != 6041) ) +
  ylim(-25, 25) +
  geom_point(aes(x = (mycosCratio), y = (nmNratio), color = Actual_fungus_by_compartment))




fcplot = ggplot(data = isowithmeta) +
  theme_classic() +
  theme(axis.text.x = element_text(face = "italic")) +
  geom_boxplot(aes(x = as.factor(enriched), y = d13C))
  # geom_point(aes(color = host, y = d15N_coarse, x = mycorrhized)) +

pdf("plots/did_13C_enrichmentwork.pdf", width = 7, height = 5)
fcplot
dev.off()

fcplot = ggplot(data = isowithmeta) +
  theme_classic() +
  theme(axis.text.x = element_text(face = "italic")) +
  geom_boxplot(aes(x = as.factor(enriched), y = log(d15N)))
# geom_point(aes(color = host, y = d15N_coarse, x = mycorrhized)) +

pdf("plots/did_15N_enrichmentwork.pdf", width = 7, height = 5)
fcplot
dev.off()

fcplot = ggplot(data = isowithmeta) +
  theme_classic() +
  theme(axis.text.x = element_text(face = "italic")) +
  geom_boxplot(aes(x = as.factor(enriched), y = log(d15N)))
# geom_point(aes(color = host, y = d15N_coarse, x = mycorrhized)) +

pdf("plots/did_15N_enrichmentwork.pdf", width = 7, height = 5)
fcplot
dev.off()

fcplot = ggplot(data = subset(isowithmeta, enriched == "1")) +
  theme_classic() +
  theme(axis.text.x = element_text(face = "italic")) +
  geom_point(aes(x = log(d13C + 50), y = log(d15N+10), color = tissue, shape = N_level))
# geom_point(aes(color = host, y = d15N_coarse, x = mycorrhized)) +

pdf("plots/exploring_enrichment_by_N_and_tissue_loglog.pdf", width = 7, height = 5)
fcplot
dev.off()

fcplot = ggplot(data = subset(isowithmeta, enriched == "1" & 
                                (Actual_fungus_by_compartment == "THETE" | Actual_fungus_by_compartment == "SUIPU"))) +
  theme_classic() +
  theme(axis.text.x = element_text(face = "italic")) +
  geom_point(aes(x = log(d13C + 30), y = log(d15N), shape = tissue, color = Actual_fungus_by_compartment))
# geom_point(aes(color = host, y = d15N_coarse, x = mycorrhized)) +

pdf("plots/exploring_enrichment_byfungus.pdf", width = 7, height = 5)
fcplot
dev.off()

fcplot = ggplot(data = isowithmeta) +
  theme_classic() +
  theme(axis.text.x = element_text(face = "italic")) +
  geom_point(aes(x = log(d13C + 30), y = log(d15N), shape = tissue, color = enriched))
# geom_point(aes(color = host, y = d15N_coarse, x = mycorrhized)) +

pdf("plots/didenrichmentwork.pdf", width = 7, height = 5)
fcplot
dev.off()

mysumm = isowithmeta %>% group_by(Plant, Side)

isowithmeta$`Receives 15N label?`

smallerset = select(isowithmeta, Plant,
                    Side, tissue, Amount,
                    d13C, pctC, d15N, pctN,
                    CNratio, N_level, `Receives 15N label?`,
                    Batch, Actual_fungus_by_compartment, enriched)

small = select(isowithmeta, Plant, Side, tissue, d13C)
# small$uniqueId = c(1:nrow(small))
smallspread = small %>% spread(tissue, d13C)

# So, this finally worked (now that I have only
# one entry per seedling/side/tissue)...
# But I'm realizing these delta values are really nothing like
# what I need. I need atom % or mmolar  enrichment
# to actually be able to compare with all these
# ratios that I'm interested in.
