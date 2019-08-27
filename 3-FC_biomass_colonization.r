# Analyzing FC biomass and colonization data

# setwd("~/Documents/2018-2019/Fungal competition/")

require(tidyverse)
require(cowplot)
require(agricolae)

### Need to change fungi codes to match Tt/Sp from manuscript

fungcomp = read_csv("./FCdata/Fungal_competition_FC_seedling_masses_final.csv")
metadata = read_csv("./FCdata/Fungal_competition_plant_tracking.csv")

fungcomp = select(fungcomp, Page:Notes2)

# Cleaning up coding for data processing
fungcomp$for_percent_col = numeric(nrow(fungcomp))
fungcomp$for_percent_col[grep("%", fungcomp$Tissue)] = 1

fungcomp$coded_tissue = numeric(nrow(fungcomp))
fungcomp$coded_tissue = "other" # default assume this is a crust or rhizomorph or something weird
fungcomp$coded_tissue[grep("TFR", fungcomp$Tissue)] = "roots" # "terminal fine roots."
fungcomp$coded_tissue[grep("root", fungcomp$Tissue)] = "roots"
fungcomp$coded_tissue[grep("no mycos", fungcomp$Tissue)] = "roots"
fungcomp$coded_tissue[grep("other %", fungcomp$Tissue)] = "roots"
fungcomp$coded_tissue[grep("NM", fungcomp$Tissue)] = "roots"
fungcomp$coded_tissue[grep("extra", fungcomp$Tissue)] = "roots"
fungcomp$coded_tissue[grep("mycos", fungcomp$Tissue)] = "mycorrhizas"
fungcomp$coded_tissue[grep("suipu", fungcomp$Tissue)] = "mycorrhizas"
fungcomp$coded_tissue[grep("SUIPU", fungcomp$Tissue)] = "mycorrhizas"
fungcomp$coded_tissue[grep("dead", fungcomp$Tissue)] = "dead_tissue"
fungcomp$coded_tissue[grep("scenescent", fungcomp$Tissue)] = "dead_tissue"
fungcomp$coded_tissue[grep("before", fungcomp$Tissue)] = "shoot" # needles taken pre-13C label
fungcomp$coded_tissue[grep("shoot", fungcomp$Tissue)] = "shoot"

# Some weird mass entries, too:
fungcomp$Mass_g[grep("not weighed", fungcomp$Mass_g)] = NA

# Make classes merge-able.
fungcomp$Mass_g = as.numeric(fungcomp$Mass_g)

# Imported some NA lines for some reason. Remove.
fungcomp = fungcomp[!is.na(fungcomp$Plant_number),]

# Check for duplicate entries
shouldbeunique = paste(fungcomp$Plant_number, fungcomp$Side, fungcomp$Tissue)

sum(duplicated(shouldbeunique)) # 1 putative duplicate

fungcomp[duplicated(shouldbeunique),]
target = shouldbeunique[duplicated(shouldbeunique)] # 6105 B mycos %
suspect = fungcomp[shouldbeunique == target,] # Arvie and Caroline
# both weighed this. Masses are different, but only by about 6 mg.
# Oh! Actually, these were part of a test to see if 
# isotopes changed after storage in water for 
# different amounts of time. They really were different samples.
# I guess they can stay!

# Metadata cleaning:
metadata = metadata[!is.na(metadata$Plant),] # read in a bunch of NA rows again for some reason.
metadata$enriched = numeric(nrow(metadata)) # code in enrichment vs not.
metadata$enriched[grep("FC", metadata$Batch)] = 1
metadata$Actual_fungi_at_harvest[metadata$Plant == 6033] = "SUIPU/SUIPU" # This was entered wrong for some reason.
metadata$Actual_fungi_at_harvest[metadata$Plant == 6057] = "MIXED/THETE" # This was entered wrong for some reason.
metadata$Failed_split[is.na(metadata$Failed_split)] = "N"
metadata$Failed_split[metadata$Plant == 6009 | metadata$Plant == 6081] = "Y"

wrongfung = as.character(metadata$Fungal_treatment) != as.character(metadata$Actual_fungi_at_harvest)
# This is worth knowing, but I think for now it makes sense
# to treat each plant as belonging to the treatment
# reflecting the fungi it actually had at harvest.

# I have a more detailed treatment of this in
# Fungal_competition_analysis_for_committee.r

#### BIOMASS ####

plant_mass = fungcomp %>% group_by(Plant_number) %>% summarize(total_biomass = sum(Mass_g, na.rm = TRUE))

# I'd like to look at mycorrhizal growth response,
# but I think this really only makes sense
# to do WITHIN nitrogen treatments.

plantleveldata = select(metadata, Plant:N_level, `Harvest date`:Actual_fungi_at_harvest, Harvest_notes, enriched)
plantleveldata = plantleveldata[!duplicated(plantleveldata$Plant),]

plant_mass = rename(plant_mass, Plant = Plant_number)

plantleveldata = left_join(plantleveldata, plant_mass)

NMbaseline = subset(plantleveldata, Actual_fungi_at_harvest == "NM/NM")

highnm = subset(NMbaseline, N_level == "High") # literally one plant
lownm = subset(NMbaseline, N_level == "Low") # two plants

highnmmass = highnm$total_biomass
lownmmass = mean(lownm$total_biomass)

plantleveldata$plant_response = numeric(nrow(plantleveldata))
for (i in 1:nrow(plantleveldata)) {
  if (plantleveldata$N_level[i] == "High") {
    plantleveldata$plant_response[i] = log(plantleveldata$total_biomass[i]/highnmmass)
  } else if (plantleveldata$N_level[i] == "Low") {
    plantleveldata$plant_response[i] = log(plantleveldata$total_biomass[i]/lownmmass)
  }
}

test = plantleveldata[is.na(plantleveldata$plant_response),]
# All of these are missing, for whatever reason. Probably
# the plants died prematurely.

plantleveldata = plantleveldata[!is.na(plantleveldata$plant_response),]
plantleveldata = rename(plantleveldata, Fungi = Actual_fungi_at_harvest)

plantleveldata = subset(plantleveldata, N_level != "None")

# write_csv(plantleveldata, "./FCdata/percent_col_and_mass_by_plant_including_mixed.csv")

plantleveldata = plantleveldata[-grep("MIXED", plantleveldata$Fungi),]
plantleveldata = plantleveldata[-grep("FAILED", plantleveldata$Fungi),]
plantleveldata = subset(plantleveldata, Fungi != "SUIPU" & Fungi != "THETE") # failed splits

plantleveldata$Fungi = factor(plantleveldata$Fungi,
                                 levels = c("NM/NM",
                                            "SUIPU/NM",
                                            "SUIPU/SUIPU",
                                            "THETE/SUIPU",
                                            "THETE/NM",
                                            "THETE/THETE"))

plantleveldata$Fungi = recode(plantleveldata$Fungi,
               "NM/NM" = "None/None",
               "SUIPU/NM" = "Sp/None",
               "SUIPU/SUIPU" = "Sp/Sp",
               "THETE/SUIPU" = "Tt/Sp",
               "THETE/NM" = "Tt/None",
               "THETE/THETE" = "Tt/Tt")


#### REPLICATION SUMMARY ####
summarytable = plantleveldata %>% group_by(Fungi, N_level, enriched) %>% summarize(count = n())
# write_csv(summarytable, "harvest_replication_summary.csv")
nrow(plantleveldata)
nrow(plantleveldata[plantleveldata$N_level == "High",])
nrow(plantleveldata[plantleveldata$N_level == "Low",])
nrow(plantleveldata[plantleveldata$enriched == 1,])

#### PLANT RESPONSE TO COLONIZATION ####

forplot = subset(plantleveldata, Fungi != "None/None")


tx = with(forplot, interaction(N_level, Fungi))
anovaforplot = aov(plant_response ~ tx, data = forplot)

# from "agricolae" package
mylabels = HSD.test(anovaforplot, "tx", group = TRUE)

annotations = data.frame(x = c((1:5), (1:5)),
                        y = c(1, 1, 1.3, 1.3, 1.7, 1, 1, 1, 1, 1),
                        N_level = c(rep("High", 5), rep("Low", 5)),
                        labs = c(paste(c("bc", "abc", "ab", "bc", "a")), paste(c("c", "bc", "bc", "bc", "bc"))))


labels = c(High = "High N", Low = "Low N")

responseplot = ggplot(data = forplot) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = Fungi, y = plant_response)) +
  geom_jitter(width = 0.20,
              aes(x = Fungi, y = plant_response)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Plant response to fungal\ncolonization (log response ratio)") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  geom_text(data = annotations, aes(x, y, label = labs)) +
  xlab("Fungal treatment")


pdf("plots/Plant_response_boxplot.pdf", width = 7, height = 5)
responseplot
dev.off()

responseanova = aov(plant_response ~ N_level * Fungi, data = forplot)
summary(responseanova)
responseTukey = TukeyHSD(responseanova)

write.csv(responseTukey$N_level, "Statistical_tables/Plant_response_Tukey_Nlevel.csv")
write.csv(responseTukey$Fungi, "Statistical_tables/Plant_response_Tukey_Fungi.csv")
write.csv(responseTukey$`N_level:Fungi`, "Statistical_tables/Plant_response_Tukey_Nlevel-Fungi.csv")

# Do these differ from zero?
summary(forplot$Fungi)

# High N:
hightx = subset(forplot, N_level == "High")
t.test(hightx$plant_response[hightx$Fungi == "Sp/None"]) # not enough observations
t.test(hightx$plant_response[hightx$Fungi == "Sp/Sp"]) # p = 0.1185
t.test(hightx$plant_response[hightx$Fungi == "Tt/Sp"]) # p = 0.02562
t.test(hightx$plant_response[hightx$Fungi == "Tt/None"]) # p = 0.1547
t.test(hightx$plant_response[hightx$Fungi == "Tt/Tt"]) # p = 7.478e-11

# Low N:
lowtx = subset(forplot, N_level == "Low")
t.test(lowtx$plant_response[lowtx$Fungi == "Sp/None"]) # p = 0.1654
t.test(lowtx$plant_response[lowtx$Fungi == "Sp/Sp"]) # p = 0.6937
t.test(lowtx$plant_response[lowtx$Fungi == "Tt/Sp"]) # p = 0.2881
t.test(lowtx$plant_response[lowtx$Fungi == "Tt/None"]) # p = 0.399
t.test(lowtx$plant_response[lowtx$Fungi == "Tt/Tt"]) # p = 0.1192

# Alpha should be 0.05/10 tests = 0.005, per Bonferroni correction. So the only significant 
# non-zero growth effect difference is Tt/Tt in high N





#### BIOMASS ####
# This function maybe does automatic letters?
tx = with(plantleveldata, interaction(N_level, Fungi))
anovaforplot = aov(total_biomass ~ tx, data = plantleveldata)

# from "agricolae" package
mylabels = HSD.test(anovaforplot, "tx", group = TRUE)
# Oh thank goodness this matches my prior
# results and simplifies my life a lot.

anothertry = data.frame(x = c((1:6), (1:6)),
                        y = c(4, 4, 4, 6, 6, 7.5, 4, 4, 4, 4, 4, 4),
                        N_level = c(rep("High", 6), rep("Low", 6)),
                        labs = c(paste(c("bc", "c", "bc", "b", "b", "a")), paste(c("c", "c", "c", "c", "c", "c"))))

labels = c(High = "High N", Low = "Low N")
massplot = ggplot(data = plantleveldata) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = Fungi, y = total_biomass)) +
  geom_jitter(width = 0.20,
             aes(x = Fungi, y = total_biomass)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Total plant biomass (g)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  geom_text(data = anothertry, aes(x, y, label = labs)) +
  xlab("Fungal treatment")

pdf("plots/Biomass_boxplot.pdf", width = 9, height = 5)
massplot
dev.off()

biomassanova = aov(total_biomass ~ N_level * Fungi, data = plantleveldata)
summary(biomassanova)
biomassTukey = TukeyHSD(biomassanova)

write.csv(biomassTukey$N_level, "Statistical_tables/Biomass_Tukey_Nlevel.csv")
write.csv(biomassTukey$Fungi, "Statistical_tables/Biomass_Tukey_Fungi.csv")
write.csv(biomassTukey$`N_level:Fungi`, "Statistical_tables/Biomass_Tukey_Nlevel-Fungi.csv")



#### COLONIZATION ####


protocol = select(fungcomp, Plant_number, coded_tissue, Mass_g, for_percent_col)
protocol = subset(protocol, for_percent_col == 1)

oneentrypertissue = protocol %>% group_by(Plant_number, coded_tissue) %>% summarize(mass = sum(Mass_g, na.rm = TRUE))

widedata = oneentrypertissue %>% spread(coded_tissue, mass)

widedata$mycorrhizas[is.na(widedata$mycorrhizas)] = 0 # if na, there were no mycorrhizas

# widedata$percentcoldenom = sum(widedata$mycorrhizas, widedata$roots, widedata$other, na.rm = TRUE)

widedata$percent_col = numeric(nrow(widedata))
for (i in 1:nrow(widedata)) {
    percentcoldenom = sum(widedata$mycorrhizas[i], widedata$roots[i], widedata$other[i], na.rm = TRUE)
    widedata$percent_col[i] = 100*(widedata$mycorrhizas[i]/percentcoldenom)
}

widedata = rename(widedata, Plant = Plant_number)

alldata = left_join(plantleveldata, widedata)

colforplot = subset(alldata, Fungi != "None/None")
labels = c(High = "High N", Low = "Low N")

tx = with(colforplot, interaction(N_level, Fungi))
anovaforplot = aov(percent_col ~ tx, data = colforplot)

# from "agricolae" package
mylabels = HSD.test(anovaforplot, "tx", group = TRUE)

# anothertry = data.frame(x = c((1:6), (1:6)),
#                         y = c(4, 4, 4, 6, 6, 7.5, 5, 4, 5, 4, 4, 4),
#                         N_level = c(rep("High", 6), rep("Low", 6)),
#                         labs = c(paste(c("c", "c", "c", "b", "b", "a")), paste(c("bc", "c", "bc", "c", "c", "c"))))

# Okay, actually there are no significant pairwise differences here.
# Maybe no need for the Tukey labels.


collabels = data.frame(N_level = c("High", "Low"),
                       x1 = c(1, 1), x2 = c(5, 5), y1 = c(80, 50), y2 = c(81, 51), 
                       xstar = c(3, 3), ystar = c(88, 58), 
                       lab = c("a", "b"))


colplot = ggplot(data = colforplot) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = Fungi, y = percent_col)) +
  geom_jitter(width = 0.20,
              aes(x = Fungi, y = percent_col)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Percent fungal colonization of\nroot system (by mass)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  xlab("Fungal treatment") +
  geom_text(data = collabels, aes(x = xstar,  y = ystar, label = lab)) +
  geom_segment(data = collabels, aes(x = x1, xend = x2, 
                                y = y2, yend = y2),
               colour = "black")


pdf("plots/Colonization_boxplot.pdf", width = 7, height = 5)
colplot
dev.off()

plot_grid(massplot, responseplot, colplot, ncol = 3, align = "h")

Figure2 = plot_grid(massplot, colplot, ncol = 2, align = "h",
                    labels = c("A", "B"))
save_plot("plots/MAIN_Mass_and_colonization_two_panel_boxplot.pdf", 
          Figure2, ncol = 2)

colanova = aov(percent_col ~ N_level * Fungi, data = colforplot)
coloutput = summary(colanova) # N level v. significant, fungi only marginal
coltukey = TukeyHSD(colanova)

high = subset(colforplot, N_level == "High")
low = subset(colforplot, N_level == "Low")
median(high$percent_col)
median(low$percent_col)


write.csv(coltukey$N_level, "Statistical_tables/Colonization_Tukey_output_Nlevel.csv")
write.csv(coltukey$Fungi, "Statistical_tables/Colonization_Tukey_output_Fungi.csv")
write.csv(coltukey$`N_level:Fungi`, "Statistical_tables/Colonization_Tukey_output_Nlevel-Fungi.csv")

# So, colonization rates aren't that interesting
# unto themselves... but I think they could
# have important bearing on N uptake
# on a PER COMPARTMENT basis. 


# How can I get those data?

prepfortable = select(fungcomp, Plant_number, Side, coded_tissue, Mass_g, for_percent_col)
prepfortable = subset(prepfortable, for_percent_col == 1)

onepertiss = prepfortable %>% group_by(Plant_number, Side, coded_tissue) %>% summarize(mass = sum(Mass_g, na.rm = TRUE))

wide_bycompt = onepertiss %>% spread(coded_tissue, mass)

wide_bycompt$mycorrhizas[is.na(wide_bycompt$mycorrhizas)] = 0 # if na, there were no mycorrhizas

wide_bycompt$percent_col = numeric(nrow(wide_bycompt))
for (i in 1:nrow(wide_bycompt)) {
  percentcoldenom = sum(wide_bycompt$mycorrhizas[i], wide_bycompt$roots[i], wide_bycompt$other[i], na.rm = TRUE)
  wide_bycompt$percent_col[i] = 100*(wide_bycompt$mycorrhizas[i]/percentcoldenom)
}

wide_bycompt = rename(wide_bycompt, Plant = Plant_number)
tomerge = select(metadata, Plant, Side, N_level, Batch, compartment_fungus = Actual_fungus_by_compartment, competitors = Actual_fungi_at_harvest, enriched)

wide_bycompt = left_join(wide_bycompt, tomerge)

write_csv(wide_bycompt, "./FCdata/percent_colonization_and_mass_data_by_compartment.csv")

#### Did plant repress colonization by less helpful fungus? ####

wide_bycompt$compartment_fungus

bycomptoplot = subset(wide_bycompt, 
                      compartment_fungus != "MIXED" & 
                        is.na(compartment_fungus) == FALSE &
                        compartment_fungus != "OTHER" & 
                        competitors != "FAILED" &
                        N_level != "None")

bycomptoplot = bycomptoplot[-grep("MIXED", bycomptoplot$competitors),]
bycomptoplot = subset(bycomptoplot, compartment_fungus != "NM")

ggplot(data = subset(bycomptoplot, competitors != "THETE/SUIPU")) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = N_level, y = percent_col, color = compartment_fungus)) +
  geom_jitter(width = 0.20,
              aes(x = N_level, y = percent_col, color = compartment_fungus))




TtSp = subset(bycomptoplot, competitors == "THETE/SUIPU")

ggplot(data = TtSp) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = N_level, y = percent_col, color = compartment_fungus)) +
  geom_jitter(width = 0.20,
              aes(x = N_level, y = percent_col, color = compartment_fungus))

require(ggpubr)
ggpaired(data = subset(TtSp, N_level == "High"), x = "compartment_fungus", y = "percent_col", 
         id = "Plant")
ggpaired(data = subset(TtSp, N_level == "Low"), x = "compartment_fungus", y = "percent_col", 
         id = "Plant")

ggplot(data = TtSp) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = N_level, y = percent_col, color = compartment_fungus)) +
  geom_point(width = 0.20,
              aes(x = N_level, y = percent_col, color = compartment_fungus)) +
  geom_line(aes(group = Plant))


test = aov(percent_col ~ compartment_fungus*N_level, data = TtSp)
summary(test)
TukeyHSD(test)

test = aov(percent_col ~ N_level*compartment_fungus, data = TtSp)
summary(test)

TtSp$diff = numeric(nrow(TtSp))
TtSp$logprop = numeric(nrow(TtSp))

for (i in 1:nrow(TtSp)) {
  for (j in 1:nrow(TtSp)) {
    if (TtSp$Plant[i] == TtSp$Plant[j]) {
      if (TtSp$compartment_fungus[i] == "THETE" &
          TtSp$compartment_fungus[j] == "SUIPU") {
        TtSp$diff[i] = TtSp$percent_col[i]-TtSp$percent_col[j]
        TtSp$logprop[i] = log(TtSp$percent_col[i]/(TtSp$percent_col[j]+1))
      }
    }
  }
}

TtSp_compare = TtSp[TtSp$diff != 0,]

ggplot(data = TtSp_compare) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = N_level, y = diff)) +
  geom_jitter(width = 0.20,
              aes(x = N_level, y = diff))

t.test(diff ~ N_level, TtSp_compare)

#### Does colonization predict biomass? ####

alldata$Fungus = recode(alldata$Fungi,
                        "None/None" = "None",
                        "Sp/None" = "Sp",
                        "Sp/Sp" = "Sp",
                        "Tt/Sp" = "Tt/Sp",
                        "Tt/None" = "Tt",
                        "Tt/Tt" = "Tt")

write_csv(alldata, "percent_col_and_mass_data_by_plant.csv")


massbycol_plot = ggplot(data = alldata) +
  geom_point(aes(x = percent_col,
                 y = total_biomass,
                 shape = Fungus,
                 color = N_level)) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
    geom_smooth(method = "lm",
                formula = y ~ log(x + 1),
              aes(x = percent_col,
                  y = total_biomass),
              color = "black",
              size = 0.5) +
  scale_shape_manual(values = c(23, 17, 16, 15)) +
  ylab("Total plant biomass (g)") +
  xlab("Percent colonization of root system") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

# Wow, nope, not well.
# Looks sigmoid/saturating to me, though?

massbycol = lm(total_biomass ~ log(percent_col + 1), data = alldata)
summary(massbycol)
plot(massbycol)
# It is better when I log transform percent col, for sure

## EXCLUDING OUTLIER FOR TAD ##
alldata[alldata$percent_col == max(alldata$percent_col),] # plant 6040 is one Tad suggests might be an outlier
noout_alldata = subset(alldata, percent_col < 60)

massbycol_noout = lm(total_biomass ~ percent_col, data = noout_alldata)
summary(massbycol_noout)
plot(massbycol_noout)
# When removing the outlier, it's no longer convincingly better 
# to log-transform.


massbycol_noout = lm(total_biomass ~ log(percent_col + 1), data = noout_alldata)
summary(massbycol_noout)
plot(massbycol_noout)

massbycol_plot_noout = ggplot(data = noout_alldata) +
  geom_point(aes(x = percent_col,
                 y = total_biomass,
                 shape = Fungus,
                 color = N_level)) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  geom_smooth(method = "lm",
              formula = y ~ x,
              aes(x = percent_col,
                  y = total_biomass),
              color = "black",
              size = 0.5) +
  scale_shape_manual(values = c(23, 17, 16, 15)) +
  ylab("Total plant biomass (g)") +
  xlab("Percent colonization of root system") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))



nonm = subset(alldata, Fungi != "None/None")
respbycol = lm(plant_response ~ log(percent_col + 1), data = nonm)
plot(respbycol) # I am not sure the log transformation improves this, but it doesn't make it worse.
summary(respbycol)

respbycol_plot = ggplot(data = nonm) +
  geom_point(aes(x = percent_col,
                 y = plant_response,
                 shape = Fungus,
                 color = N_level)) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  geom_smooth(method = "lm",
              formula = y ~ log(x + 1),
              aes(x = percent_col,
                  y = plant_response),
              color = "black",
              size = 0.5) +
  scale_shape_manual(values = c(17, 16, 15)) +
  ylab("Plant response to fungal\ncolonization (log response ratio)") +
  xlab("Percent colonization of root system") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
  
  # geom_hline(yintercept = 0, linetype = "dashed", size = 0.5)

## REMOVING EXTREME VALUE FOR TAD ##
noout_nonm = subset(noout_alldata, Fungi != "None/None")
respbycol_nonm = lm(plant_response ~ log(percent_col + 1), data = noout_nonm)
plot(respbycol_nonm) # I am not sure the log transformation improves this

respbycol_nonm = lm(plant_response ~ percent_col, data = noout_nonm)
plot(respbycol_nonm) 

summary(respbycol_nonm)

noout_respbycol_plot = ggplot(data = noout_nonm) +
  geom_point(aes(x = percent_col,
                 y = plant_response,
                 shape = Fungus,
                 color = N_level)) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  geom_smooth(method = "lm",
              formula = y ~ x,
              aes(x = percent_col,
                  y = plant_response),
              color = "black",
              size = 0.5) +
  scale_shape_manual(values = c(17, 16, 15)) +
  ylab("Plant response to fungal\ncolonization (log response ratio)") +
  xlab("Percent colonization of root system") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

massandrespbycol = plot_grid(massbycol_plot, respbycol_plot, labels = c("A", "B"), align = "h", ncol = 2)
save_plot("plots/Mass_and_plant_response_by_percent_colonization.pdf", massandrespbycol, 
          ncol = 2,
          base_aspect_ratio = 1.3)

massandrespbycol_noout = plot_grid(massbycol_plot_noout, noout_respbycol_plot, labels = c("A", "B"), align = "h", ncol = 2)
save_plot("plots/FOR_SUPP_mass_and_plant_response_by_percent_colonization.pdf", massandrespbycol_noout, 
          ncol = 2,
          base_aspect_ratio = 1.3)

nonmnomix = subset(nonm, Fungus != "Tt/Sp")
sponly = subset(alldata, Fungus == "Sp")
ggplot(data = sponly) +
  geom_point(aes(x = percent_col,
                 y = plant_response,
                 shape = Fungus,
                 color = N_level)) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  geom_smooth(method = "lm", 
              aes(x = percent_col,
                  y = plant_response),
              color = "black",
              size = 0.5) +
  scale_shape_manual(values = c(17, 16, 15)) +
  ylab("Plant response to\ncolonization (log ratio)") +
  xlab("Percent colonization of root system") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

ttonly = subset(alldata, Fungus == "Tt")
ggplot(data = ttonly) +
  geom_point(aes(x = percent_col,
                 y = plant_response,
                 shape = Fungus,
                 color = N_level)) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  geom_smooth(method = "lm", 
              aes(x = percent_col,
                  y = plant_response),
              color = "black",
              size = 0.5) +
  scale_shape_manual(values = c(17, 16, 15)) +
  ylab("Plant response to\ncolonization (log ratio)") +
  xlab("Percent colonization of root system") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

sponly_model = lm(plant_response ~ percent_col, data = sponly)
summary(sponly_model) # No relationship when looking just at sp.

highNonlytt = subset(ttonly, N_level == "High")
ttonly_model = lm(total_biomass ~ percent_col, data = highNonlytt)
summary(ttonly_model) #

massbycol = lm(total_biomass ~ log(percent_col + 1), data = alldata)
summary(massbycol)
plot(massbycol)

acctforeverything_plantresp = lm(plant_response ~ N_level*log(percent_col + 1)*Fungus, data = nonm)
summary(acctforeverything_plantresp)
summary(aov(acctforeverything_plantresp))

acctforeverything_biomass = lm(total_biomass ~ N_level*log(percent_col + 1)*Fungus, data = alldata)
summary(acctforeverything_biomass)
summary(aov(acctforeverything_biomass))

acctforeverything_biomass_switchorder = lm(total_biomass ~ log(percent_col + 1)*N_level*Fungus, data = alldata)
summary(aov(acctforeverything_biomass_switchorder))

acctforeverything_biomass_switchorder2 = lm(total_biomass ~ Fungus*log(percent_col + 1)*N_level, data = alldata)
summary(aov(acctforeverything_biomass_switchorder2))


plantresponset1ss = aov(plant_response ~ log(percent_col + 1)*N_level, data = nonm)
summary(plantresponset1ss) # Gaaah okay an ANOVA here shows
# both factors as highly significant



alldata$propcol = alldata$percent_col/100

atest = glm(propcol ~ total_biomass, family = binomial, data = alldata)
# This doesn't feel like the right direction
# for a relationship, and is not significant.

alinearoption = glm(total_biomass ~ propcol, data = alldata)
# significant when done in linear manner.

# What about plant response?

ggplot(data = alldata) +
  geom_point(aes(x = percent_col,
                 y = plant_response,
                 shape = N_level,
                 color = Fungi))

#### TRANSITIONS TO WRONG FUNGI ####

transitions_tbl = select(plantleveldata, plant = Plant,
                         N_level,
                         attempted = Fungal_treatment, 
                         realized = Fungi)

transsummary = transitions_tbl %>% group_by(attempted, realized) %>% summarize(count = n())

write_csv(transsummary, "Fungal_transitions_summary_table.csv")

# How did I lose 18 plants?

summary(metadata$Failed_split == "Y") # 8/2 = 4 failed splits
hmm = metadata[metadata$Failed_split == "Y",]


sum(grepl("MIXED|OTHER", metadata$Actual_fungi_at_harvest)) #16/2 = 8 microcosms where we had a mixed result
more = metadata[grepl("MIXED|OTHER|FAILED", metadata$Actual_fungi_at_harvest),]

# I'm still missing six plants that need explaining

letslook = metadata[!metadata$Plant %in% plantleveldata$Plant,]
letslook = letslook[!letslook$Plant %in% more$Plant,]
letslook = letslook[!letslook$Plant %in% hmm$Plant,]

# > length(unique(more$Plant))
# [1] 9
# > length(unique(hmm$Plant))
# [1] 6
# > length(unique(letslook$Plant))
# [1] 4 # but one of these doesn't count because it was the "no N" control plant

#### BRINGING TOGETHER MASS, COLONIZATION, AND PLANT RESPONSE

massandcol = plot_grid(massplot, colplot, labels = c("A", "B"), align = "h")
save_plot("plots/Mass_and_colonization_two_panel_plot.pdf", massandcol, ncol = 2)

# The below does not look as good as the one that I saved with the pdf device directly.
# save_plot("plots/Plant_response_to_colonization_boxplot.pdf", responseplot)

#### Plant response by fungal species and nitrogen level ####

annotations = data.frame(x = c((1:2), (1:2)),
                         y = c(1, 1.5, 1, 1),
                         N_level = c(rep("High", 2), rep("Low", 2)),
                         labs = c("ab", "a", "b", "b"))

justself = subset(nonm, Fungi == "Sp/Sp" | Fungi == "Tt/Tt")
nomix = subset(nonm, Fungi != "Tt/Sp")

ggplot(data = nomix) +
  geom_boxplot(aes(x = Fungus,
                   y = plant_response)) +
  geom_point(aes(x = Fungus,
                 y = plant_response)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  xlab("Fungi on root system") +
  ylab("Plant response to colonization\n(log response ratio)") +
  geom_text(data = annotations, aes(x, y, label = labs))
  
responsebyfungi_selfonly = aov(plant_response ~ Fungus*N_level, data = nomix)
plot(responsebyfungi_selfonly) # okay, I guess?
summary(responsebyfungi_selfonly) # both factors highly significant.

responsebyfungi_selfonly_Tukey = TukeyHSD(responsebyfungi_selfonly)
write.csv(responsebyfungi_selfonly_Tukey$N_level, "Statistical_tables/responsebyfungi_selfonlyNlevel_Tukey_output_Nlevel.csv")
write.csv(responsebyfungi_selfonly_Tukey$Fungus, "Statistical_tables/responsebyfungi_selfonlyNlevel_Tukey_output_Fungi.csv")
write.csv(responsebyfungi_selfonly_Tukey$`N_level:Fungus`, "Statistical_tables/responsebyfungi_selfonlyNlevel_Tukey_output_Nlevel-Fungi.csv")

tx = with(nomix, interaction(N_level, Fungus))

forlabels = aov(plant_response ~ tx, data = nomix)
mylabels = HSD.test(forlabels, "tx", group = TRUE)

annotations = data.frame(x = c((1:3), (1:3)),
                         y = c(1, 1.5, 1.5, 1, 1, 1),
                         N_level = c(rep("High", 3), rep("Low", 3)),
                         labs = c("abc", "ab", "a", "c", "bc", "bc"))

responsebyfungi_plot = ggplot(data = nonm) +
  geom_boxplot(aes(x = Fungus,
                   y = plant_response)) +
  geom_point(aes(x = Fungus,
                 y = plant_response)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  xlab("Fungi on root system") +
  ylab("Plant response to colonization\n(log response ratio)") +
  geom_text(data = annotations, aes(x, y, label = labs))

pdf("plots/Plant_response_by_fungi_and_N_level_simplified.pdf", width = 7, height = 5)
responsebyfungi_plot
dev.off()

# Overall, though, it must be said that this plot tells EXACTLY 
# the same story as my current plant response plot:
# Thelephora helps plants significantly more in high N
# than low; Suillus does not.
# I think I should actually use only one of these figures. 

responsebyfungi = aov(plant_response ~ Fungus*N_level, data = nomix)
plot(responsebyfungi) # okay, I guess?
summary(responsebyfungi) # both factors highly significant.

responsebyfungi_Tukey = TukeyHSD(responsebyfungi)
write.csv(responsebyfungi_Tukey$N_level, "Statistical_tables/ResponsebyfungiNlevel_Tukey_output_Nlevel.csv")
write.csv(responsebyfungi_Tukey$Fungus, "Statistical_tables/ResponsebyfungiNlevel_Tukey_output_Fungi.csv")
write.csv(responsebyfungi_Tukey$`N_level:Fungus`, "Statistical_tables/ResponsebyfungiNlevel_Tukey_output_Nlevel-Fungi.csv")
