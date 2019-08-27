# FC isotope analyses cleaner

# setwd("~/Documents/2018-2019/Fungal competition/")

require(tidyverse)
require(cowplot)
require(agricolae) # for automatic Tukey labels

# allisotopes = read_csv("isotope_data_two_rows_per_plant.csv")
allisotopes = read_csv("./FCdata/isotope_data_two_rows_per_plant_July.csv")
# correcting a spreadsheet error
allisotopes$Actual_fungi_at_harvest[allisotopes$Plant == 6033] = "SUIPU/SUIPU"

isotopes_forN = read_csv("./FCdata/isotope_data_one_row_per_plant_July.csv")
minimally_processed_isotopes = read_csv("./FCdata/Cleaned_processed_FC_isotope_data_July.csv")
percent_col = read_csv("./FCdata/percent_colonization_and_mass_data_by_compartment.csv")
metadata_byplant = read_csv("./FCdata/percent_col_and_mass_data_by_plant.csv")
# mydata = read_csv("minimally_processed_isotope_data.csv")

allisotopes = rename(allisotopes, compartment_fungus = Actual_fungus_by_compartment)
allisotopes = select(allisotopes, everything(), -tissue)
isotopes_forN = rename(isotopes_forN, compartment_fungus = mycofungus)
isotopes_forN = select(isotopes_forN, everything(), -tissue)

minimally_processed_isotopes = rename(minimally_processed_isotopes, compartment_fungus = Actual_fungus_by_compartment)
# minimally_processed_isotopes = subset(minimally_processed_isotopes, Failed_split != "Y")

percent_col$Side = tolower(percent_col$Side)

together = left_join(allisotopes, percent_col)

together = rename(together, Fungi = Actual_fungi_at_harvest)

together = subset(together, Failed_split != "Y")

together$Fungi = recode(together$Fungi,
                         "NM/NM" = "None/None",
                         "SUIPU/NM" = "Sp/None",
                         "SUIPU/SUIPU" = "Sp/Sp",
                         "THETE/SUIPU" = "Tt/Sp",
                         "THETE/NM" = "Tt/None",
                         "THETE/THETE" = "Tt/Tt")
# 
# together$Fungus = recode(together$Fungi,
#                         "None/None" = "None",
#                         "Sp/None" = "Sp",
#                         "Sp/Sp" = "Sp",
#                         "Tt/Sp" = "Tt/Sp",
#                         "Tt/None" = "Tt",
#                         "Tt/Tt" = "Tt")

together$compartment_fungus = recode(together$compartment_fungus,
                         "NM" = "None",
                         "SUIPU" = "Sp",
                         "THETE" = "Tt")

together$mycofungus = recode(together$mycofungus,
                                     "NM" = "None",
                                     "SUIPU" = "Sp",
                                     "THETE" = "Tt")

write_csv(together, "isotopes_together_as_analyzed.csv")

#### Checking replication so I can reanalyze samples ####

# Any plant which doesn't have at least an NM root sample
# on both sides is missing something.
# Also, any compartment that is Sp or Tt and doesn't
# have myco data is missing something.



formissing = subset(together, compartment_fungus != "OTHER" & compartment_fungus != "MIXED")
checkingcompleteness = as.data.frame(table(formissing$Plant))
checkingcompleteness = rename(checkingcompleteness, Plant = Var1, numbercompartments = Freq)

checkingcompleteness$incomp = checkingcompleteness$numbercompartments < 2
myincompletes = checkingcompleteness[checkingcompleteness$incomp == TRUE,]

fullincompletes = formissing[formissing$Plant %in% myincompletes$Plant,]
fullincompletes = select(fullincompletes, Plant, Side, N_level, Fungi, compartment_fungus, mycorrhizas.APE13C, uncolonized_roots.APE13C)
fullincompletes = left_join(fullincompletes, select(metadata_byplant, Plant, Batch))
fullincompletes = select(fullincompletes, Batch, everything())


# write_csv(fullincompletes, "plants_needing_other_side.csv")

missingnm = formissing[is.na(formissing$uncolonized_roots.APE13C),]
missingnm$missing = "uncolonized"
missingmyco = subset(formissing, compartment_fungus == "Tt" | compartment_fungus == "Sp")
missingmyco = missingmyco[is.na(missingmyco$mycorrhizas.APE13C),]
missingmyco$missing = "myco"
missingsomething = full_join(missingnm, missingmyco)

tosearch = select(missingsomething, Plant, Side, Fungi, N_level, compartment_fungus, missing)

batchinfo = left_join(tosearch, select(metadata_byplant, Plant, Batch))
batchinfo = select(batchinfo, Batch, everything())

# write_csv(batchinfo, "samples_to_tin_and_submit.csv")

notmissing = formissing[!formissing$Plant %in% batchinfo$Plant,]
notmissing = notmissing[!duplicated(notmissing$Plant),]
# 64 complete plants?

currentreplication = notmissing %>% group_by(Fungi, N_level) %>% summarize(totalreps = n())
# write_csv(currentreplication, "isotope_replication_summary_July.csv")

#### How much did the N-15 label leak across mesh panels? ####

unenriched = subset(minimally_processed_isotopes, enriched == 0 & 
                      (tissue == "uncolonized_roots"|tissue == "mycorrhizas"))

nmunenriched = subset(unenriched, tissue == "uncolonized_roots")

mean(unenriched$APE15N)
sd(unenriched$APE15N)

# How were unenriched plants distributed across treatments?

formerge = select(metadata_byplant, Plant, N_level)
unwithn = left_join(unenriched, formerge)
unwithn_nodup = unwithn[!duplicated(unwithn$Plant),]

unwithn_nodup %>% group_by(Actual_fungi_at_harvest, N_level) %>% summarize(total = n())
# Oh gosh, I've only got three of these. Shouldn't matter a 
# lot, though, since they're just a baseline -- as long
# as all plants have the SAME baseline for this study, it's
# low stakes.

unwithn %>% group_by(tissue, N_level) %>% summarize(total = n())


# How many unenriched plants did I have overall?
summary(as.factor(metadata_byplant$Batch))
# Okay, this is still just three. I have more in my plant spreadsheet.

nm15N = subset(minimally_processed_isotopes, 
               compartment_fungus == "NM" &
                 `Receives 15N label?` == "Y")

mean(nm15N$APE15N)
sd(nm15N$APE15N)


t.test(nm15N$APE15N, unenriched$APE15N)


# This is encouraging:
# > t.test(nm15N$APE15N, unenriched$APE15N)
# 
# Welch Two Sample t-test
# 
# data:  nm15N$APE15N and unenriched$APE15N
# t = -1.0766, df = 8.9805, p-value = 0.3097
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.0016680643  0.0005925751
# sample estimates:
#   mean of x     mean of y 
# -5.377446e-04  1.850373e-17 

#### Creating subsetted data without low APE 15N values ####

# I am a bit concerned about my very low APE 15N values.
# I think they represent situations in which
# the fungus may not have found our label at all.

maxunenriched = max(unenriched$APE15N)

onlyhigh15N = subset(together, 
                     mycorrhizas.APE15N > maxunenriched)

# forplot = as_tibble(rbind(nm15N, unenriched))
# forplot$enriched = as.factor(forplot$enriched)
# leakagecheck = ggplot(data = forplot) +
#   geom_boxplot(aes(x = enriched, y = APE15N)) +
#   geom_point(aes(x = enriched, y = APE15N)) +
#   ylab(expression("Atom percent excess "^15*N))+
#   scale_x_discrete(name = "Treatment",
#                    breaks = c(0, 1),
#                    labels = c("Unenriched roots\nand mycorrhizas", "Uncolonized roots from\nN-15-receiving no-fungus compartments"))
# 
# pdf("plots/No_15N_leakage_boxplot.pdf", width = 5, height = 3)
# leakagecheck
# dev.off()

#### Investigating the outlier ####

# Dealing with compartment 6041b, the consistent outlier:
# My harvest note for this plant
# indicate that it was supposed to be
# SUIPU/NM, and ended up SUIPU/THETE.
# Furthermore, on the THETE side, for which
# we have the labeling data, my notes say
# "not much fungus; looks young/new."
# I could make an argument that this
# sample appears to represent a very different,
# metabolically active developmental stage,
# and exclude it from analysis for now.

outlier = percent_col[percent_col$Plant == 6041,]
# Very low colonization on the side receiving the label,
# but took up a ton.
# I think this fungus was ACTIVELY growing and courting the plant.
# It was new on that side (should have been NM).
# Including it makes it a lot harder to understand anything about
# what is happening with the rest of the plants.
# Exclude for now.


#### How well did root N track myco N? ####

nitrogeninfo = subset(together, compartment_fungus != "None" &
                        compartment_fungus != "MIXED" &
                        compartment_fungus != "OTHER" &
                        received15N == "Y")
looktheseup = (together[together$compartment_fungus == "MIXED" & together$received15N == "Y",])

nitrogeninfo = nitrogeninfo[!is.na(nitrogeninfo$mycorrhizas.APE15N),]
nitrogeninfo = nitrogeninfo[!is.na(nitrogeninfo$uncolonized_roots.APE15N),]

# For log transformation, need predictor axis that is 
# all positive values.
forcefactor = -min(nitrogeninfo$mycorrhizas.APE15N) # small negative: - 0.007
nitrogeninfo$forced.mycorrhizas.APE15N = nitrogeninfo$mycorrhizas.APE15N + forcefactor + 0.000001 # Prevent ratios with zero in denominator

forcefactor_uncolroots = -min(nitrogeninfo$uncolonized_roots.APE15N) # small negative (-0.002), but not as much as mycos.
# Use myco value, then, as the linear transformation.
nitrogeninfo$forced.uncolonized.APE15N = nitrogeninfo$uncolonized_roots.APE15N + forcefactor_uncolroots + 0.000001 # Prevent ratios with zero in denominator

nitrogeninfo_nooutlier = subset(nitrogeninfo,
                                Plant != 6041)


rootNformycoN_plot = ggplot(data = nitrogeninfo_nooutlier) +
  geom_point(aes(x = mycorrhizas.APE15N,
                 y = uncolonized_roots.APE15N, 
                 color = N_level,
                 shape = mycofungus)) +
  geom_smooth(method = "lm",
              aes(x = mycorrhizas.APE15N,
                  y = uncolonized_roots.APE15N),
              color = "black",
              size = 0.5) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  ylab(bquote(atop("Uncolonized roots "^15*N, "(atom percent excess)"))) +
  xlab(expression("Mycorrhizal "^15*"N (atom percent excess)")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")

# STATS FOR PAPER: Root N correlates with mycorrhiza N
rootNformycoN_linear = lm(uncolonized_roots.APE15N ~ mycorrhizas.APE15N, data = nitrogeninfo_nooutlier)
plot(rootNformycoN_linear)
summary(rootNformycoN_linear)

save_plot("plots/Regression_root_N_for_myco_N.pdf", 
          rootNformycoN_plot, 
          base_aspect_ratio = 1.4)

nitrogeninfo[nitrogeninfo$mycorrhizas.APE15N == max(nitrogeninfo$mycorrhizas.APE15N),]
# Plant 6041b is being weird AGAIN. Got a TON of C
# in the hyphae, and a TON of N at the myco.

rootNformycoN_log = glm(forced.uncolonized.APE15N ~ log(forced.mycorrhizas.APE15N), data = nitrogeninfo)
plot(rootNformycoN_log) # very hump-shaped residuals

# INCLUDING OUTLIER:

rootNformycoN_plot_withoutlier = ggplot(data = nitrogeninfo) +
  geom_point(aes(x = mycorrhizas.APE15N,
                 y = uncolonized_roots.APE15N, 
                 color = N_level,
                 shape =compartment_fungus)) +
  geom_smooth(method = "lm",
              aes(x = mycorrhizas.APE15N,
                  y = uncolonized_roots.APE15N),
              color = "black",
              size = 0.5) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  ylab(bquote(atop("Uncolonized roots "^15*N, "(atom percent excess)"))) +
  xlab(expression("Mycorrhizal "^15*"N (atom percent excess)")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")

# STATS FOR PAPER: root N corresponds to myco N including outlier

rootNformycoN_linear_withoutlier = lm(uncolonized_roots.APE15N ~ mycorrhizas.APE15N, data = nitrogeninfo)
plot(rootNformycoN_linear_withoutlier)
summary(rootNformycoN_linear_withoutlier)

rootNformycoN_log_withoutlier = lm(uncolonized_roots.APE15N ~ log(forced.mycorrhizas.APE15N), data = nitrogeninfo)
plot(rootNformycoN_log_withoutlier) # not as good as linear
summary(rootNformycoN_log_withoutlier) # also not as good as linear

save_plot("plots/Regression_root_N_for_myco_N_with_outlier.pdf", 
          rootNformycoN_plot_withoutlier, 
          base_aspect_ratio = 1.4)


#### How well did hyphal C track myco C? ####
carboninfo = subset(together, compartment_fungus != "None" &
                         compartment_fungus != "MIXED" &
                         compartment_fungus != "OTHER")

carboninfo = carboninfo[!is.na(carboninfo$hyphae.APE13C),]
carboninfo = carboninfo[!is.na(carboninfo$mycorrhizas.APE13C),]
# For any analysis involving ONLY C, it makes sense to include
# 1) carbon data from root compartments that
# didn't necessarily receive N label (20 of these, 22 labeled compartments here with hyphal C info), and
# 2) data from plants with failed splits (two of these, 40 successfully split)

justTt = subset(carboninfo, compartment_fungus == "Tt")

carboninfo_nooutlier = carboninfo[!carboninfo$hyphae.APE13C == max(carboninfo$hyphae.APE13C),] # omit outlier 6024b


hyphalCformycoC_plot = ggplot(data = carboninfo_nooutlier) +
  geom_point(aes(x = mycorrhizas.APE13C,
                 y = hyphae.APE13C, 
                 color = N_level,
                 shape =compartment_fungus)) +
  geom_smooth(method = "lm", aes(x = mycorrhizas.APE13C,
                                 y = hyphae.APE13C),
              color = "black",
              size = 0.5) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  ylab(expression("Hyphal "^13*"C (atom percent excess)")) +
  xlab(expression("Mycorrhizal "^13*"C (atom percent excess)")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")

# STATS for supplement: hyphae C tracks myco C, excluding outlier  
myhyphaelm = lm((hyphae.APE13C) ~ mycorrhizas.APE13C, data = carboninfo_nooutlier)
plot(myhyphaelm)
summary(myhyphaelm)

save_plot("plots/Regression_hyphal_C_for_myco_C.pdf", 
          hyphalCformycoC_plot,
          ncol = 1,
          base_aspect_ratio = 1.4)

# Doing this with just Tt yields almost exactly the same result.

# INCLUDING OUTLIER

hyphalCformycoC_plot_withoutlier = ggplot(data = carboninfo) +
  geom_point(aes(x = mycorrhizas.APE13C,
                 y = hyphae.APE13C, 
                 color = N_level,
                 shape =compartment_fungus)) +
  geom_smooth(method = "lm", aes(x = mycorrhizas.APE13C,
                                 y = hyphae.APE13C),
              color = "black",
              size = 0.5) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  ylab(expression("Hyphal "^13*"C (atom percent excess)")) +
  xlab(expression("Mycorrhizal "^13*"C (atom percent excess)")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")


myhyphaelm_withoutlier = lm((hyphae.APE13C) ~ mycorrhizas.APE13C, data = carboninfo)
plot(myhyphaelm_withoutlier)
summary(myhyphaelm_withoutlier)

save_plot("plots/Regression_hyphal_C_for_myco_C_with_outlier.pdf", 
          hyphalCformycoC_plot_withoutlier,
          ncol = 1,
          base_aspect_ratio = 1.4)

#### Making multipanel plot: hyphal C tracks myco C, and root N tracks myco N ####

rootNformycoN_plot_nolegend = rootNformycoN_plot +
  theme(legend.position = "none")

NforNandCforC = plot_grid(rootNformycoN_plot_nolegend, hyphalCformycoC_plot,
          labels = c("A", "B"),
          align = "h",
          rel_widths = c(1, 1.15))

save_plot("plots/SUPP_Multipanel_regressions_NforN_and_CforC.pdf",
          NforNandCforC,
          ncol = 2,
          base_aspect_ratio = 1.4)

rootNformycoN_plot_nolegend_withoutlier = rootNformycoN_plot_withoutlier +
  theme(legend.position = "none")

NforNandCforC_withoutliers = plot_grid(rootNformycoN_plot_nolegend_withoutlier, hyphalCformycoC_plot_withoutlier,
                          labels = c("A", "B"),
                          align = "h",
                          rel_widths = c(1, 1.15))

save_plot("plots/MAIN_Multipanel_regressions_NforN_and_CforC_withoutliers.pdf",
          NforNandCforC_withoutliers,
          ncol = 2,
          base_aspect_ratio = 1.4)

#### C for N regression plot: Do mycorrhizas get more C when they bring back more N? ####

forcurveplot_withoutliers = nitrogeninfo
forcurveplot = subset(nitrogeninfo_nooutlier, Plant != 6024 & Plant != 6041)

# Including outliers

curveplot_withoutliers = ggplot(data = forcurveplot_withoutliers) +
  geom_point(aes(x = forced.mycorrhizas.APE15N,
                 y = mycorrhizas.APE13C, 
                 color = N_level,
                 shape = compartment_fungus)) +
  geom_smooth(method = "lm", 
              formula = y ~ log(x), 
              aes(x = forced.mycorrhizas.APE15N,
                  y = mycorrhizas.APE13C),
              color = "black",
              size = 0.5) +
  ylab(bquote(atop("Mycorrhizal "^13*"C", "(atom percent excess)"))) +
  xlab(bquote(atop("Mycorrhizal "^15*"N", "(atom percent excess, forced positive)"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

save_plot("plots/SUPP_Mycorrhiza_curveplot_withoutliers.pdf",
          curveplot_withoutliers,
          base_aspect_ratio = 1.4)

linear_plot_withoutliers = ggplot(data = forcurveplot_withoutliers) +
  geom_point(aes(x = mycorrhizas.APE15N,
                 y = mycorrhizas.APE13C, 
                 color = N_level,
                 shape = compartment_fungus)) +
  geom_smooth(method = "lm", 
              formula = y ~ x, 
              aes(x = mycorrhizas.APE15N,
                  y = mycorrhizas.APE13C),
              color = "black",
              size = 0.5) +
  ylab(bquote(atop("Mycorrhizal "^13*"C", "(atom percent excess)"))) +
  xlab(bquote(atop("Mycorrhizal "^15*"N (atom percent excess)"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

pdf("plots/Mycorrhiza_linearplot_withoutliers_fortalk.pdf", width = 7, height = 4)
linear_plot_withoutliers
dev.off()

save_plot("plots/MAIN_Mycorrhiza_linearplot_withoutliers.pdf",
          linear_plot_withoutliers,
          base_aspect_ratio = 1.4)

# STATS for supplement: curve fit not as good
logmodel_withoutliers = lm(mycorrhizas.APE13C ~ log(forced.mycorrhizas.APE15N), data = forcurveplot_withoutliers)
plot(logmodel_withoutliers)
summary(logmodel_withoutliers) 

# STATS for main text
linearmodel_withoutliers = lm(mycorrhizas.APE13C ~ (mycorrhizas.APE15N), data = forcurveplot_withoutliers)
plot(linearmodel_withoutliers)
summary(linearmodel_withoutliers)


# Excluding outliers
curveplot = ggplot(data = forcurveplot) +
  geom_point(aes(x = forced.mycorrhizas.APE15N,
                 y = mycorrhizas.APE13C, 
                 color = N_level,
                 shape = compartment_fungus)) +
  geom_smooth(method = "lm", 
              formula = y ~ log(x), 
              aes(x = forced.mycorrhizas.APE15N,
                  y = mycorrhizas.APE13C),
              color = "black",
              size = 0.5) +
  ylab(bquote(atop("Mycorrhizal "^13*"C", "(atom percent excess)"))) +
  xlab(bquote(atop("Mycorrhizal "^15*"N", "(atom percent excess, forced positive)"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
  

save_plot("plots/Mycorrhiza_curveplot.pdf",
          curveplot,
          base_aspect_ratio = 1.4)

linear_plot = ggplot(data = forcurveplot) +
  geom_point(aes(x = mycorrhizas.APE15N,
                 y = mycorrhizas.APE13C, 
                 color = N_level,
                 shape = compartment_fungus)) +
  geom_smooth(method = "lm", 
              formula = y ~ x, 
              aes(x = mycorrhizas.APE15N,
                  y = mycorrhizas.APE13C),
              color = "black",
              size = 0.5) +
  ylab(bquote(atop("Mycorrhizal "^13*"C", "(atom percent excess)"))) +
  xlab(bquote(atop("Mycorrhizal "^15*"N (atom percent excess)"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  theme(plot.margin = unit(c(1,1,0.5,0.5), "cm"))

save_plot("plots/MAIN_Mycorrhiza_linearplot.pdf",
          linear_plot,
          base_aspect_ratio = 1.4)

pdf("plots/Mycorrhiza_linearplot__fortalk.pdf", width = 7, height = 4)
linear_plot
dev.off()

nooutliers_linearlm = lm(mycorrhizas.APE13C ~ mycorrhizas.APE15N, data = forcurveplot)
par(mfrow = c(2,2))
plot(nooutliers_linearlm) # it'll do, I guess.
summary(nooutliers_linearlm)

linear_plot_nolegend = linear_plot +
  theme(legend.position = "none")

linear_and_curve_nooutliers = plot_grid(linear_plot_nolegend, curveplot,
                                        labels = c("A", "B"),
                                        align = "h",
                                        rel_widths = c(1, 1.15))

logmodel= lm(mycorrhizas.APE13C ~ log(forced.mycorrhizas.APE15N), data = forcurveplot)
plot(logmodel)
summary(logmodel) 

linearmodel = lm(mycorrhizas.APE13C ~ (mycorrhizas.APE15N), data = forcurveplot)
plot(linearmodel)
summary(linearmodel)

save_plot("plots/SUPP_Mycorrhizal_multipanel_linear_and_curveplot_nooutliers.pdf",
          linear_and_curve_nooutliers,
          ncol = 2,
          base_aspect_ratio = 1.4)

# FOR SUPP: Versus NM roots, including outliers

mycosvsnm_curve_withoutliers = ggplot(data = forcurveplot_withoutliers) +
  geom_point(aes(x = forced.uncolonized.APE15N,
                 y = mycorrhizas.APE13C, 
                 color = N_level,
                 shape = compartment_fungus)) +
  geom_smooth(method = "lm", 
              formula = y ~ log(x), 
              aes(x = forced.uncolonized.APE15N,
                  y = mycorrhizas.APE13C),
              color = "black",
              size = 0.5) +
  ylab(bquote(atop("Mycorrhizal "^13*"C", "(atom percent excess)"))) +
  xlab(bquote(atop("Uncolonized roots "^15*"N", "(atom percent excess, forced positive)"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))


mycosvsnm_linear_withoutliers = ggplot(data = forcurveplot_withoutliers) +
  geom_point(aes(x = uncolonized_roots.APE15N,
                 y = mycorrhizas.APE13C, 
                 color = N_level,
                 shape = compartment_fungus)) +
  geom_smooth(method = "lm", 
              formula = y ~ x, 
              aes(x = uncolonized_roots.APE15N,
                  y = mycorrhizas.APE13C),
              color = "black",
              size = 0.5) +
  ylab(bquote(atop("Mycorrhizal "^13*"C", "(atom percent excess)"))) +
  xlab(bquote(atop("Uncolonized roots "^15*"N (atom percent excess)"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

mycosvsnm_linear_withoutliers_nolegend = mycosvsnm_linear_withoutliers +
  theme(legend.position = "none")

mycosvsnm_linear_and_curve = plot_grid(mycosvsnm_linear_withoutliers_nolegend, mycosvsnm_curve_withoutliers,
                                        labels = c("A", "B"),
                                        align = "h",
                                        rel_widths = c(1, 1.15))

save_plot("plots/SUPP_MycoC_vs_uncolN_linear_plot.pdf",
          mycosvsnm_linear_withoutliers,
          base_aspect_ratio = 1.4)

mycoc_vs_uncoln_model = lm(mycorrhizas.APE13C ~ uncolonized_roots.APE15N, data = forcurveplot_withoutliers)
plot(mycoc_vs_uncoln_model) # not awesome but will have to do.
summary(mycoc_vs_uncoln_model)

save_plot("plots/Mycovsnm_multipanel_linear_and_curveplot_WITHoutliers.pdf",
          mycosvsnm_linear_and_curve,
          ncol = 2,
          base_aspect_ratio = 1.4)


mycosvsnm_logmodel_withoutliers = lm(mycorrhizas.APE13C ~ log(forced.uncolonized.APE15N), data = forcurveplot_withoutliers)
plot(mycosvsnm_logmodel_withoutliers)
summary(mycosvsnm_logmodel_withoutliers) 

mycosvsnm_linearmodel_withoutliers = lm(mycorrhizas.APE13C ~ uncolonized_roots.APE15N, data = forcurveplot_withoutliers)
plot(mycosvsnm_linearmodel_withoutliers)
summary(mycosvsnm_linearmodel_withoutliers)
# Linear still way better here.

# Hyphae C versus NM roots N (strictest fungal N for plant C comparison), including outliers

hyphaevsnm_curve_withoutliers = ggplot(data = forcurveplot_withoutliers) +
  geom_point(aes(x = forced.uncolonized.APE15N,
                 y = hyphae.APE13C, 
                 color = N_level,
                 shape = compartment_fungus)) +
  geom_smooth(method = "lm", 
              formula = y ~ log(x), 
              aes(x = forced.uncolonized.APE15N,
                  y = hyphae.APE13C),
              color = "black",
              size = 0.5) +
  ylab(bquote(atop("Hyphal "^13*"C", "(atom percent excess)"))) +
  xlab(bquote(atop("Uncolonized roots "^15*"N", "(atom percent excess, forced positive)"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))


hyphaevsnm_linear_withoutliers = ggplot(data = forcurveplot_withoutliers) +
  geom_point(aes(x = uncolonized_roots.APE15N,
                 y = hyphae.APE13C, 
                 color = N_level,
                 shape = compartment_fungus)) +
  geom_smooth(method = "lm", 
              formula = y ~ x, 
              aes(x = uncolonized_roots.APE15N,
                  y = hyphae.APE13C),
              color = "black",
              size = 0.5) +
  ylab(bquote(atop("Mycorrhizal "^13*"C", "(atom percent excess)"))) +
  xlab(bquote(atop("Uncolonized roots "^15*"N (atom percent excess)"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

hyphaevsnm_linear_withoutliers_nolegend = hyphaevsnm_linear_withoutliers +
  theme(legend.position = "none")

hyphaevsnm_linear_and_curve = plot_grid(hyphaevsnm_linear_withoutliers_nolegend, hyphaevsnm_curve_withoutliers,
                                       labels = c("A", "B"),
                                       align = "h",
                                       rel_widths = c(1, 1.15))

save_plot("plots/Hyphaevsnm_multipanel_linear_and_curveplot_WITHoutliers.pdf",
          hyphaevsnm_linear_and_curve,
          ncol = 2,
          base_aspect_ratio = 1.4)


hyphaevsnm_logmodel_withoutliers = lm(hyphae.APE13C ~ log(forced.uncolonized.APE15N), data = forcurveplot_withoutliers)
plot(hyphaevsnm_logmodel_withoutliers)
summary(hyphaevsnm_logmodel_withoutliers) # oof this is terrible 

hyphaevsnm_linearmodel_withoutliers = lm(hyphae.APE13C ~ uncolonized_roots.APE15N, data = forcurveplot_withoutliers)
plot(hyphaevsnm_linearmodel_withoutliers)
summary(hyphaevsnm_linearmodel_withoutliers) # oof also terrible.

# I think maybe my chase periods did not set things up well for 
# strong correspondence between hyphal C and root N.

# parabola_plot = ggplot(data = forcurveplot) +
#   geom_point(aes(x = forced.mycorrhizas.APE15N,
#                  y = mycorrhizas.APE13C, 
#                  color = N_level,
#                  shape = compartment_fungus)) +
#   geom_smooth(method = "lm", 
#               formula = y ~ poly(x, 2), 
#               aes(x = forced.mycorrhizas.APE15N,
#                   y = mycorrhizas.APE13C)) +
#   ylab(bquote(atop("Mycorrhizal "^13*"C", "(atom percent excess)"))) +
#   xlab(bquote(atop("Mycorrhizal "^15*"N (atom percent excess)"))) +
#   scale_color_manual(values = c("steelblue4", "steelblue1"),
#                      name = "N level") +
#   scale_shape_manual(values = c(17, 15),
#                      name = "Fungus") +
#   theme(plot.margin = unit(c(1,1,1,1), "cm"))
# 
# save_plot("plots/Mycorrhiza_parabolaplot.pdf",
#           parabola_plot,
#           base_aspect_ratio = 1.4)

# logmodel = lm(mycorrhizas.APE13C ~ log(forced.mycorrhizas.APE15N), data = forcurveplot)
# plot(logmodel) # Not SUPER great...
# summary(logmodel) # R^2 = 0.17
# 
# linearmodel = lm(mycorrhizas.APE13C ~ (mycorrhizas.APE15N), data = forcurveplot)
# plot(linearmodel) # Also not great.
# summary(linearmodel) # R^2 = 0.22
# 
# parabola_model = lm(mycorrhizas.APE13C ~ poly(forced.mycorrhizas.APE15N, 2), data = forcurveplot)
# plot(parabola_model) # Oh wow, this actually looks way better to me.
# summary(parabola_model) # Also has a lower AIC if I use the glm function
# # R^2 = 0.17


# What if I try omitting SUIPU?

Ttonly = subset(forcurveplot, compartment_fungus == "Tt")

Ttlinear = lm(mycorrhizas.APE13C ~ (forced.mycorrhizas.APE15N), data = Ttonly)
Ttlog = lm(mycorrhizas.APE13C ~ log(forced.mycorrhizas.APE15N), data = Ttonly)

par(mfrow = c(2,2))
plot(Ttlinear)
plot(Ttlog)

plot(linearmodel)
plot(logmodel)

ggplot(data = Ttonly) +
  geom_point(aes(x = forced.mycorrhizas.APE15N,
                 y = mycorrhizas.APE13C, 
                 color = N_level,
                 shape = compartment_fungus)) +
  geom_smooth(method = "lm", 
              formula = y ~ log(x), 
              aes(x = forced.mycorrhizas.APE15N,
                  y = mycorrhizas.APE13C)) +
  ylab(bquote(atop("Mycorrhizal "^13*"C", "(atom percent excess)"))) +
  xlab(bquote(atop("Mycorrhizal "^15*"N (atom percent excess)"))) +
  scale_color_manual(values = c("steelblue4", "steelblue1"),
                     name = "N level") +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))




#### Did the fungi get different exchange rates at diff N levels? ####
labels = c(High = "High N", Low = "Low N")

exchangerates = nitrogeninfo


# exchangerates_onlyhighN = subset(onlyhigh15N, compartment_fungus != "None")
# # This is where the negative APE values no longer work.
# min(exchangerates$mycorrhizas.APE13C) # small positive
# forcefactor = -min(exchangerates$mycorrhizas.APE15N) # small negative
# 
# # How many of these samples have negative values, anyway?
# exchangerates$mycorrhizas.APE15N[exchangerates$mycorrhizas.APE15N <0]
# # Oh my. 21 is a lot.
# 
# exchangerates$waslow = numeric(nrow(exchangerates))
# exchangerates$waslow[exchangerates$mycorrhizas.APE15N < 0] = 1
# 
# 
# exchangerates$forced.mycorrhizas.APE15N = exchangerates$mycorrhizas.APE15N + forcefactor + 0.000001 # Prevent ratios with zero in denominator
# 
exchangerates$forced.mycoC13forN15 = exchangerates$mycorrhizas.APE13C/exchangerates$forced.mycorrhizas.APE15N
# 
# test = exchangerates[!exchangerates$waslow == 1,]

myoutliers = c(6032, 6090) # just outlandishly high C for N ratios; >= 1 order of magnitude greater than other values

toexamine = exchangerates[exchangerates$Plant %in% myoutliers,]

toexamine$mycorrhizas.APE13C
toexamine$mycorrhizas.APE15N # Okay, both of these
# had the same APE value, because (looking back at
# original UCSC data) they had the same d15N: 3.86.


exchangerates_withoutliers = exchangerates

exchangerates_nooutliers = exchangerates[!exchangerates$Plant %in% myoutliers,]

mymax = max(exchangerates_nooutliers$forced.mycoC13forN15)
mythirdquartile = as.numeric(summary(exchangerates_nooutliers$forced.mycoC13forN15)[5])
myIQR = IQR(exchangerates_nooutliers$forced.mycoC13forN15)

# outlierfudgefactor = mythirdquartile + 1.5*myIQR
outlierfudgefactor = mymax


# I know the fungi were getting a lot of C for very little N
# in these outliers.
# The ratio makes them artificially high values.
# Let's coerce them to the maximum value observed in other samples?
exchangerates$forced.mycoC13forN15[exchangerates$Plant %in% myoutliers] = outlierfudgefactor
# exchangerates = exchangerates[!exchangerates$Plant %in% myoutliers,]
# exchangerates = subset(exchangertes, Plant != 6024)

# Oh man, these only are true if I exclude the outliers:
annotations = data.frame(x = c((1:2), (1:2)),
                         y = c(5.5, 4, 4, 5.5),
                         N_level = c(rep("High", 2), rep("Low", 2)),
                         labs = c("a", "b", "b", "ab"))

exchangerate_plot = ggplot(data = exchangerates) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = compartment_fungus, 
                   y = log(forced.mycoC13forN15))) + 
  geom_jitter(aes(x = compartment_fungus,
                 y = log(forced.mycoC13forN15)),
                width = 0.25) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  xlab("Fungus receiving nitrogen label") +
  ylab(bquote(atop("Log exchange rate", "("^13*"C to "^15*"N excess in mycorrhizas)"))) +
  geom_text(data = annotations, aes(x, y, label = labs))
  

exchangerate_test = aov(log(forced.mycoC13forN15) ~ compartment_fungus*N_level, data = exchangerates)
plot(exchangerate_test) #blech

# 25 and 15 look like real outliers.
# These are 6055 and 6031.
# Okay, but removing those, now 18, 22, and 24 are problems
# These are 6043, 6050, 6064 # Haha okay now all the suillus high N are gone.
# Sticking with only removing the crazy outliers
# that were >2 orders of magnitude removed from
# the median rest of the data. (e.g. >1000 ratio values)

# Definitely not quite normal, but anova is pretty robust to that.
exchangerate_test_lm = lm(log(forced.mycoC13forN15) ~ compartment_fungus*N_level, data = exchangerates)
exchangerate_test_orderflipped = aov(log(forced.mycoC13forN15) ~ N_level*compartment_fungus, data = exchangerates)

summary(exchangerate_test) # significant interaction.
summary(exchangerate_test_orderflipped) # Same result as when I put comp fungus first. Interaction sinificant.
summary(exchangerate_test_lm) # everything's significant if I just use lm

TukeyHSD(exchangerate_test)

median(exchangerates$forced.mycoC13forN15[exchangerates$compartment_fungus == "Sp" & exchangerates$N_level == "High"])
median(exchangerates$forced.mycoC13forN15[exchangerates$compartment_fungus == "Tt" & exchangerates$N_level == "High"])


exchangerate_tukey = TukeyHSD(exchangerate_test)
write.csv(exchangerate_tukey$N_level, "Statistical_tables/exchangerate_Tukey_output_Nlevel.csv")
write.csv(exchangerate_tukey$compartment_fungus, "Statistical_tables/exchangerate_Tukey_output_fungus.csv")
write.csv(exchangerate_tukey$`compartment_fungus:N_level`, "Statistical_tables/exchangerate_Tukey_output_Nlevel-Fungus.csv")

# pdf("plots/exchange_rate_boxplot_nooutliers.pdf", width = 7, height = 5)
# exchangerate_plot
# dev.off()

save_plot("plots/MAIN_Exchange_rate_boxplot_by_fungus.pdf", exchangerate_plot)


pdf("plots/exchange_rate_boxplot_fortalk.pdf", width = 8, height = 4)
exchangerate_plot
dev.off()

save_plot("plots/MAIN_Exchange_rate_boxplot_by_fungus.pdf", 
          exchangerate_plot,
          base_aspect_ratio = 1.4)


tx = with(exchangerates, interaction(N_level, compartment_fungus))

forlabels = aov(log(forced.mycoC13forN15) ~ tx, data = exchangerates)
mylabels = HSD.test(forlabels, "tx", group = TRUE)

#### Exchange rates by competition treatment ####

exchangerates$competition_treatment = numeric(nrow(exchangerates))
for (i in 1:nrow(exchangerates)) {
  if (exchangerates$Fungi[i] == "Sp/None" |
      exchangerates$Fungi[i] == "Tt/None" |
      exchangerates$Fungi[i] == "THETE") {
      exchangerates$competition_treatment[i] = "None"
  } else if (exchangerates$Fungi[i] == "Sp/Sp" | exchangerates$Fungi[i] == "Tt/Tt") {
      exchangerates$competition_treatment[i] = "Self"
  } else if (exchangerates$Fungi[i] == "MIXED/THETE" | 
             exchangerates$Fungi[i] == "MIXED/SUIPU" |
             exchangerates$Fungi[i] == "Tt/Sp") {
      exchangerates$competition_treatment[i] = "Other"
  }
}

onlyTt = subset(exchangerates, compartment_fungus == "Tt")

competition_exchangerate_plot = ggplot(data = exchangerates) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = competition_treatment,
                   y = log(forced.mycoC13forN15))) +
  geom_jitter(width = 0.20,
              aes(x = competition_treatment,
                  y = log(forced.mycoC13forN15),
                 shape = compartment_fungus)) +
  scale_shape_manual(values = c(17, 15),
                     name = "Fungus") +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  xlab("Competition treatment") +
  ylab(bquote(atop("Log exchange rate", "("^13*"C to "^15*"N excess in mycorrhizas)")))

save_plot("plots/MAIN_Competition_and_N_treatment_exchange_rate_plot.pdf",
          competition_exchangerate_plot,
          base_aspect_ratio = 1.4)

exchangeratebycomp = aov(log(forced.mycoC13forN15) ~ competition_treatment*N_level, data = exchangerates)
plot(exchangeratebycomp)
summary(exchangeratebycomp)

allfactorsexchangerate = lm(log(forced.mycoC13forN15) ~ competition_treatment*N_level*compartment_fungus, data = exchangerates)
plot(allfactorsexchangerate)
summary(allfactorsexchangerate)

justhighN = subset(exchangerates, N_level == "High")
Nanova = aov(log(forced.mycoC13forN15) ~ competition_treatment, data = justhighN)
plot(Nanova)
summary(Nanova)

onlyTt_selfvsother = subset(onlyTt, competition_treatment != "None")
onlyTt_selfvsother$competition_treatment = as.factor(onlyTt_selfvsother$competition_treatment)

require(ggsignif)

seginfo = data.frame(x = 0.9, 
                         xend = 2.1,
                         y = 5.3, 
                         yend = 5.3, 
                         color = "black",
                         N_level = "High")

starinfo = data.frame(x = c((1.5), (1.5)),
                         y = c(5.7, 5.7),
                         N_level = c(rep("High", 1), rep("Low", 1)),
                         labs = c(".", ""))



onlyTt_selfvsother_plot = ggplot(data = onlyTt_selfvsother) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = competition_treatment,
                   y = log(forced.mycoC13forN15))) +
  geom_jitter(width = 0.25,
              aes(x = competition_treatment,
                 y = log(forced.mycoC13forN15))) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  geom_segment(data=seginfo,
           aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE) +
  geom_text(data = starinfo, aes(x, y, label = labs), size = 8) +
  xlab("Competition treatment") +
  theme(plot.margin = unit(c(1,1,0.5,0.5), "cm")) +
  ylab(bquote(atop("Log exchange rate: Tt", "("^13*"C to "^15*"N excess in mycorrhizas)")))

pdf("plots/MAIN_Thelephora_self_vs_other_t_test_plot_fortalk.pdf", width = 8, height = 4)
onlyTt_selfvsother_plot
dev.off()

save_plot("plots/MAIN_Thelephora_self_vs_other_t_test_plot_forpaper.pdf",
          onlyTt_selfvsother_plot,
          base_aspect_ratio = 1.8)

median(onlyTt_selfvsother$forced.mycoC13forN15[onlyTt_selfvsother$competition_treatment == "Self" & only_Tt_selfvsother$N_level == "High"])
median(onlyTt_selfvsother$forced.mycoC13forN15[onlyTt_selfvsother$competition_treatment == "Other" & only_Tt_selfvsother$N_level == "High"])


exchangeratebycomp_onlytt = aov(log(forced.mycoC13forN15) ~ competition_treatment*N_level, data = onlyTt)
summary(exchangeratebycomp_onlytt)
# None of this is even a little bit significant with a simple ANOVA.
# But, to be fair, I did have a clear a priori hypothesis.

self_high = subset(onlyTt, competition_treatment == "Self" & N_level == "High")
other_high = subset(onlyTt, competition_treatment == "Other" & N_level == "High")

mean(log(self_high$forced.mycoC13forN15))
sd(log(self_high$forced.mycoC13forN15))

mean(log(other_high$forced.mycoC13forN15))
sd(log(other_high$forced.mycoC13forN15))


t.test(log(self_high$forced.mycoC13forN15), log(other_high$forced.mycoC13forN15))
# This appears to be the only way to pull out a "significant" difference
# in exchange rate depending on competitive setting.
self_low = subset(onlyTt, competition_treatment == "Self" & N_level == "Low")
other_low = subset(onlyTt, competition_treatment == "Other" & N_level == "Low")
t.test(log(self_low$forced.mycoC13forN15), log(other_low$forced.mycoC13forN15))

mean(log(self_low$forced.mycoC13forN15))
sd(log(self_low$forced.mycoC13forN15))

mean(log(other_low$forced.mycoC13forN15))
sd(log(other_low$forced.mycoC13forN15))



#### Looking at exchange rates for ONLY plants that DEFINITELY took up label ####

# Okay, here we only have high N Tt plants.
onlyhighN_plot = ggplot(data = exchangerates_onlyhighN) +
  geom_boxplot(aes(x = Fungi, 
                   y = (mycoC13forN15))) + 
  geom_point(aes(x = Fungi,
                 y = (mycoC13forN15))) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  xlab("Fungal treatment") +
  ylab(expression("Ratio of "^13*"C to "^15*"N excess in Tt mycorrhizas"))

onlyhighN_model = aov(mycoC13forN15 ~ Fungi, data = exchangerates_onlyhighN)
plot(onlyhighN_model)
summary(onlyhighN_model)
# Not significant.

t.test(onlyhigh15N$mycoC13forN15[onlyhigh15N$Fungi == "Tt/Sp"], onlyhigh15N$mycoC13forN15[onlyhigh15N$Fungi == "Tt/Tt"])
# Still NS, though would be marginal if did a one-tailed test.

pdf("plots/boxplot_C_for_N_by_N_level.pdf", width = 7, height = 5)
CforNbyN_level
dev.off()

noNM = subset(together, compartment_fungus != "NM")

exchangeratebyn = glm(mycoC13forN15 ~ N_level*compartment_fungus*percent_col, data = noNM)
plot(exchangeratebyn) # Not great
summary(exchangeratebyn)
# percent_col, N level: fungus, N level: percent_col,
# fungus: percent col, and fungus: N level: percent col.

test = noNM[is.na(noNM$mycoC13forN15),]

ggplot(data = noNM) +
  geom_point(aes(x = percent_col,
             y = mycoC13forN15,
             color = compartment_fungus,
             shape = N_level))

# WOW that is not informative because there
# are a couple huge outliers.

noNM[noNM$mycoC13forN15 == max(noNM$mycoC13forN15),] # 6050b
noNM[noNM$mycoC13forN15 == min(noNM$mycoC13forN15),] # 6029b

noNM_nooutliers = subset(noNM, Plant != 6050 & Plant != 6029)

ggplot(data = noNM_nooutliers) +
  geom_point(aes(x = percent_col,
                 y = mycoC13forN15,
                 color = compartment_fungus,
                 shape = N_level))

ggplot(data = noNM) +
  geom_point(aes(x = percent_col,
                 y = uncolonized_roots.APE15N,
                 color = compartment_fungus,
                 shape = N_level))

# Makes no difference how colonized you were, got same amount of N

ggplot(data = noNM) +
  geom_point(aes(x = percent_col,
                 y = mycorrhizas.APE13C,
                 color = compartment_fungus,
                 shape = N_level))

# Also does not affect concentration of 13C in mycorrhizas

# Note that we are considering percent colonization at the compartment
# level right now.

#### How does the allocation ratio to Tt change with competition? ####
# Hypothesis: Relative allocation to Tt should
# be greater when Tt is vs. Sp or vs. None.
# When vs. Tt, log alloc ratio should be zero.

allocratios = together[-grep("MIXED", together$Fungi),]
allocratios = allocratios[-grep("FAILED", together$Fungi),]

allocratios$allocratio = numeric(nrow(allocratios))
smallconstant = 0.000473 # Need to add this to all C values to get around a negative value below.
allocratios$mycorrhizas.APE13C = allocratios$mycorrhizas.APE13C + smallconstant
allocratios$competition_treatment[grep("OTHER", allocratios$competition_treatment)] = "Other"
allocratios$competition_treatment[grep("SELF", allocratios$competition_treatment)] = "Self"
allocratios$competition_treatment[grep("NM", allocratios$competition_treatment)] = "None"


for (i in 1:nrow(allocratios)) {
  for (j in 1:nrow(allocratios)) {
    if (allocratios$Plant[i] == allocratios$Plant[j] &
        allocratios$Side[i] != allocratios$Side[j]) {
      if (allocratios$compartment_fungus[i] == "Tt") {
        if (allocratios$compartment_fungus[j] != "None") {
          allocratios$allocratio[i] = allocratios$mycorrhizas.APE13C[i]/allocratios$mycorrhizas.APE13C[j]
          allocratios$allocratio[j] = NA # to ensure one entry per plant
        } else if (allocratios$compartment_fungus[j] == "None") {
          allocratios$allocratio[i] = allocratios$mycorrhizas.APE13C[i]/allocratios$uncolonized_roots.APE13C[j]
        }
      } else if (allocratios$compartment_fungus[i] == "Sp") {
        if (allocratios$compartment_fungus[j] == "Sp") {
          allocratios$allocratio[i] = allocratios$mycorrhizas.APE13C[i]/allocratios$mycorrhizas.APE13C[j]
          allocratios$allocratio[j] = NA # to ensure one entry per plant
        } else if (allocratios$compartment_fungus[j] == "None") {
          allocratios$allocratio[i] = allocratios$mycorrhizas.APE13C[i]/allocratios$uncolonized_roots.APE13C[j]
        }
      }
    }
  }
}

spalloc = allocratios[grep("Sp", allocratios$compartment_fungus),]
spalloc = spalloc[!is.na(spalloc$allocratio),]
spalloc$logallocratio = log(spalloc$allocratio + 0.07)

allocratios = allocratios[grep("Tt", allocratios$compartment_fungus),]
allocratios = allocratios[!is.na(allocratios$allocratio),]
# Maybe this is working?

# allocratios$logallocratio = log(allocratios$allocratio + 0.07) 
# whoops though! also have a negative value.
# It's a Tt/Tt plant. 6073. I guess at least
# one side wasn't really enriched for C?
# just6073 = together[together$Plant == 6073,]
# That value is -0.0004728348.
# Let's just add 0.000473 to everything upstream.

# minimum allocation ratio observed, besides 0, is 0.1428159.
# I'll add 0.07, half this value, to everything that is zero so it can be logged.
allocratios$logallocratio = log(allocratios$allocratio + 0.07) 

labels = c(High = "High N", Low = "Low N")

nonehigh = subset(allocratios, competition_treatment == "None" & N_level == "High")
t.test(nonehigh$logallocratio) # marginal: p = 0.06
otherhigh = subset(allocratios, competition_treatment == "Other" & N_level == "High")
t.test(otherhigh$logallocratio) # marginal: p = 0.09
otherlow = subset(allocratios, competition_treatment == "Other" & N_level == "Low")
t.test(otherlow$logallocratio) # significant: p = 0.02


ggplot(data = allocratios) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = competition_treatment, y = logallocratio)) +
  geom_jitter(aes(x = competition_treatment, y = logallocratio),
              width = 0.25) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels))

ggplot(data = spalloc) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = competition_treatment, y = logallocratio)) +
  geom_jitter(aes(x = competition_treatment, y = logallocratio),
              width = 0.25) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels))

# VERY interesting.
# Tt gets more C relative to Sp under high N;
# basically same amount vs itself.
# Tt gets LESS N relative to Sp under low N,
# basically same amount vs itself.

justvsmycos = subset(allocratios, competition_treatment != "None")

# STATS for main text
t.test(logallocratio ~ competition_treatment, data = subset(justvsmycos, N_level == "High"))
t.test(logallocratio ~ competition_treatment, data = subset(justvsmycos, N_level == "Low"))

justvsmycos %>% group_by(N_level, competition_treatment) %>% summarize(mu = mean(logallocratio), stdev = sd(logallocratio), med = median(logallocratio))


anovaoverall = aov(logallocratio ~ competition_treatment*N_level, data = justvsmycos)
summary(anovaoverall)
# Significant interaction, not either thing by itself.
anovaflipped = aov(logallocratio ~ N_level*competition_treatment, data = justvsmycos)
summary(anovaflipped) # Basically same result, phew!

TukeyHSD(anovaoverall)
# Other low vs other high
# self high vs other high

figfortalk = ggplot(data = justvsmycos) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = competition_treatment, y = logallocratio)) +
  geom_jitter(aes(x = competition_treatment, y = logallocratio),
              width = 0.25) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  xlab("Competition treatment") +
  ylab(bquote(atop("Log exchange rate", "("^13*"C to "^15*"N excess in mycorrhizas)")))

pdf("plots/Allocation_ratio_boxplot_fortalk.pdf", width = 8, height = 4)
figfortalk
dev.off()

median(justvsmycos$allocratio[justvsmycos$competition_treatment == "Other" &
                              justvsmycos$N_level == "High"])
# 1.904
median(justvsmycos$allocratio[justvsmycos$competition_treatment == "Self" &
                              justvsmycos$N_level == "High"])
# 1.0887

median(justvsmycos$allocratio[justvsmycos$competition_treatment == "Other" &
                              justvsmycos$N_level == "Low"])
# 0.4771 -- got 48% as much C as Suillus. 

median(justvsmycos$allocratio[justvsmycos$competition_treatment == "Self" &
                                justvsmycos$N_level == "Low"])
# 1.0814





starinfo = data.frame(x = c((1.5), (1.5)),
                      y = c(2.4, 2.1),
                      N_level = c(rep("High", 1), rep("Low", 1)),
                      labs = c(".", "*"))

seginfo = data.frame(x = c(0.9, 0.9), 
                     xend = c(2.1, 2.1),
                     y = c(2, 2), 
                     yend = c(2, 2), 
                     color = c("black", "black"),
                     N_level = c("High", "Low"))


figforpaper = ggplot(data = justvsmycos) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = competition_treatment, y = logallocratio)) +
  geom_jitter(aes(x = competition_treatment, y = logallocratio),
              width = 0.25) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  geom_segment(data=seginfo,
               aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE) +
  geom_text(data = starinfo, aes(x, y, label = labs), size = 8) +
  theme(plot.margin = unit(c(1,1,0.5,0.5), "cm")) +
  xlab("Competition treatment") +
  ylab(expression(atop("Log carbon allocation ratio","("*italic("Thelephora terrestris"*")"))))

save_plot("plots/MAIN_Allocation_ratio_for_paper.pdf",
          figforpaper,
          base_aspect_ratio = 1.8)

pdf("plots/Allocation_ratio_boxplot_withannots_fortalk.pdf", width = 8, height = 4)
figforpaper
dev.off()

#### Were rhizomorphs from low N treatments more C-13 enriched? ####

nitrogeninfo[nitrogeninfo$N_level == "Low",]


ggplot(data = nitrogeninfo) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = N_level, y = hyphae.APE13C)) +
  geom_jitter(aes(x = N_level, y = hyphae.APE13C),
              width = 0.25)

ggplot(data = nitrogeninfo) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = N_level, y = mycorrhizas.APE13C)) +
  geom_jitter(aes(x = N_level, y = mycorrhizas.APE13C),
              width = 0.25)

ggplot(data = nitrogeninfo) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = N_level, y = hyphae.APE15N)) +
  geom_jitter(aes(x = N_level, y = hyphae.APE15N),
              width = 0.25)

#### How did plant C:N and fungal C:N compare? ####

fortissueplot = minimally_processed_isotopes %>% group_by(Plant, tissue)
ggplot(data = fortissueplot) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = tissue, y = CNratio)) +
  geom_jitter(aes(x = tissue, y = CNratio),
              width = 0.25)

fortissueplot %>% group_by(tissue) %>% summarize(meanCN = mean(CNratio), sdCN = sd(CNratio))

byfungus = minimally_processed_isotopes %>% group_by(Plant, compartment_fungus, tissue)

ggplot(data = subset(byfungus, compartment_fungus != "MIXED" & compartment_fungus != "OTHER")) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = tissue, y = CNratio, color = compartment_fungus)) +
  geom_jitter(aes(x = tissue, y = CNratio, color = compartment_fungus),
              width = 0.25) 

# Yep, plain ol' roots have much higher C:N than hyphae or mycos!

ggplot(data = subset(fortissueplot, tissue == "mycorrhizas" & enriched == 1 & `Receives 15N label?` == "Y")) +
         geom_point(aes(x = pctN, y = pctC))
test = lm(pctC ~ pctN, data = subset(fortissueplot, tissue == "mycorrhizas" & enriched == 1 & `Receives 15N label?` == "Y"))
# k great, this almost over-dispersed scatter plot indeed
# represents no relationship.

ggplot(data = subset(fortissueplot, tissue == "mycorrhizas" & enriched == 1 & `Receives 15N label?` == "Y")) +
  geom_point(aes(x = atmpct15N, y = atmpct13C))

test2 = lm(APE13C ~ APE15N, data = subset(fortissueplot, tissue == "mycorrhizas" & enriched == 1 & `Receives 15N label?` == "Y"))
summary(test2)
# And, indeed, this is still a good relationship
# even when I haven't subsetted to
# exclude mixed cultures.

ggplot(data = subset(fortissueplot, tissue == "mycorrhizas" & enriched == 1 & `Receives 15N label?` == "Y")) +
  geom_point(aes(x = APE15N, y = APE13C))

ggplot(data = subset(fortissueplot, tissue == "mycorrhizas" & enriched == 1)) +
  geom_point(aes(x = CNratio, y = APE13C))
# Less fungi (higher CN) seems to mean less carbon from plant.

ggplot(data = subset(fortissueplot, tissue == "mycorrhizas" & enriched == 1 & `Receives 15N label?` == "Y")) +
  geom_point(aes(x = CNratio, y = APE15N))
# Basically no relationship here. Maybe more N label
# in the mycos with more fungi (lower CN)

ggplot(data = subset(fortissueplot, tissue == "hyphae" & enriched == 1 & `Receives 15N label?` == "Y")) +
  geom_point(aes(x = pctN, y = pctC))
# more consistent CN ratio for hyphae

ggplot(data = subset(fortissueplot, tissue == "hyphae" & enriched == 1 & `Receives 15N label?` == "Y")) +
  geom_point(aes(x = APE15N, y = APE13C))
# Generally positive.

ggplot(data = subset(fortissueplot, tissue == "uncolonized_roots" & enriched == 1 & `Receives 15N label?` == "Y")) +
  geom_point(aes(x = pctN, y = pctC))
# Very inconsistent CN ratio here for uncolonized roots.
# Generally 47% C or so, but N % varies a lot (1-3)

ggplot(data = subset(fortissueplot, tissue == "uncolonized_roots" & enriched == 1 & `Receives 15N label?` == "Y")) +
  geom_point(aes(x = APE15N, y = APE13C))
# More spread, generally positive?

#### Did Sp get more 13C than Tt, generally? ####
justmycos = select(together, Plant, Side, mycorrhizas.APE13C, mycofungus, N_level:enriched)
justmycos = subset(justmycos, mycofungus != "OTHER" & 
                     compartment_fungus != "MIXED" &
                     mycofungus != "None" &
                     is.na(mycorrhizas.APE13C) == FALSE)
justmycos %>% group_by(mycofungus) %>% summarize(mean.13C = mean(mycorrhizas.APE13C), sd.13C = sd(mycorrhizas.APE13C))

# Yes.

ggplot(justmycos) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = mycofungus, y = mycorrhizas.APE13C)) +
  geom_jitter(aes(x = mycofungus, y = mycorrhizas.APE13C)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels))

# But this plot is a bit pseudoreplicated, because
# any plant that had two +EMF compartments is shown twice.
# Can I fix this?

getridpseudo = subset(justmycos, competition_treatment == "NM" | competition_treatment == "SELF")
getridpseudo$avgmyco13C = numeric(nrow(getridpseudo))
for (i in 1:nrow(getridpseudo)) {
  for (j in 1:nrow(getridpseudo)) {
    if (getridpseudo$competition_treatment[i] != "NM" &
        getridpseudo$Plant[i] == getridpseudo$Plant[j] &
        getridpseudo$Side[i] != getridpseudo$Side[j]) {
      getridpseudo$avgmyco13C[i] = mean(getridpseudo$mycorrhizas.APE13C[i], getridpseudo$mycorrhizas.APE13C[j])
    } else if (getridpseudo$competition_treatment[i] == "NM") {
      getridpseudo$avgmyco13C[i] = getridpseudo$mycorrhizas.APE13C[i]
    }
  }
}
getridpseudo = getridpseudo[!duplicated(getridpseudo$Plant),]

ggplot(getridpseudo) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = mycofungus, y = avgmyco13C)) +
  geom_jitter(aes(x = mycofungus, y = avgmyco13C)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels))

  
