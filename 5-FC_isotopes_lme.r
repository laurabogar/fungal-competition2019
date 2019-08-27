# Trying to analyze FC isotopes in a more cohesive way
# Attempting a linear mixed model a la Argüello et al. 2016
# August 2019

# This tutorial may be helpful: https://ourcodingclub.github.io/2017/03/15/mixed-models.html
# This is also good: http://www.bodowinter.com/tutorial/bw_LME_tutorial.pdf

# Libraries needed:
require(tidyverse)
require(lme4)

# Data:
together = read_csv("./FCdata/isotopes_together_as_analyzed.csv")
metadata = read_csv("./FCdata/Fungal_competition_plant_tracking.csv")

batchtomerge = select(metadata_byplant, Plant, Batch)

together = left_join(together, batchtomerge)

# Filtering

together$Batch = as.factor(together$Batch)

together = together[-grep("FAILED", together$competitors),]

together$versus = numeric(nrow(together))

for (i in 1:nrow(together)) {
  if (together$compartment_fungus[i] == "Sp") {
    if (together$competitors[i] == "SUIPU/NM") {
      together$versus[i] = "None"
    } else if (together$competitors[i] == "SUIPU/SUIPU") {
      together$versus[i] = "Sp"
    } else if (together$competitors[i] == "THETE/SUIPU") {
      together$versus[i] = "Tt"
    } else if (grepl("MIXED", together$competitors[i])){
      together$versus[i] = "Mixed"
    }
  } else if (together$compartment_fungus[i] == "Tt") {
    if (together$competitors[i] == "THETE/NM") {
      together$versus[i] = "None"
    } else if (together$competitors[i] == "THETE/THETE") {
      together$versus[i] = "Tt"
    } else if (together$competitors[i] == "THETE/SUIPU") {
      together$versus[i] = "Sp"
    } else if (grepl("MIXED", together$competitors[i])) {
      together$versus[i] = "Mixed"
    }
  } else if (together$compartment_fungus[i] == "None") {
    if (together$competitors[i] == "SUIPU/NM") {
      together$versus[i] = "Sp"
    } else if (together$competitors[i] == "NM/NM") {
      together$versus[i] = "None"
    } else if (together$competitors[i] == "THETE/NM") {
      together$versus[i] = "Tt"
    }
  } else if (together$compartment_fungus[i] == "MIXED") {
      if (together$competitors[i] == "MIXED/SUIPU") {
        together$versus[i] = "Sp"
      } else if (together$competitors[i] == "MIXED/THETE") {
        together$versus[i] = "Tt"
     }
  }
}

# Shouldn't model mycorrhiza-specific phenomena with NM plants

nonm = together[!together$compartment_fungus == "None",]

# Should I exclude microcosms with mixed cultures?

excluding_mixed = together[-grep("MIXED", together$competitors),]

#### C-13 enrichment of mycos by species ####
# Does the C-13 enrichment of mycorrhizas depend on the species
# of fungus forming the mycorrhiza, controlling for competitor
# identity?

# Argüello model: 
# log(14C in hyphal compartment A) ~ 
# AMF fungus side A*AMF fungus side B + 
# plant species identity + random pot effect
# with C labeling group "as a covariate." Does this mean as a random effect?

# Plant ID and C-13 labeling batch should  be random effects here.
mymodel = lmer((mycorrhizas.APE13C) ~ compartment_fungus + versus + N_level + (1|Plant) + (1|Batch), data = nonm)

mymodel = lmer((mycorrhizas.APE13C) ~ compartment_fungus * versus * N_level + (1|Plant) + (1|Batch), data = nonm)
summary(mymodel)

excludingmixedmodel = lmer((mycorrhizas.APE13C) ~ compartment_fungus * versus * N_level + (1|Plant) + (1|Batch), data = excluding_mixed)
