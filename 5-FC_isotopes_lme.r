# Trying to analyze FC isotopes in a more cohesive way
# Attempting a linear mixed model a la Argüello et al. 2016
# August 2019

# This tutorial may be helpful: https://ourcodingclub.github.io/2017/03/15/mixed-models.html
# This is also good: http://www.bodowinter.com/tutorial/bw_LME_tutorial.pdf
# Hmm... https://stats.stackexchange.com/questions/35071/what-is-rank-deficiency-and-how-to-deal-with-it

# My model is consistently rank deficient.
# This seems unnecessarily annoying, given
# that the patterns I can see appear to be reasonably clear.

# Consider Bayesian approach? https://stats.stackexchange.com/questions/368272/why-does-logistic-regression-doesnt-give-the-same-result-as-a-chi-square-test-w

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

together$mycoC13ppmexcess = together$mycorrhizas.APE13C * (10^4)

together = together[together$enriched != 0,]

min(together$mycoC13ppmexcess[!is.na(together$mycoC13ppmexcess)])

together$transmycoC13 = (log(together$mycoC13ppmexcess))

# Shouldn't model mycorrhiza-specific phenomena with NM plants

nonm = together[!is.na(together$mycorrhizas.APE13C),]
nonm = subset(nonm, compartment_fungus != "None")
# Should I exclude microcosms with mixed cultures?

excluding_mixed = nonm[-grep("MIXED", nonm$competitors),]
# excluding_mixed$versus = relevel(excluding_mixed$versus, levels = c("None", "Sp", "Tt"))

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
# Ideally, I'd have a random slope for each,
# but a random intercept is all I can do with my dataset
# and is probably reasonable.

# Let's try an LRT

## trying with transformed data
c13.null = lmer(transmycoC13 ~ compartment_fungus +(1|Plant) + (1|Batch), 
                REML = FALSE,
                data = excluding_mixed)
c13.withversus = lmer(transmycoC13 ~ compartment_fungus + versus +(1|Plant) + (1|Batch), 
                      REML = FALSE,
                      data = excluding_mixed)
c13.withversusandint = lmer(transmycoC13 ~ compartment_fungus * versus +(1|Plant) + (1|Batch), 
                      REML = FALSE,
                      data = excluding_mixed)
c13.vsandintandn = lmer(transmycoC13 ~ compartment_fungus * versus + N_level +(1|Plant) + (1|Batch), 
                      REML = FALSE,
                      data = excluding_mixed)

c13.versusandnint = lmer(transmycoC13 ~ compartment_fungus * versus * N_level +(1|Plant) + (1|Batch), 
                      REML = FALSE,
                      data = excluding_mixed)
c13.versusandnint = lmer(transmycoC13 ~ compartment_fungus * versus * N_level + (1|Batch) + (1|Plant), 
                         REML = FALSE,
                         data = excluding_mixed)
# Okay actually even the above is working now.
# I think my consistent failures before were
# due to some weird problem with data filtering.
# Hooray?

anova(c13.null, c13.withversus)
anova(c13.withversus, c13.versusandn)

# Trying random slopes

c13.null = lmer(mycoC13ppmexcess ~ compartment_fungus + (mycoC13ppmexcess|Batch), 
                REML = FALSE,
                data = excluding_mixed)
# Error in eval_f(x, ...) : Downdated VtV is not positive definite

#### Bayesian? ####
require(blme)

glmb = blmer(transmycoC13 ~ compartment_fungus + (1|Batch), data=excluding_mixed, 
              fixef.prior = normal(cov = diag(9,3)))
summary(glmb)

glmb = blmer(transmycoC13 ~ compartment_fungus * versus * N_level + (1|Batch), data=excluding_mixed, 
             fixef.prior = normal(cov = diag(9,3)))
