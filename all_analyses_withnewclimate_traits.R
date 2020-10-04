################################
##Trait / Environment Analyses##
################################

setwd("~/Google Drive/ProteaPaper/Andrew_analyses_makeover")

library(dplyr)
library(lme4)
library(lmerTest)
library(R.utils)
library(fields)
 
# load trait data with Rachel's MDS values included
d = read.csv("repens_drydown_alldata_mds.csv")
# get rid of junk rows
d = d[,-grep("X.", names(d))]
# get rid of rows with no population value (not measured)
d = d[!is.na(d$pop),]

# load environmental data
d1_pxt = read.csv("../Melis_analyses_makeover/voomd1_popxtreatGenesTraits.csv", header=TRUE)
envt_data = d1_pxt[,c("pop", "elevation", "MAP", "JanmaxT", "JulyminT", "summerRain", "rain_interannual_variability", "TempRange")]
envt_data = envt_data[match(unique(envt_data$pop), envt_data$pop),]

# merge trait and environmental data 
d_envt = merge(d, envt_data, by="pop")
dim(d_envt)
  
# Separate data for sampling time 1 (day 6) and sampling time 2 (day 12)
d1 = filter(d_envt, round==1)
d2 = filter(d_envt, round==2)



##############################
## with treatments TOGETHER ##


###################################################
# comparing traits to environmental variables #
###################################################

# Function to fit combinations of traits and environment, and return the values for Tables S6 and S7 
fit_envt_trait_models <- function(traits, envt_vars, trt, pop, MaternalLine, PC, use.PC, use.fdrtool) {
  # traits is the traits dataset, envt_vars is the environment data set, pop is population codes, MaternalLine is maternal line codes, trt is treatment level.
  # use.PC indicates whether to add PC/MDS values for neutral genetic structure
  # use.fdrtool indicates whether to use fdrtool to get qvalues (default is p.adjust)
  trait_labels = envt_labels = df_envt = df_trt = df_int = chisq_envt = chisq_trt = chisq_int = pval_envt = pval_trt = pval_int = {} # storage variables
for (i in 1:ncol(traits)) { # look through trait response variables
  x = scale(traits[,i])
  for (j in 1:ncol(envt_vars)) { # loop through environmental variables
    y = scale(envt_vars[,j])
    dtemp = data.frame(trait=x, envt=y, trt=trt, pop = pop, MaternalLine=MaternalLine, PC=PC) # temporary data fram for analysis containing this iteration's trait & envt
    trait_labels = c(trait_labels, names(traits)[i])
    envt_labels = c(envt_labels, names(envt_vars)[j])
    
    # Fit the models that include nested random effects for Maternal Line within population, but not the genetic structure PC
    if (use.PC==FALSE) {
      # remove any rows containing NAs in y or x
      missings = which(is.na(y) | is.na(x))
      if (length(missings)>0) dtemp = dtemp[-missings,]
      # fit the full and reduced models
      fit_full = lmer(trait ~ envt*trt + (1 | pop/MaternalLine), data=dtemp, REML=FALSE)
      fit_no_envt = lmer(trait ~ trt + (1 | pop/MaternalLine), data=dtemp, REML=FALSE)
      fit_no_trt = lmer(trait ~ envt + (1 | pop/MaternalLine), data=dtemp, REML=FALSE)
      fit_no_int = lmer(trait ~ envt + trt + (1 | pop/MaternalLine), data=dtemp, REML=FALSE)
    }
    
    # Alternatively, fit the models that include the genetic structure PC as well as nested random effects
    if (use.PC==TRUE) {
      # remove any rows containing NAs in y or x or the genetic structure PC
      missings = which(is.na(y) | is.na(x) | is.na(PC))
      if (length(missings)>0) dtemp = dtemp[-missings,]
      # fit the full and reduced models
      fit_full = lmer(trait ~ envt*trt + PC + (1 | pop), data=dtemp, REML=FALSE)
      fit_no_envt = lmer(trait ~ trt + PC + (1 | pop), data=dtemp, REML=FALSE)
      fit_no_trt = lmer(trait ~ envt + PC + (1 | pop), data=dtemp, REML=FALSE)
      fit_no_int = lmer(trait ~ envt + trt + PC + (1 | pop), data=dtemp, REML=FALSE)
    }
    
    # Model tests: full versus removing one variable at a time (LRT using Chisq)
    anova_no_envt = anova(fit_no_envt, fit_full)
    anova_no_trt = anova(fit_no_trt, fit_full)
    anova_no_int = anova(fit_no_int, fit_full)
    # Pull out values from the model comparisons for the results table
    df_envt = c(df_envt, anova_no_envt$`Chi Df`[2])
    df_trt = c(df_trt, anova_no_trt$`Chi Df`[2])
    df_int = c(df_int, anova_no_int$`Chi Df`[2])
    chisq_envt = c(chisq_envt, anova_no_envt$Chisq[2])
    chisq_trt = c(chisq_trt, anova_no_trt$Chisq[2])
    chisq_int = c(chisq_int, anova_no_int$Chisq[2])
    pval_envt = c(pval_envt, anova_no_envt$`Pr(>Chisq)`[2])
    pval_trt = c(pval_trt, anova_no_trt$`Pr(>Chisq)`[2])
    pval_int = c(pval_int, anova_no_int$`Pr(>Chisq)`[2])
    }
  }
  # gather results
  traits_envt_results = data.frame(trait=trait_labels, envt_var = envt_labels, df_envt, df_trt, df_int, chisq_envt, chisq_trt, chisq_int, pval_envt, pval_trt, pval_int)
  
  # correct for multiple comparisons using p.adjust or fdrtool
  # correction is for ALL variables simultaneously -- puts them all together, 
  #     applies fdr correction, then splits them back apart by explanatory variable
  if (!use.fdrtool) {
    n.models = length(pval_envt)
    pval_vec = c(pval_envt, pval_trt, pval_int)
    qval_vec = p.adjust(pval_vec, method="fdr")
    traits_envt_results$qval_envt = qval_vec[1:n.models]
    traits_envt_results$qval_trt = qval_vec[(n.models+1):(2*n.models)]
    traits_envt_results$qval_int = qval_vec[(2*n.models+1):(3*n.models)]
  }
  if (use.fdrtool) {
    n.models = length(pval_envt)
    pval_vec = c(pval_envt, pval_trt, pval_int)
    qval_vec = fdrtool(pval_vec, statistic="pvalue")$qval
    traits_envt_results$qval_envt = qval_vec[1:n.models]
    traits_envt_results$qval_trt = qval_vec[(n.models+1):(2*n.models)]
    traits_envt_results$qval_int = qval_vec[(2*n.models+1):(3*n.models)]
  }
  return(traits_envt_results)
}

#########################################
### Use this function to do the analyses 

############
### DAY 1 

popd1=d1$pop
MaternalLined1 = d1$MaternalLine
trtd1=d1$drought_water
PCd1 = d1$rachel_mds_2
# Specify which traits to analyze
traitsd1 = d1[,c("height_cm", "conductance", "stomatal_size", "stomatal_density", "SLA", "LWR", "area_cm2", "H4")]
# Specify which environmental variables to analyze 
envt_varsd1 = d1[,c("elevation", "MAP", "JanmaxT", "JulyminT", "summerRain", "rain_interannual_variability", "TempRange")]

# Mixed model with nested random efffects but no neutral structure PC
d1_traits_envt_results = fit_envt_trait_models(traits=traitsd1, envt_vars=envt_varsd1, trt=trtd1, pop=popd1, MaternalLine=MaternalLined1, PC=PCd1, use.PC=FALSE, use.fdrtool = FALSE)
head(d1_traits_envt_results)
filter(d1_traits_envt_results, qval_envt<0.05)

write.csv(d1_traits_envt_results, file="d1_traits_envt_results.csv")

# Mixed model with site random effects and neutral structure PC
d1_traits_envt_results_mds = fit_envt_trait_models(traits=traitsd1, envt_vars=envt_varsd1, trt=trtd1, pop=popd1, MaternalLine=MaternalLined1, PC=PCd1, use.PC=TRUE, use.fdrtool=FALSE)
head(d1_traits_envt_results_mds)
filter(d1_traits_envt_results_mds, qval_envt<0.05)

write.csv(d1_traits_envt_results_mds, file="d1_traits_envt_results_mds.csv")


############
### DAY 2 

popd2=d2$pop
MaternalLined2 = d2$MaternalLine
trtd2=d2$drought_water
PCd2 = d2$rachel_mds_2
# Specify which traits to analyze
traitsd2 = d2[,c("height_cm", "conductance", "stomatal_size", "stomatal_density", "SLA", "LWR", "area_cm2", "H4", "H4_diff", "root_length")]
# Specify which environmental variables to analyze 
envt_varsd2 = d2[,c("elevation", "MAP", "JanmaxT", "JulyminT", "summerRain", "rain_interannual_variability", "TempRange")]

# Mixed model with nested random efffects but no neutral structure PC
d2_traits_envt_results = fit_envt_trait_models(traits=traitsd2, envt_vars=envt_varsd2, trt=trtd2, pop=popd2, MaternalLine=MaternalLined2, PC=PCd2, use.PC=FALSE, use.fdrtool=FALSE)
head(d2_traits_envt_results)
filter(d2_traits_envt_results, qval_envt<0.05 | qval_int<0.05)

write.csv(d2_traits_envt_results, file="d2_traits_envt_results.csv")

# Mixed model with nested random efffects AND also neutral structure PC
d2_traits_envt_results_mds = fit_envt_trait_models(traits=traitsd2, envt_vars=envt_varsd2, trt=trtd2, pop=popd2, MaternalLine=MaternalLined2, PC=PCd2, use.PC=TRUE, use.fdrtool=FALSE)
head(d2_traits_envt_results_mds)
filter(d2_traits_envt_results_mds, qval_envt<0.05 | qval_int<0.05)

write.csv(d2_traits_envt_results_mds, file="d2_traits_envt_results_mds.csv")


### Visually compare results for the 2 models
options(scipen=2)
filter(d2_traits_envt_results[,c(1,2,12, 14)], qval_envt<0.05 | qval_int<0.05)
filter(d2_traits_envt_results_mds[,c(1,2,12, 14)], qval_envt<0.05 | qval_int<0.05)


# Make formatted table to match the numerical columns of Tables S6 and S7 
format_model_comp <- function(x) { # x is a table in format of traits_envt_results data frame outputted from the fit_envt_trait_models() function in this script
  envt_vars = as.vector(unique(x$envt_var), mode="character")
  traits = as.vector(unique(x$trait), mode="character")
  n.envt = length(envt_vars)
  n.traits = length(traits)
  # make the label columns of the table
  Variable = Fixed_effects = d.f = Chi.sq = p.val = q.val = {}
  Fixef_labels = insert(envt_vars, (1:n.envt)+1, values="Treatment")
  Fixef_labels = insert(Fixef_labels, ((1:n.envt)*2)+1, values="Interaction")
  for (i in 1:n.traits) {
    Variable = c(Variable, rep(traits[i], n.envt*3))
    Fixed_effects = c(Fixed_effects, Fixef_labels)
    for (j in 1:n.envt) {
      xtemp = filter(x, trait==traits[i] & envt_var==envt_vars[j])
      d.f = c(d.f, xtemp[1,"df_envt"], xtemp[1,"df_trt"], xtemp[1,"df_int"])
      Chi.sq = c(Chi.sq, xtemp[1,"chisq_envt"], xtemp[1,"chisq_trt"], xtemp[1,"chisq_int"])
      p.val = c(p.val, xtemp[1,"pval_envt"], xtemp[1,"pval_trt"], xtemp[1,"pval_int"])
      q.val = c(q.val, xtemp[1,"qval_envt"], xtemp[1,"qval_trt"], xtemp[1,"qval_int"])
    }
  }
  # Put table together
  table_output = data.frame(Variable, Fixed_effects, d.f, Chi.sq, p.val, q.val)
  table_output$Chi.sq = round(table_output$Chi.sq, 2)
  table_output$p.val = round(table_output$p.val, 3)
  table_output$q.val = round(table_output$q.val, 3)
  # clean up text, using the specific labels we want to put into the table
  # NOTE -- crappy code here, it works, but can only be used on the specific traits included in this script and analysis, not a general set of traits.
  table_output$Variable = unlist(sapply(as.character(table_output$Variable), switch, height_cm="Stem height", conductance="Stomatal conductance", stomatal_size="Stomatal pore size", stomatal_density="Stomatal density", SLA="SLA", LWR="Leaf length-width ratio", area_cm2 = "Leaf area", H4="Stem pigmentation", H4_diff="Stem pigment accumulation", root_length="Root length", soluble_shoot_all_g="Shoot soluble sugars", soluble_root_all_g="Root soluble sugars", starch_shoot_all_g="Shoot starch", starch_root_all_g="Root starch", growth_no_neg="Growth rate (height change day 6 to day 12)"))
  table_output$Fixed_effects = unlist(sapply(as.character(table_output$Fixed_effects), switch, Treatment="Treatment", Interaction="Interaction", elevation="Elevation", MAP="Mean annual precip.", JanmaxT="Jan. max. temp.", JulyminT="July min. temp", summerRain="Precip. seasonality", rain_interannual_variability="Interannual CV precip.", TempRange = "Annual temp. range"))
  return(table_output)
}

# Format output
table6_output_d1 = format_model_comp(d1_traits_envt_results)
table6_output_d2 = format_model_comp(d2_traits_envt_results)

write.csv(table6_output_d1, "table6_output_d1.csv")
write.csv(table6_output_d2, "table6_output_d2.csv")



############################################################################
# Fit models and make Table 7 -- relating performance traits to environment
############################################################################

popperf=d2$pop
MaternalLineperf = d2$MaternalLine
trtperf=d2$drought_water
PCperf = d2$rachel_mds_2
# Specify which traits to analyze
traitsperf = d2[,c("soluble_shoot_all_g", "soluble_root_all_g", "starch_shoot_all_g", "starch_root_all_g", "growth_no_neg", "growth_no_neg_normalized_initialHeight")]
# Specify which environmental variables to analyze 
envt_varsperf = d2[,c("elevation", "MAP", "JanmaxT", "JulyminT", "summerRain", "rain_interannual_variability", "TempRange")]

# Mixed model with nested random efffects but no neutral structure PC
perf_traits_envt_results = fit_envt_trait_models(traits=traitsperf, envt_vars=envt_varsperf, trt=trtperf, pop=popperf, MaternalLine=MaternalLineperf, PC=PCperf, use.PC=FALSE, use.fdrtool=FALSE)
head(perf_traits_envt_results)
filter(perf_traits_envt_results, qval_envt<0.05 | qval_int<0.05)

write.csv(perf_traits_envt_results, file="perf_traits_envt_results.csv")

# Mixed model with site random efffects and neutral structure 
perf_traits_envt_results_mds = fit_envt_trait_models(traits=traitsperf, envt_vars=envt_varsperf, trt=trtperf, pop=popperf, MaternalLine=MaternalLineperf, PC=PCperf, use.PC=TRUE, use.fdrtool=FALSE)
head(perf_traits_envt_results_mds)
filter(perf_traits_envt_results_mds, qval_envt<0.05 | qval_int<0.05)

write.csv(perf_traits_envt_results_mds, file="perf_traits_envt_results_mds.csv")

# format output
table7_output = format_model_comp(perf_traits_envt_results)
write.csv(table7_output, "table7_output.csv")



# Check trait correlations
alltraits = c("height_cm", "conductance", "stomatal_size", "stomatal_density", "SLA", "LWR", "area_cm2", "H4", "H4_diff", "root_length", "soluble_shoot_all_g", "soluble_root_all_g", "starch_shoot_all_g", "starch_root_all_g", "growth_no_neg")
cor(d2[d2$drought_water=="drought",alltraits], use="pairwise.complete.obs")
cor(d2[d2$drought_water=="water",alltraits], use="pairwise.complete.obs")



##########################################################
## Get raw correlations of traits and environment for the heatmap figure

# Trait and environment data are already contained in d1 and d2
# So we need to load the eigengene expression data and merge
eig.d1 <- read.csv("../Melis_analyses_makeover/voomd1_popxtreatGenesTraits.csv")
eig.d2 <- read.csv("../Melis_analyses_makeover/voomd2_popxtreatGenesTraits.csv")
head(eig.d1)

# Calculate trait-environment correlations 
# Choose variables -- make them in the same order as the existing Figure 4
envt_vars.cor <- c("elevation", "MAP", "JanmaxT", "JulyminT", "rain_interannual_variability", "summerRain", "TempRange")

alltraits.cor <- c("height_cm", "conductance", "stomatal_size", "stomatal_density", "SLA", "LWR", "area_cm2", "H4", "growth_no_neg", "root_length", "soluble_shoot_all_g", "soluble_root_all_g", "starch_shoot_all_g", "starch_root_all_g")
eig_names.cor <- c("dataME13", "dataME12","dataME11","dataME10","dataME9","dataME8","dataME7","dataME6","dataME5","dataME4","dataME3","dataME2","dataME1","dataME0")

# Both treatments together
# trait-environment correlations
traitcor.d1 <- cor(d1[,alltraits.cor], d1[,envt_vars.cor], use="pairwise.complete")
traitcor.d2 <- cor(d2[,alltraits.cor], d2[,envt_vars.cor], use="pairwise.complete")
# expression-environment correlations
eigcor.d1 <- cor(eig.d1[,eig_names.cor], eig.d1[,envt_vars.cor], use="pairwise.complete")
eigcor.d2 <- cor(eig.d2[,eig_names.cor], eig.d2[,envt_vars.cor], use="pairwise.complete")
# together
cor.d1.all <- rbind(eigcor.d1, traitcor.d1)
cor.d2.all <- rbind(eigcor.d2, traitcor.d2)

## Water treatment only
d1.wet <- filter(d1, drought_water=="water")
d2.wet <- filter(d2, drought_water=="water")
eig.d1.wet <- filter(eig.d1, treatment=="wet")
eig.d2.wet <- filter(eig.d2, treatment=="wet")
# trait-environment correlations
traitcor.d1.wet <- cor(d1.wet[,alltraits.cor], d1.wet[,envt_vars.cor], use="pairwise.complete")
traitcor.d2.wet <- cor(d2.wet[,alltraits.cor], d2.wet[,envt_vars.cor], use="pairwise.complete")
# expression-environment correlations
eigcor.d1.wet <- cor(eig.d1.wet[,eig_names.cor], eig.d1.wet[,envt_vars.cor], use="pairwise.complete")
eigcor.d2.wet <- cor(eig.d2.wet[,eig_names.cor], eig.d2.wet[,envt_vars.cor], use="pairwise.complete")
# together
cor.d1.wet <- rbind(eigcor.d1.wet, traitcor.d1.wet)
cor.d2.wet <- rbind(eigcor.d2.wet, traitcor.d2.wet)

## Drought treatment only
d1.dry <- filter(d1, drought_water=="drought")
d2.dry <- filter(d2, drought_water=="drought")
eig.d1.dry <- filter(eig.d1, treatment=="dry")
eig.d2.dry <- filter(eig.d2, treatment=="dry")
# trait-environment correlations
traitcor.d1.dry <- cor(d1.dry[,alltraits.cor], d1.dry[,envt_vars.cor], use="pairwise.complete")
traitcor.d2.dry <- cor(d2.dry[,alltraits.cor], d2.dry[,envt_vars.cor], use="pairwise.complete")
# expression-environment correlations
eigcor.d1.dry <- cor(eig.d1.dry[,eig_names.cor], eig.d1.dry[,envt_vars.cor], use="pairwise.complete")
eigcor.d2.dry <- cor(eig.d2.dry[,eig_names.cor], eig.d2.dry[,envt_vars.cor], use="pairwise.complete")
# together
cor.d1.dry <- rbind(eigcor.d1.dry, traitcor.d1.dry)
cor.d2.dry <- rbind(eigcor.d2.dry, traitcor.d2.dry)

rownames(cor.d1.all) = rownames(cor.d2.all) = rownames(cor.d1.wet) = rownames(cor.d2.wet) = rownames(cor.d1.dry) = rownames(cor.d2.dry) = c(eig_names.cor, alltraits.cor)
colnames(cor.d1.all) = colnames(cor.d2.all) = colnames(cor.d1.wet) = colnames(cor.d2.wet) = colnames(cor.d1.dry) = colnames(cor.d2.dry) = envt_vars.cor

write.csv(cor.d1.all, "cor_d1_all.csv")
write.csv(cor.d2.all, "cor_d2_all.csv")
write.csv(cor.d1.wet, "cor_d1_wet.csv")
write.csv(cor.d2.wet, "cor_d2_wet.csv")
write.csv(cor.d1.dry, "cor_d1_dry.csv")
write.csv(cor.d2.dry, "cor_d2_dry.csv")

image.plot(cor.d2.dry)

head(cor.d1.all)


