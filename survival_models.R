# Survival data analysis

setwd("~/Google Drive/ProteaPaper/survival_related")
library(survival)

# load survival data
d = read.table("survival_4Andrew.txt", header=T)
head(d)
d = d[d$pop != "A",]
d$pop = droplevels(d$pop)


# load trait data
plast = read.csv("final_RDPI_output_survival_averages.csv")
head(plast)
#traits = read.csv("../datasets/trait_data_MEs_PCaxes.csv")
traits = read.csv("../datasets/repens_drydown_alldata.csv")
# Focus only on the drought plants
traits = traits[traits$drought_water=="drought",]
head(traits)
dim(traits)
pc = read.csv("trait pca day 1 and day2.csv")


# Fix plant id codes in the pc data set 
pc$code = sapply(pc$code, gsub, pattern=".", replacement="_", fixed=TRUE)
# merge with survival data
pcsurv = merge(pc, d, by="code")
head(pcsurv)

# Test for associations between pcs and survival 
s = survreg(Surv(pcsurv$time1, pcsurv$time2, type="interval2")~day1_Prin1+frailty(pop.x), data=pcsurv, dist="logistic")
summary(s) # no association with PC1
s = survreg(Surv(pcsurv$time1, pcsurv$time2, type="interval2")~day1_Prin2+frailty(pop.x), data=pcsurv, dist="logistic")
summary(s) # marginal positive association p=0.062 -- higher PC2 score associated weakly with longer survival
plot(days_til_dead~day1_Prin2, pcsurv, col=1, pch=16)
abline(coef(lm(days_til_dead~day1_Prin2, pcsurv)), col="darkgray")
text(y=pcsurv$days_til_dead, x=pcsurv$day1_Prin2, labels=as.character(pcsurv$pop.x)) # weak 

s = survreg(Surv(pcsurv$time1, pcsurv$time2, type="interval2")~day1_Prin3+frailty(pop.x), data=pcsurv, dist="logistic")
summary(s) # no association with PC3



#################################
# survival analysis for trait means at the population level

# TO DO: First average up to maternal lines, then average from there up to the population mean. 


# read in data
#traits = read.csv("../datasets/trait_data_MEs_PCaxes.csv")
traits = read.csv("../datasets/repens_drydown_alldata.csv")
survdata = read.table("survival_4Andrew.txt", header=T)

# focus only on drought plants
traits = traits[traits$drought_water=="drought",]
survdata = survdata[survdata$drought_water=="drought",]

# Make sure there aren't any plants from population A
unique(traits$pop)
traits = traits[traits$pop != "A" & !is.na(traits$pop), ]
traits$pop = droplevels(traits$pop)
survdata = survdata[survdata$pop != "A",]
survdata$pop = droplevels(survdata$pop)

# Average trait data up to the population level 
#continuous_cols =  c(12:43)
continuous_cols = c(10:48, 56:71)
factor_cols = "pop"
traits_pop = aggregate(traits[,continuous_cols], by=list(traits$pop), FUN=mean, na.rm=T)
traits_pop = cbind(traits_pop, aggregate(traits[,factor_cols], by=list(traits$pop), FUN=f<-function(x){return(x[1])}))
head(traits_pop) 
names(traits_pop)[57] = "pop"

# Average survival data up to the population level
surv_pop = aggregate(survdata$days_til_dead, by=list(survdata$pop), FUN=mean, na.rm=T)
names(surv_pop) = c("pop", "days_til_dead")

# Merge
traits_surv = merge(traits_pop, surv_pop, by="pop")

# Test for relationships between trait means and survival 

# Carbohydrate traits
summary(lm(days_til_dead~SolubleSugarsRoots_conc, traits_surv))
summary(lm(days_til_dead~StarchRoots_conc, traits_surv))
summary(lm(days_til_dead~SolubleSugarsShoots_conc, traits_surv))
# p = 0.0618
plot(days_til_dead~SolubleSugarsShoots_conc, traits_surv, pch=16)
abline(coef(lm(days_til_dead~SolubleSugarsShoots_conc, traits_surv)), col="darkgray")
summary(lm(days_til_dead~StarchShoots_conc, traits_surv))

# Eigengene expression levels
summary(lm(days_til_dead~dataME1, traits_surv))
summary(lm(days_til_dead~dataME2, traits_surv))
summary(lm(days_til_dead~dataME3, traits_surv))
summary(lm(days_til_dead~dataME4, traits_surv))
summary(lm(days_til_dead~dataME5, traits_surv))
summary(lm(days_til_dead~dataME6, traits_surv))
summary(lm(days_til_dead~dataME7, traits_surv))
summary(lm(days_til_dead~dataME8, traits_surv))
summary(lm(days_til_dead~dataME9, traits_surv))
summary(lm(days_til_dead~dataME10, traits_surv))
summary(lm(days_til_dead~dataME11, traits_surv))
summary(lm(days_til_dead~dataME12, traits_surv))
summary(lm(days_til_dead~dataME13, traits_surv))
# nothing



# Traits available from the first measurement time (day 1)
traits1 = c("code", "pop", "meas_time", "conductance", "height_cm", "H4", "area_cm2", "LWR", "SLA", "root_length","stomatal_density", "stomatal_size", "SPI")
traits_day1 = traits[,traits1]
traits_day1 = traits[traits$meas_time=="1",]

testdata_day1 = merge(traits_day1, survdata, by="pop")

n_traits = length(traits1)
traitmodels_day1 = list()
for (i in 4:n_traits) {
  modform = as.formula(paste("days_til_dead~", traits1[i], " + (1|pop)", sep=""))
  traitmodels_day1[[i-3]] = lmer(formula=modform, data=testdata_day1)
}
traitmodels_day1
lapply(traitmodels_day1, summary)
# Nothing significant

# Traits available from the second measurement time (day 2)
traits2 = c("code", "pop", "meas_time", "conductance", "height_cm", "H4", "area_cm2", "LWR", "SLA", "root_length","stomatal_density", "stomatal_size", "SPI")
traits_day2 = traits[,traits2]
traits_day2 = traits[traits$meas_time=="2",]

testdata_day2 = merge(traits_day2, survdata, by="pop")

n_traits = length(traits2)
traitmodels_day2 = list()
for (i in 4:n_traits) {
  modform = as.formula(paste("days_til_dead~", traits1[i], " + (1|pop)", sep=""))
  traitmodels_day2[[i-3]] = lmer(formula=modform, data=testdata_day2)
}
lapply(traitmodels_day2, summary)

# Nothing significant

