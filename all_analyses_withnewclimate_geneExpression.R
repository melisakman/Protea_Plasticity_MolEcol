################
##DGE Analyses##
################
library(edgeR)
reads = read.csv("ProExp_counts.csv", header=TRUE) # first csv should have the contig names in the first column,
head(reads)
dim(reads)

#define count data
counts = reads[,c(2:121)]
# In eXpression website they recommend using rounded eff_counts, so we will first round up the counts
counts = round(counts)
head(counts)
barplot(colSums(counts)*1e-6, names =1:120, ylab="Library size (millions)")
dim(counts)
contigs=reads[,1]
counts_d1= reads[,c(2:65)]
counts_d2= reads[,c(66:121)]

# Indicates which population each individual is from
pops = as.vector(sapply(names(counts), substr, start=15, stop=15)) # population ids in order of the columns of the reads matrix
timepoint = as.vector(sapply(names(counts), substr, start=11, stop=11)) # pull out timepoint
treatment = as.vector(sapply(names(counts), substr, start=13, stop=13)) # pull out timepoint
batch= as.vector(sapply(names(counts), substr, start=8, stop=8))
group= paste(timepoint, pops,treatment, sep="_")
pops_d1 = as.vector(sapply(names(counts_d1), substr, start=15, stop=15)) # population ids in order of the columns of the reads matrix
treatment_d1 = as.vector(sapply(names(counts_d1), substr, start=13, stop=13)) # pull out timepoint
batch_d1= as.vector(sapply(names(counts_d1), substr, start=8, stop=8))
group_d1= paste(pops_d1,treatment_d1, sep='_')

pops_d2 = as.vector(sapply(names(counts_d2), substr, start=15, stop=15)) # population ids in order of the columns of the reads matrix
treatment_d2 = as.vector(sapply(names(counts_d2), substr, start=13, stop=13)) # pull out timepoint
batch_d2= as.vector(sapply(names(counts_d2), substr, start=8, stop=8))
group_d2= paste(pops_d2,treatment_d2, sep='_')
# create dge object
dge = DGEList(counts, genes=contigs, group= group)
dge_d1 = DGEList(counts_d1, genes=contigs, group=group_d1)
dge_d2 = DGEList(counts_d2, genes=contigs, group=group_d2)

# remove contigs with low counts (Keep genes with least 1 count-per-million reads (cpm) in at least half of samples samples)
# it makes more sensefor this to be done seperately for afternoons and mornings 
expr = rowSums(cpm(dge)>1)>=96
expr_d1 = rowSums(cpm(dge_d1)>1)>=51
expr_d2 = rowSums(cpm(dge_d2)>1)>=45

dge = dge[expr, keep.lib.sizes=FALSE]
dge_d1 = dge_d1[expr_d1, keep.lib.sizes=FALSE]
dge_d2 = dge_d2[expr_d2, keep.lib.sizes=FALSE]

dim(dge)
dim(dge_d1)
dim(dge_d2)

#normalization for lib size before voom
reads.norm = calcNormFactors(dge)
reads.norm_d1 = calcNormFactors(dge_d1)
reads.norm_d2 = calcNormFactors(dge_d2)

##########################
## Pairwise comparisons ##
##########################

design_d1= model.matrix(~batch_d1+group_d1, data= dge_d1$samples)
design_d1
design_d2= model.matrix(~batch_d2+group_d2, data= dge_d2$samples)
design_d2

dge_d1= estimateDisp(dge_d1, design_d1)  #this might take a while
dge_d2= estimateDisp(dge_d2, design_d2)
dge= estimateTagwiseDisp(dge, design)

# Let's plot dispersion or BCV
plotBCV(dge_d1)
plotBCV(dge_d2)

# Now we fit our negative binomial model 

fit=glmFit(dge_d1,design_d1)
fit=glmFit(dge_d2,design_d2)

# Perform likelihood ratio tests d1
lrtA=glmLRT(fit, coef=3) # pop A
lrtB=glmLRT(fit, contrast=c(0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0)) #popB
lrtC=glmLRT(fit, contrast=c(0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0)) #popC
lrtF=glmLRT(fit, contrast=c(0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0)) #popF
lrtG=glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0)) #popG
lrtR=glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0)) #popR
lrtS=glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0)) #popS
lrtV=glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1))#popV

# Perform likelihood ratio tests d2
lrtB=glmLRT(fit, coef=3) # pop B
lrtC=glmLRT(fit, contrast=c(0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0)) #popC
lrtF=glmLRT(fit, contrast=c(0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0)) #popF
lrtG=glmLRT(fit, contrast=c(0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0)) #popG
lrtR=glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0)) #popR
lrtS=glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0)) #popS
lrtV=glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1)) #popV


topTags(lrtV)
# and now get the top genes
a=topTags(lrtV, n=7000)
write.table(a, file="d2_DE_popV.csv", sep=",")

##################################################################
## Differentially regulated genes among populations (popxtreat) ##
##################################################################

design_d1= model.matrix(~batch_d1+pops_d1*treatment_d1, data= dge_d1$samples)
colnames(design_d1)
design_d2= model.matrix(~batch_d2+pops_d2*treatment_d2, data= dge_d2$samples)
design_d2

dge_d1= estimateDisp(dge_d1, design_d1)  #this might take a while
dge_d2= estimateDisp(dge_d2, design_d2)

# Let's plot dispersion or BCV
plotBCV(dge_d1)
plotBCV(dge_d2)

# Now we fit our negative binomial model 

fit_d1=glmFit(dge_d1,design_d1)
fit_d2=glmFit(dge_d2,design_d2)

# Perform likelihood ratio tests
lrt_d1=glmLRT(fit_d1, coef=11:17) 
lrt_d2=glmLRT(fit_d2, coef=10:15) 


topTags(lrt_d1)
topTags(lrt_d2)
# and now get the top genes
a=topTags(lrt_d1, n=7000)
b=topTags(lrt_d2, n=7000)

write.table(a, file="DE_d1.csv", sep=",")
write.table(b, file="DE_d2.csv", sep=",")

#######################################################################
## Extracting gene expression patterns of popxtreatment effect genes ##
#######################################################################
## In order to plot these I am using voom-transformed values

library(ggplot2)
library(reshape2)
d1_pxt = read.csv("voomd1_popxtreatGenes.csv", header=TRUE)
d2_pxt = read.csv("voomd2_popxtreatGenes.csv", header=TRUE)

ggplot(data=d1_pxt, aes(x=pop, y=contig39647, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d1_pxt, aes(x=pop, y=contig58964, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d1_pxt, aes(x=pop, y=contig56108, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d1_pxt, aes(x=pop, y=contig77386, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d1_pxt, aes(x=pop, y=contig58961, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d1_pxt, aes(x=pop, y=contig93102, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d1_pxt, aes(x=pop, y=contig5728, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d1_pxt, aes(x=pop, y=contig24218, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d1_pxt, aes(x=pop, y=contig44885, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d1_pxt, aes(x=pop, y=contig58962, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d1_pxt, aes(x=pop, y=contig30375, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d1_pxt, aes(x=pop, y=contig11561, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d1_pxt, aes(x=pop, y=contig116675, fill=factor(treatment))) + geom_boxplot()

ggplot(data=d2_pxt, aes(x=pop, y=contig60140, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d2_pxt, aes(x=pop, y=contig89198, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d2_pxt, aes(x=pop, y=contig22152, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d2_pxt, aes(x=pop, y=contig9620, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d2_pxt, aes(x=pop, y=contig94608, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d2_pxt, aes(x=pop, y=contig7425, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d2_pxt, aes(x=pop, y=contig19237, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d2_pxt, aes(x=pop, y=contig9970, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d2_pxt, aes(x=pop, y=contig18074, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d2_pxt, aes(x=pop, y=contig9868, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d2_pxt, aes(x=pop, y=contig42012, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d2_pxt, aes(x=pop, y=contig97840, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d2_pxt, aes(x=pop, y=contig19472, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d2_pxt, aes(x=pop, y=contig29616, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d2_pxt, aes(x=pop, y=contig32983, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d2_pxt, aes(x=pop, y=contig49, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d2_pxt, aes(x=pop, y=contig12722, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d2_pxt, aes(x=pop, y=contig90959, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d2_pxt, aes(x=pop, y=contig25073, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d2_pxt, aes(x=pop, y=contig115461, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d2_pxt, aes(x=pop, y=contig3503, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d2_pxt, aes(x=pop, y=contig5102, fill=factor(treatment))) + geom_boxplot()
ggplot(data=d2_pxt, aes(x=pop, y=contig31354, fill=factor(treatment))) + geom_boxplot()


##############################
## with treatments TOGETHER ##
#######################################################################
## Correlations between popxtreatment effect genes, envir and traits ##
#######################################################################

# In order to look at correlations I am going to use voom-transformed genes. If results are not good then I will use actual normalized gene counts
install.packages("lme4")
install.packages("lmerTest")
library(lme4)
library(lmerTest)

d1T = read.csv("voomd1_popxtreatGenesTraits.csv", header=TRUE, na.strings = "NA")
d2T = read.csv("voomd2_popxtreatGenesTraits.csv", header=TRUE, na.strings = "NA")

###############################################
##function to be used for looping for traits###
###############################################

fit.lme= function(dtemp){
  fit=lmer(x ~ y + trt + (1|batch) + (1|pop), REML=F, data = dtemp)
  print(summary(fit))
  fit.null=lmer(x ~ trt + (1|batch) + (1|pop), REML=F, data = dtemp)
  anova(fit.null, fit)
}

fit.lme.int= function(dtemp){
  fit=lmer(x ~ y + trt + y*trt + (1|batch) + (1|pop), REML=F, data = dtemp)
  print(summary(fit))
  fit.null=lmer(x~ y + trt + (1|batch) + (1|pop), REML=F, data = dtemp)
  anova(fit.null, fit)
}

fit.lme.trt= function(dtemp){
  fit=lmer(x ~ y + trt + (1|batch) + (1|pop), REML=F, data = dtemp)
  print(summary(fit))
  fit.null=lmer(x ~ y + (1|batch) + (1|pop), REML=F, data = dtemp)
  anova(fit.null, fit)
}

###########################################################################################
## new functions including maternal effects Rachel MDS to be used for looping for traits###
###########################################################################################

fit.lme= function(dtemp){
  fit=lmer(x ~ y + trt + MaternalLineMDS + (1|batch), REML=F, data = dtemp)
  print(summary(fit))
  fit.null=lmer(x ~ trt +  MaternalLineMDS + (1|batch), REML=F, data = dtemp)
  anova(fit.null, fit)
}

fit.lme.int= function(dtemp){
  fit=lmer(x ~ y + trt + y*trt + MaternalLineMDS + (1|batch), REML=F, data = dtemp)
  print(summary(fit))
  fit.null=lmer(x~ y + trt + MaternalLineMDS + (1|batch), REML=F, data = dtemp)
  anova(fit.null, fit)
}

fit.lme.trt= function(dtemp){
  fit=lmer(x ~ y + trt + MaternalLineMDS + (1|batch), REML=F, data = dtemp)
  print(summary(fit))
  fit.null=lmer(x ~ y + MaternalLineMDS + (1|batch), REML=F, data = dtemp)
  anova(fit.null, fit)
}

################################################
# comparing d1 popxtreatment contigs to traits #
################################################

names(d1T)
pop=d1T$pop
trt=d1T$treatment
MaternalLineMDS = d1T$rachel_mds_2
batch=d1T$seqBatch
genesd1T=d1T[,c(6:32)]
n.genesd1T=dim(genesd1T)[2]
traitsd1T=d1T[,c(33:39)]
n.traitsd1T=dim(traitsd1T)[2]

#fit=lmer(genesd1T[,1] ~ traitsd1T[,2] + trt + y*trt + (1|batch) + (1|pop), REML=F, data = dtemp)
#summary(fit)
#coef(fit)

sink('d1genestraits.txt')
lme.traits={}
for (i in 1:n.genesd1T) {
  x = genesd1T[,i]
  for (j in 1:n.traitsd1T) {
    y = traitsd1T[,j]
    dtemp = data.frame(x=x, y=y, trt = d1T$treatment, batch = d1T$seqBatch, pop = d1T$pop)
    # remove any NAs from y, and also remove corresponding rows of other variables
    missings = which(is.na(y))
    if (length(missings)>0) dtemp = dtemp[-missings,]
    lme.traits=rbind(lme.traits, fit.lme(dtemp))
  }
}
sink()
write.csv(lme.traits, file="d1_traits.csv")

sink('d1genestraits_interaction.txt')
lme.traits.int={}
for (i in 1:n.genesd1T) {
  x = genesd1T[,i]
  for (j in 1:n.traitsd1T) {
    y = traitsd1T[,j]
    dtemp = data.frame(x=x, y=y, trt = d1T$treatment, batch = d1T$seqBatch, pop = d1T$pop)
    # remove any NAs from y, and also remove corresponding rows of other variables
    missings = which(is.na(y))
    if (length(missings)>0) dtemp = dtemp[-missings,]
    lme.traits.int=rbind(lme.traits.int, fit.lme.int(dtemp))
  }
}
sink()
write.csv(lme.traits.int, file="d1_traits_interaction.csv")

sink('d1genestraits_treatment.txt')
lme.traits.trt={}
for (i in 1:n.genesd1T) {
  x = genesd1T[,i]
  for (j in 1:n.traitsd1T) {
    y = traitsd1T[,j]
    dtemp = data.frame(x=x, y=y, trt = d1T$treatment, batch = d1T$seqBatch, pop = d1T$pop)
    # remove any NAs from y, and also remove corresponding rows of other variables
    missings = which(is.na(y))
    if (length(missings)>0) dtemp = dtemp[-missings,]
    lme.traits.trt=rbind(lme.traits.trt, fit.lme.trt(dtemp))
  }
}
sink()
write.csv(lme.traits.trt, file="d1_traits_trt.csv")


################################################################
# comparing d2 popxtreatment contigs and eigen genes to traits #
################################################################

names(d2T)
pop=d2T$pop
trt=d2T$treatment
batch=d2T$seqBatch
genesd2T=d2T[,c(6:42)]
n.genesd2T=dim(genesd2T)[2]
traitsd2T=d2T[,c(43:61)]
n.traitsd2T=dim(traitsd2T)[2]

sink('d2genestraits.txt')
lme.traits={}
for (i in 1:n.genesd2T) {
  x = genesd2T[,i]
  for (j in 1:n.traitsd2T) {
    y = traitsd2T[,j]
    dtemp = data.frame(x=x, y=y, trt = d2T$treatment, batch = d2T$seqBatch, pop = d2T$pop)
    # remove any NAs from y, and also remove corresponding rows of other variables
    missings = which(is.na(y))
    if (length(missings)>0) dtemp = dtemp[-missings,]
    lme.traits=rbind(lme.traits, fit.lme(dtemp))
  }
}
sink()
write.csv(lme.traits, file="d2_traits.csv")

sink('d2genestraits_interaction.txt')
lme.traits.int={}
for (i in 1:n.genesd2T) {
  x = genesd2T[,i]
  for (j in 1:n.traitsd2T) {
    y = traitsd2T[,j]
    dtemp = data.frame(x=x, y=y, trt = d2T$treatment, batch = d2T$seqBatch, pop = d2T$pop)
    # remove any NAs from y, and also remove corresponding rows of other variables
    missings = which(is.na(y))
    if (length(missings)>0) dtemp = dtemp[-missings,]
    lme.traits.int=rbind(lme.traits.int, fit.lme.int(dtemp))
  }
}
sink()
write.csv(lme.traits.int, file="d2_traits_interaction.csv")

sink('d2genestraits_trt.txt')
lme.traits.trt={}
for (i in 1:n.genesd2T) {
  x = genesd2T[,i]
  for (j in 1:n.traitsd2T) {
    y = traitsd2T[,j]
    dtemp = data.frame(x=x, y=y, trt = d2T$treatment, batch = d2T$seqBatch, pop = d2T$pop)
    # remove any NAs from y, and also remove corresponding rows of other variables
    missings = which(is.na(y))
    if (length(missings)>0) dtemp = dtemp[-missings,]
    lme.traits.trt=rbind(lme.traits.trt, fit.lme.trt(dtemp))
  }
}
sink()
write.csv(lme.traits.trt, file="d2_traits_trt.csv")


#########################################################################
##function to be used for looping for climate and different day traits ##
#########################################################################

fit.lme= function(dtemp){
  fit=lmer(x ~ y + trt  + (1|pop) + (1|batch), REML=F, data = dtemp)
  print(summary(fit))
  fit.null=lmer(x~ trt + (1|pop) + (1|batch), REML=F, data = dtemp)
  anova(fit.null, fit)
}

fit.lme.int= function(dtemp){
  fit=lmer(x ~ y + trt + trt*y + (1|pop) + (1|batch), REML=F, data = dtemp)
  print(summary(fit))
  fit.null=lmer(x ~ y + trt + (1|pop) + (1|batch), REML=F, data = dtemp)
  anova(fit.null, fit)
}

fit.lme.trt= function(dtemp){
  fit=lmer(x ~ y + trt + (1|pop) + (1|batch), REML=F, data = dtemp)
  print(summary(fit))
  fit.null=lmer(x ~ y + (1|pop) + (1|batch), REML=F, data = dtemp)
  anova(fit.null, fit)
}

###############################################################################################
##function to be used for looping for climate and different day traits with Rachels MDS axes ##
###############################################################################################

fit.lme= function(dtemp){
  fit=lmer(x ~ y + trt  + MaternalLineMDS + (1|batch), REML=F, data = dtemp)
  print(summary(fit))
  fit.null=lmer(x~ trt + MaternalLineMDS + (1|batch), REML=F, data = dtemp)
  anova(fit.null, fit)
}

fit.lme.int= function(dtemp){
  fit=lmer(x ~ y + trt + trt*y + MaternalLineMDS + (1|batch), REML=F, data = dtemp)
  print(summary(fit))
  fit.null=lmer(x ~ y + trt + MaternalLineMDS + (1|batch), REML=F, data = dtemp)
  anova(fit.null, fit)
}

fit.lme.trt= function(dtemp){
  fit=lmer(x ~ y + trt + MaternalLineMDS + (1|batch), REML=F, data = dtemp)
  print(summary(fit))
  fit.null=lmer(x ~ y + MaternalLineMDS + (1|batch), REML=F, data = dtemp)
  anova(fit.null, fit)
}

#####################################################################################################
##function to be used for looping for climate and different day traits with maternal nested in pop ##
#####################################################################################################

fit.lme= function(dtemp){
  fit=lmer(x ~ y + trt  + (1|pop/MaternalLine) + (1|batch), REML=F, data = dtemp)
  print(summary(fit))
  fit.null=lmer(x~ trt + (1|pop/MaternalLine) + (1|batch), REML=F, data = dtemp)
  anova(fit.null, fit)
}

fit.lme.int= function(dtemp){
  fit=lmer(x ~ y + trt + trt*y + (1|pop/MaternalLine) + (1|batch), REML=F, data = dtemp)
  print(summary(fit))
  fit.null=lmer(x ~ y + trt + (1|pop/MaternalLine) + (1|batch), REML=F, data = dtemp)
  anova(fit.null, fit)
}

fit.lme.trt= function(dtemp){
  fit=lmer(x ~ y + trt + (1|pop/MaternalLine) + (1|batch), REML=F, data = dtemp)
  print(summary(fit))
  fit.null=lmer(x ~ y + (1|pop/MaternalLine) + (1|batch), REML=F, data = dtemp)
  anova(fit.null, fit)
}


#########################################################################################################################################
##function to be used for looping for climate and different day traits with maternal nested in pop and MDS pop average as fixed factor ##
#########################################################################################################################################

fit.lme= function(dtemp){
  fit=lmer(x ~ y + trt  + MaternalLineMDS + (1|pop/MaternalLine) + (1|batch), REML=F, data = dtemp)
  print(summary(fit))
  fit.null=lmer(x~ trt + MaternalLineMDS + (1|pop/MaternalLine) + (1|batch), REML=F, data = dtemp)
  anova(fit.null, fit)
}

fit.lme.int= function(dtemp){
  fit=lmer(x ~ y + trt + trt*y + popMDS + (1|pop/MaternalLine) + (1|batch), REML=F, data = dtemp)
  print(summary(fit))
  fit.null=lmer(x ~ y + trt + popMDS + (1|pop/MaternalLine) + (1|batch), REML=F, data = dtemp)
  anova(fit.null, fit)
}

fit.lme.trt= function(dtemp){
  fit=lmer(x ~ y + trt + popMDS + (1|pop/MaternalLine) + (1|batch), REML=F, data = dtemp)
  print(summary(fit))
  fit.null=lmer(x ~ y + popMDS + (1|pop/MaternalLine) + (1|batch), REML=F, data = dtemp)
  anova(fit.null, fit)
}
###############################################################################
# comparing d1 genes and eigen genes to climate and traits from the other day #
###############################################################################
d1TE = read.csv("voomd1_popxtreatGenesTraits.csv", header=TRUE, na.strings = "NA")
d2TE = read.csv("voomd2_popxtreatGenesTraits.csv", header=TRUE, na.strings = "NA")

names(d1TE)
pop=d1TE$pop
trt=d1TE$treatment
batch=d1TE$seqBatch
MaternalLineMDS = d1TE$rachel_mds_2
popMDS = d1TE$popMDS
MaternalLine = d1TE$MaternalLine
genesd1TE=d1TE[,c(6:32)]
n.genesd1TE=dim(genesd1TE)[2]
climated1TE=d1TE[,c(40:46)]
n.climated1TE=dim(climated1TE)[2]

sink('d1genesVSclimate_maternalMDS_maternalLine_nestedpop.txt')
lme.env={}
for (i in 1:n.genesd1TE) {
  x = genesd1TE[,i]
  for (j in 1:n.climated1TE) {
    y = climated1TE[,j]
    dtemp = data.frame(x=x, y=y, trt = d1TE$treatment, pop = d1TE$pop, popMDS = d1TE$popMDS, MaternalLine= d1TE$MaternalLine, batch=d1TE$seqBatch)
    # remove any NAs from y, and also remove corresponding rows of other variables
    missings = which(is.na(y))
    if (length(missings)>0) dtemp = dtemp[-missings,]
    lme.env=rbind(lme.env, fit.lme(dtemp))
  }
}
sink()
write.csv(lme.env, file="d1_gene_climate_maternalMDS_maternalLine_nestedpop.csv")


sink('d1genesVSclimate_interaction_popMDS_maternal_nested_pop.txt')
lme.env.int={}
for (i in 1:n.genesd1TE) {
  x = genesd1TE[,i]
  for (j in 1:n.climated1TE) {
    y = climated1TE[,j]
    dtemp = data.frame(x=x, y=y, trt = d1TE$treatment, pop = d1TE$pop, batch = d1TE$seqBatch)
    # remove any NAs from y, and also remove corresponding rows of other variables
    missings = which(is.na(y))
    if (length(missings)>0) dtemp = dtemp[-missings,]
    lme.env.int=rbind(lme.env.int, fit.lme.int(dtemp))
  }
}
sink()
write.csv(lme.env.int, file="d1_gene_climate_interaction_popMDS_maternal_nested_pop.csv")

sink('d1genesVSclimate_trt_popMDS_maternal_nested_pop.txt')
lme.env.trt={}
for (i in 1:n.genesd1TE) {
  x = genesd1TE[,i]
  for (j in 1:n.climated1TE) {
    y = climated1TE[,j]
    dtemp = data.frame(x=x, y=y, trt = d1TE$treatment, pop = d1TE$pop, batch = d1TE$seqBatch)
    # remove any NAs from y, and also remove corresponding rows of other variables
    missings = which(is.na(y))
    if (length(missings)>0) dtemp = dtemp[-missings,]
    lme.env.trt=rbind(lme.env.trt, fit.lme.trt(dtemp))
  }
}
sink()
write.csv(lme.env.trt, file="d1_gene_climate_trt_popMDS_maternal_nested_pop.csv")

###############################################################################
# comparing d2 genes and eigen genes to climate and traits from the other day #
###############################################################################

names(d2TE)
pop=d2TE$pop
trt=d2TE$treatment
batch=d2TE$seqBatch
genesd2TE=d2TE[,c(6:42)]
n.genesd2TE=dim(genesd2TE)[2]
climated2TE=d2TE[,c(62:68)]
n.climated2TE=dim(climated2TE)[2]

sink('d2genesVSclimate.txt')
lme.env={}
for (i in 1:n.genesd2TE) {
  x = genesd2TE[,i]
  for (j in 1:n.climated2TE) {
    y = climated2TE[,j]
    dtemp = data.frame(x=x, y=y, trt = d2TE$treatment, pop = d2TE$pop, batch = d2TE$seqBatch)
    # remove any NAs from y, and also remove corresponding rows of other variables
    missings = which(is.na(y))
    if (length(missings)>0) dtemp = dtemp[-missings,]
    lme.env=rbind(lme.env, fit.lme(dtemp))
  }
}
sink()
write.csv(lme.env, file="d2_gene_climate.csv")

sink('d2genesVSclimate_int.txt')
lme.env.int={}
for (i in 1:n.genesd2TE) {
  x = genesd2TE[,i]
  for (j in 1:n.climated2TE) {
    y = climated2TE[,j]
    dtemp = data.frame(x=x, y=y, trt = d2TE$treatment, pop = d2TE$pop, batch = d2TE$seqBatch)
    # remove any NAs from y, and also remove corresponding rows of other variables
    missings = which(is.na(y))
    if (length(missings)>0) dtemp = dtemp[-missings,]
    lme.env.int=rbind(lme.env.int, fit.lme.int(dtemp))
  }
}
sink()
write.csv(lme.env.int, file="d2_gene_climate_interaction.csv")

sink('d2genesVSclimate_trt.txt')
lme.env.trt={}
for (i in 1:n.genesd2TE) {
  x = genesd2TE[,i]
  for (j in 1:n.climated2TE) {
    y = climated2TE[,j]
    dtemp = data.frame(x=x, y=y, trt = d2TE$treatment, pop = d2TE$pop, batch = d2TE$seqBatch)
    # remove any NAs from y, and also remove corresponding rows of other variables
    missings = which(is.na(y))
    if (length(missings)>0) dtemp = dtemp[-missings,]
    lme.env.trt=rbind(lme.env.trt, fit.lme.trt(dtemp))
  }
}
sink()
write.csv(lme.env.trt, file="d2_gene_climate_trt.csv")

##################
##fdr correction##
##################
install.packages("fdrtool")
library(fdrtool)
d1_p = read.csv("d2_contig_climate_interaction_trt.csv")
p_values = d1_p[,12]
q_values = p.adjust(p_values, method="fdr")
write.csv(q_values, file="d2_contig_climate_interaction_trt_qvalues.csv")

##############################
## with treatments SEPARATE ##
#######################################################################
## Correlations between popxtreatment effect genes, envir and traits ##
#######################################################################

# In order to look at correlations I am going to use voom-transformed genes. If results are not good then I will use actual normalized gene counts
install.packages("lme4")
install.packages("lmerTest")
library(lme4)
library(lmerTest)

d1TD = read.csv("voomd1_popxtreatGenesTraits_dry.csv", header=TRUE, na.strings = "NA")
d1TW = read.csv("voomd1_popxtreatGenesTraits_wet.csv", header=TRUE, na.strings = "NA")
d2TD = read.csv("voomd2_popxtreatGenesTraits_dry.csv", header=TRUE, na.strings = "NA")
d2TW = read.csv("voomd2_popxtreatGenesTraits_wet.csv", header=TRUE, na.strings = "NA")
d1ED = read.csv("voomd1_popxtreatGenesEnv_ave_dry.csv", header=TRUE, na.strings = "NA")
d1EW = read.csv("voomd1_popxtreatGenesEnv_ave_wet.csv", header=TRUE, na.strings = "NA")
d2ED = read.csv("voomd2_popxtreatGenesEnv_ave_dry.csv", header=TRUE, na.strings = "NA")
d2EW = read.csv("voomd2_popxtreatGenesEnv_ave_wet.csv", header=TRUE, na.strings = "NA")

###############################################
##function to be used for looping for traits###
###############################################

fit.lme= function(dtemp){
  fit=lmer(x ~ y + (1|batch) + (1|pop), REML=F, data = dtemp)
  print(summary(fit))
  fit.null=lmer(x~ (1|batch) + (1|pop), REML=F, data = dtemp)
  anova(fit.null, fit)
}

################################################
# comparing d1 popxtreatment contigs to traits #
################################################

names(d1TD)
pop=d1TD$pop
batch=d1TD$seqBatch
genesd1TD=d1TD[,c(6:32)]
n.genesd1TD=dim(genesd1TD)[2]
traitsd1TD=d1TD[,c(33:39)]
n.traitsd1TD=dim(traitsd1TD)[2]

fit=lmer(genesd1TD[,1] ~ traitsd1TD[,2] + (1|batch) + (1|pop), REML=F)
summary(fit)
#coef(fit)

sink('d1contigstraits_dry.txt')
lme.traits={}
for (i in 1:n.genesd1TD) {
  x = genesd1TD[,i]
  for (j in 1:n.traitsd1TD) {
    y = traitsd1TD[,j]
    dtemp = data.frame(x=x, y=y, batch = d1TD$seqBatch, pop = d1TD$pop)
    # remove any NAs from y, and also remove corresponding rows of other variables
    missings = which(is.na(y))
    if (length(missings)>0) dtemp = dtemp[-missings,]
    lme.traits=rbind(lme.traits, fit.lme(dtemp))
  }
}
sink()
write.csv(lme.traits, file="d1_traits_final_DRY.csv")

names(d1TW)
pop=d1TW$pop
batch=d1TW$seqBatch
genesd1TW=d1TW[,c(6:32)]
n.genesd1TW=dim(genesd1TW)[2]
traitsd1TW=d1TW[,c(33:39)]
n.traitsd1TW=dim(traitsd1TW)[2]

sink('d1contigstraits_wet.txt')
lme.traits={}
for (i in 1:n.genesd1TW) {
  x = genesd1TW[,i]
  for (j in 1:n.traitsd1TW) {
    y = traitsd1TW[,j]
    dtemp = data.frame(x=x, y=y, batch = d1TW$seqBatch, pop = d1TW$pop)
    # remove any NAs from y, and also remove corresponding rows of other variables
    missings = which(is.na(y))
    if (length(missings)>0) dtemp = dtemp[-missings,]
    lme.traits=rbind(lme.traits, fit.lme(dtemp))
  }
}
sink()
write.csv(lme.traits, file="d1_traits_final_WET.csv")

################################################################
# comparing d2 popxtreatment contigs and eigen genes to traits #
################################################################

names(d2TD)
pop=d2TD$pop
batch=d2TD$seqBatch
genesd2TD=d2TD[,c(6:42)]
n.genesd2TD=dim(genesd2TD)[2]
traitsd2TD=d2TD[,c(43:61)]
n.traitsd2TD=dim(traitsd2TD)[2]

sink('d2contigstraits_dry.txt')
lme.traits={}
for (i in 1:n.genesd2TD) {
  x = genesd2TD[,i]
  for (j in 1:n.traitsd2TD) {
    y = traitsd2TD[,j]
    dtemp = data.frame(x=x, y=y, trt = d2TD$treatment, batch = d2TD$seqBatch, pop = d2TD$pop)
    # remove any NAs from y, and also remove corresponding rows of other variables
    missings = which(is.na(y))
    if (length(missings)>0) dtemp = dtemp[-missings,]
    lme.traits=rbind(lme.traits, fit.lme(dtemp))
  }
}
sink()
write.csv(lme.traits, file="d2_traits_final_DRY.csv")

names(d2TW)
pop=d2TW$pop
batch=d2TW$seqBatch
genesd2TW=d2TW[,c(6:42)]
n.genesd2TW=dim(genesd2TW)[2]
traitsd2TW=d2TW[,c(43:61)]
n.traitsd2TW=dim(traitsd2TW)[2]

sink('d2contigstraits_wet.txt')
lme.traits={}
for (i in 1:n.genesd2TW) {
  x = genesd2TW[,i]
  for (j in 1:n.traitsd2TW) {
    y = traitsd2TW[,j]
    dtemp = data.frame(x=x, y=y, trt = d2TW$treatment, batch = d2TW$seqBatch, pop = d2TW$pop)
    # remove any NAs from y, and also remove corresponding rows of other variables
    missings = which(is.na(y))
    if (length(missings)>0) dtemp = dtemp[-missings,]
    lme.traits=rbind(lme.traits, fit.lme(dtemp))
  }
}
sink()
write.csv(lme.traits, file="d2_traits_final_WET.csv")

#########################################################################
##function to be used for looping for climate and different day traits ##
#########################################################################

fit.lm= function(x, y){
  fit=lm(x~y)		
  a=summary(fit)
  print(coef(a)[,4])
}

###############################################################################
# comparing d1 genes and eigen genes to climate and traits from the other day #
###############################################################################

names(d1ED)
pop=d1ED$pop
genesd1ED=d1ED[,c(5:31)]
n.genesd1ED=dim(genesd1ED)[2]
climated1Ed2TD=d1ED[,c(32:57)]
n.climated1Ed2TD=dim(climated1Ed2TD)[2]

fit.lm= function(x, y){
  fit=lm(x~y)		
  a=summary(fit)
  print(coef(a)[,4])
}

lme.climate={}
for (i in 1:n.genesd1ED) {
  x = genesd1ED[,i]
  for (j in 1:n.climated1Ed2TD) {
    y = climated1Ed2TD[,j]
    lme.climate=rbind(lme.climate, fit.lm(genesd1ED[,i],climated1Ed2TD[,j]))
  }
}
write.csv(lme.climate, file="d1_climate_final_dry.csv")

names(d1EW)
pop=d1EW$pop
genesd1EW=d1EW[,c(5:31)]
n.genesd1EW=dim(genesd1EW)[2]
climated1Ed2TW=d1EW[,c(32:57)]
n.climated1Ed2TW=dim(climated1Ed2TW)[2]

lme.climate={}
for (i in 1:n.genesd1EW) {
  x = genesd1EW[,i]
  for (j in 1:n.climated1Ed2TW) {
    y = climated1Ed2TW[,j]
    lme.climate=rbind(lme.climate, fit.lm(genesd1EW[,i],climated1Ed2TW[,j]))
  }
}
write.csv(lme.climate, file="d1_climate_final_wet.csv")



###############################################################################
# comparing d2 genes and eigen genes to climate and traits from the other day #
###############################################################################

names(d2ED)
pop=d2ED$pop
genesd2ED=d2ED[,c(5:41)]
n.genesd2ED=dim(genesd2ED)[2]
climated2Ed1TD=d2ED[,c(42:55)]
n.climated2Ed1TD=dim(climated2Ed1TD)[2]

lme.climate={}
for (i in 1:n.genesd2ED) {
  x = genesd2ED[,i]
  for (j in 1:n.climated2Ed1TD) {
    y = climated2Ed1TD[,j]
    lme.climate=rbind(lme.climate, fit.lm(genesd2ED[,i],climated2Ed1TD[,j]))
  }
}
write.csv(lme.climate, file="d2_climate_final_dry.csv")

names(d2EW)
pop=d2EW$pop
genesd2EW=d2EW[,c(5:41)]
n.genesd2EW=dim(genesd2EW)[2]
climated2Ed1TW=d2EW[,c(42:55)]
n.climated2Ed1TW=dim(climated2Ed1TW)[2]

lme.climate={}
for (i in 1:n.genesd2EW) {
  x = genesd2EW[,i]
  for (j in 1:n.climated2Ed1TW) {
    y = climated2Ed1TW[,j]
    lme.climate=rbind(lme.climate, fit.lm(genesd2EW[,i],climated2Ed1TW[,j]))
  }
}
write.csv(lme.climate, file="d2_climate_final_wet.csv")


##################
##fdr correction##
##################
install.packages("fdrtool")
library(fdrtool)

d1_p = read.csv("d2_traits_contigs_fdr.csv")
names(d1_p)
p_values = d1_p[,12]
q_values = p.adjust(p_values, method="fdr")
write.csv(q_values, file="d2_traits_contigs_fdr_qvalues.csv")

##################################################
## plotting significant contigs and eigen genes ##
##################################################
library(ggplot2)
library(grid)
library(gridExtra)
library(Rmisc)

d1T = read.csv("voomd1_popxtreatGenesTraits.csv", header=TRUE, na.strings = "NA")
data = read.csv("voomd2_popxtreatGenesTraits.csv", header=TRUE, na.strings = "NA")



listd1 = read.csv("list_ME_d2T_int.csv", header=TRUE)
d1 = as.matrix(listd1)
d1
par(mfrow = c(5, 4))  
plots = list()
for (i in 1:nrow(d1)) {
    p1 = ggplot(data=d2T, aes_string(x=d1[i,1], y=d1[i,2], color="treatment")) + 
    geom_point() +
    scale_colour_hue(l=50) + # Use a slightly darker palette than normal
    geom_smooth(method=lm,   # Add linear regression lines
                se=FALSE)    # Don't add shaded confidence region
    plots[[i]] = p1
    
}
multiplot(plotlist = plots, cols=5)


ggplot(data=d1, aes(x= contig19472 , y=total_above_sugar, color=factor(treatment))) + 
  geom_point() +
  scale_colour_hue(l=50) + # Use a slightly darker palette than normal
  geom_smooth(method=lm,   # Add linear regression lines
              se=FALSE)    # Don't add shaded confidence region


### context dependency figure ####

library(ggplot2)
require(lme4)
library(plyr)

data = read.csv("ME7_elevation_dry.csv", header=TRUE, na.string="NA")
head(data)
names(ind_mean)
dim(ind_mean)
summary(ind_mean)
pops = sapply(data$pop, substr, start=1, stop=1)
data_f = data.frame(data, pop=pops)

z = ddply (data_f, c("pops"), summarize,na.rm =TRUE,
           y = mean(dataME7, na.rm =TRUE),
           x = mean(elevation, na.rm =TRUE),
           sdy= sd(dataME7, na.rm =TRUE),
           sdx= sd(elevation, na.rm =TRUE),
           ly= length(dataME7[!is.na(dataME7)]),
           lx= length(elevation[!is.na(elevation)]),
           ySE = sdy/sqrt(ly),
           xSE = sdx/sqrt(lx))



##Jane's not used
ggplot(data = z,aes(x = x,y = y)) + xlab("Elevation") + ylab("Gene network 7")  + 
  geom_point(shape=20, size=20) +   theme_bw() + xlim(c(100, 1600))  + ylim(c(-0.1, 0.4))+  
  geom_errorbarh(data=z, aes(xmax= x+xSE, xmin= x- xSE))+
  geom_errorbar(data=z, aes(x=x, ymax= y+ySE, ymin=y-ySE)) +  
  geom_hline(yintercept=0, lty=2, size=1, color="grey40") +
  geom_text(data=z, fontface="bold", color="white", mapping=aes(x=x, y=y, label=substr(pops, 1, 1))) +
  geom_smooth(method=lm) + theme_bw()+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=25))

