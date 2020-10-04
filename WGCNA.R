##########################
## Consensus clustering ##
##########################

library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in d1 data set: thi is from all the data divided in to 2 parts d1 and d2
data_d1 = read.csv("consensus_d1.csv");
data_d2 = read.csv("consensus_d2.csv");
nSets = 2;
setLabels = c("day1", "day2")
shortLabels = c("d1", "d2")
# Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = as.data.frame(t(data_d1[-c(1)])));
names(multiExpr[[1]]$data) = data_d1$genes;
rownames(multiExpr[[1]]$data) = names(data_d1)[-c(1)];
multiExpr[[2]] = list(data = as.data.frame(t(data_d2[-c(1)])));
names(multiExpr[[2]]$data) = data_d1$genes;
rownames(multiExpr[[2]]$data) = names(data_d2)[-c(1)];
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)
multiExpr.mat = mtd.apply(multiExpr, as.matrix)

gsg = goodSamplesGenesMS(multiExpr.mat, verbose = 3); 
gsg$allOK

if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
      2
    printFlush(paste("In set", setLabels[set], "removing samples",
                     paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}
sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}

pdf(file = "Consensus_SampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
dev.off();
# Define data set dimensions
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;
save(multiExpr, nGenes, nSamples, setLabels, shortLabels,
     file = "Consensus-dataInput.RData");

### consensus clustering on the farm cluster ###
# ssh to farm and then do
srun -p bigmemh --pty --mem 81920 bash
module load R
R


library(WGCNA)
enableWGCNAThreads()
lnames = load(file = "Consensus-dataInput.RData")
lnames
nSets = checkSets(multiExpr)$nSets
# Choose a set of soft-thresholding powers
powers = c(seq(1,10,by=1), seq(12,20, by=2))
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]]);
collectGarbage();
write.table(sft,file="sft_d2.csv", sep="," )

# to select sft treshold check, 4 or 6 appropiate, let's go with 6

# constructing networks in a blockwise manner for large datasets

bnet = blockwiseConsensusModules(
  multiExpr, maxBlockSize = 25000, power = 4, minModuleSize = 30,
  deepSplit = 2, 
  pamRespectsDendro = FALSE, 
  mergeCutHeight = 0.25, numericLabels = TRUE,
  minKMEtoStay = 0,
  saveTOMs = TRUE, verbose = 5)

consMEs = bnet$multiMEs;
moduleLabels = bnet$colors;
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
consTree = bnet$dendrograms[[1]]; 
save(consMEs, moduleLabels, moduleColors, consTree, bwLabels, bwColors, file = "Consensus-NetworkConstruction-final.RData")
##
save(bnet, file="bnet_consensus.RData")
geneTree = bwnet$dendrograms[[1]];
colors=bnet$colors
moduleColors=colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1;
unmergedColors=bnet$unmergedcolors
MEs=bnet$MEs
goodSamples=bnet$goodSamples
goodGenes=bnet$goodGenes
dendrograms=bnet$dendrograms
TOMFiles=bnet$TOMFiles
blockGenes=bnet$blockGenes
blocks=bnet$blocks


save(consMEs, moduleLabels, moduleColors, consTree,colors,moduleLabels,bwColors,bwLabels, geneTree,unmergedColors,MEs,goodSamples,goodGenes,dendrograms,
     TOMFiles,blockGenes,blocks, file = "networkConstruction_consensus.RData")
q()
n
exit

## Visulaization of consensus modules ##
lnames=load(file="bnet_consensus.RData")
lnames=load(file="Consensus-NetworkConstruction-final.RData")

names(bnet)

sizeGrWindow(12,6)
# Use the layout function for more involved screen sectioning
layout(matrix(c(1:4), 2, 2), heights = c(0.8, 0.2), widths = c(1,1))
layout.show(4);
nBlocks = length(bnet$dendrograms)

# Plot the dendrogram and the module colors underneath for each block (I am not sure how the colors are selected since colors are not present in bnet)
for (block in 1:nBlocks)
  plotDendroAndColors(bnet$dendrograms[[block]], moduleColors[bnet$blockGenes[[block]]],
                      "Module colors", 
                      main = paste("Gene dendrogram and module colors in block", block), 
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      setLayout = FALSE)


###GO analysis done by GOrilla###

##Exporting results##

lnames = load(file = "Consensus-dataInput.RData");
nGenes
lnames=load(file="Consensus-NetworkConstruction-final.RData")
lnames
# Read in the probe annotation
annot = read.csv(file = "ATannotation.csv");
nSets=2
# Match probes in the data set to the probe IDs in the annotation file
probes = names(multiExpr[[1]]$data)
probes2annot = match(probes, annot$queryID)
allLLIDs = annot$subjectID[probes2annot];

consMEs.unord = multiSetMEs(multiExpr, universalColors = moduleLabels, excludeGrey = TRUE)
GS = list();
kME = list();
for (set in 1:nSets)
{
  kME[[set]] = corAndPvalue(multiExpr[[set]]$data, consMEs.unord[[set]]$data);
}
kME.metaZ = (kME[[1]]$Z + kME[[2]]$Z)/sqrt(2);
kME.metaP = 2*pnorm(abs(kME.metaZ), lower.tail = FALSE);
# Same code for kME:
kMEmat = rbind(kME[[1]]$cor, kME[[2]]$cor, kME[[1]]$p, kME[[2]]$p, kME.metaZ, kME.metaP);
MEnames = colnames(consMEs.unord[[1]]$data);
nMEs = checkSets(consMEs.unord)$nGenes
dim(kMEmat) = c(nGenes, 6*nMEs)
rownames(kMEmat) = probes;
colnames(kMEmat) = spaste(
  c("kME.set1.", "kME.set2.", "p.kME.set1.", "p.kME.set2.", "Z.kME.meta.", "p.kME.meta"),
  rep(MEnames, rep(6, nMEs)))
info = data.frame(Probe = probes, GeneSymbol = allLLIDs,
                  ModuleLabel = moduleLabels,
                  ModuleColor = labels2colors(moduleLabels),
                  kMEmat);
write.csv(info, file = "consensusAnalysis-CombinedNetworkResults.csv",
          row.names = FALSE, quote = FALSE);



# names (colors) of the modules
modNames = substring(names(bnet$multiMEs), 3)
geneModuleMembership = as.data.frame(cor(multiExpr[[1]], MEs, use = "p"));
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
annot= read.csv(file= "ATAnnotation.csv")
dim(annot)
names(annot)
gene=names(datExpr)
gene2AT=match(gene, annot$queryID)
sum(is.na(gene2AT))

geneInfo0 = data.frame(genes=names(datExpr),LocusLinkID = annot$subjectID[gene2AT],
                       moduleColor = bwnet$colors)
# Order modules by their significance 
modOrder = order(-abs(cor(bwnet$MEs, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership,
                         MMPvalue);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames, sep=""),
                       paste("p.MM.", modNames, sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor);
geneInfo = geneInfo0[geneOrder, ]
head(geneInfo)
write.csv(geneInfo, file="geneInfo_all.csv")



#### Plotting top 50 genes with smallest p values for meta consesnsus ME membership###

data = read.csv("consensus_d2MM_average_ME7_ele.csv")
dim(data)
names(data)
head(data)
genes=data$genes
counts= data[,2:8]
#all genes
EG =read.csv("ME7_set2_abs_top100.csv", header=FALSE)
EG=EG$V1
EG
EGcounts= counts[match(EG, genes),]
EGcounts

EGcounts.t=t(EGcounts)
EG.scaled= scale(abs(EGcounts.t), center=TRUE, scale=apply(EGcounts.t,2,sd))

head(EG.scaled)
names(EG.scaled)
library(scales)
matplot(EG.scaled, pch=20, cex=0.2,type="b", col=alpha("black", 0.3), ylab="Scaled Gene Expression (EG7 module d2)", xaxt="n", xlab=NA)
points(c(EG.scaled[,1]), col="red", cex=1.75, pch=16)
lines(1:7,EG.scaled[,1], col="red2", cex=1.75, pch=16, lwd=2)
text(seq(1,14,by=1), par("usr")[3]-0.5,srt=90, cex=1, pos=1, xpd=TRUE,labels=substr(names(EGcounts),11,13))

#### Plotting all gene network averages###

data = read.csv("consensus_d2MM_average.csv")
dim(data)
names(data)
head(data)
genes=data$genes
counts= data[,2:15]
#all genes
EG =read.table("GN6_first50Genes.txt", header=FALSE)
EG=EG$V1
EG
EGcounts= counts[match(EG, genes),]
EGcounts

EGcounts.t=t(EGcounts)
EG.scaled= scale(abs(EGcounts.t), center=TRUE, scale=apply(EGcounts.t,2,sd))

head(EG.scaled)
names(EG.scaled)
library(scales)
matplot(EG.scaled, pch=20, cex=0.2,type="b", col=alpha("black", 0.1), ylab="Scaled Gene Expression (GN6 module d2)", xaxt="n", xlab=NA)
points(c(EG.scaled[,1]), col="red", cex=1.75, pch=16)
lines(1:7,EG.scaled[,1], col="red2", cex=1.75, pch=16, lwd=2)
text(seq(1,14,by=1), par("usr")[3]-0.5,srt=90, cex=1, pos=1, xpd=TRUE,labels=substr(names(EGcounts),11,13))
