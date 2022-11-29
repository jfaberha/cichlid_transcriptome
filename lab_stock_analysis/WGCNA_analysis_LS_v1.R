#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("GO.db")
#BiocManager::install("impute")
#install.packages("WGCNA")
library("WGCNA")
library(gridExtra)
library(emmeans)
library(lmtest)
library(ggplot2)
library(limmaDE2)
library("SimSeq")
library("plyr")
library(doParallel);

setwd("/Volumes/Renn_RNAseq_2017/transcriptome")

options(stringsAsFactors = FALSE);

#Load gene(/transcript) count matrix and labels
countData_LS <- read.csv("gene_count_matrix_filt_refseq_vst_LS.csv", row.names="gene_id")
colData_LS <- read.csv("PHENO_DATA_filt_LS_wgcna.csv")

dim(countData_LS)
names(countData_LS)

# We work with one set:
nSets = 1;
# For easier labeling of plots, create a vector holding descriptive names of the sets.
setLabels = c("Lab stock")
shortLabels = c("LS")
# Form multi-set expression data
multiExprLS = vector(mode = "list", length = nSets)
multiExprLS[[1]] = list(data = as.data.frame(t(countData_LS)));
names(multiExprLS[[1]]$data) = row.names(countData_LS);
rownames(multiExprLS[[1]]$data) = names(countData_LS);
exprSizeLS = checkSets(multiExprLS)

# Check that all genes and samples have sufficiently low numbers of missing values.
gsgLS = goodSamplesGenesMS(multiExprLS, verbose = 3);
gsgLS$allOK

#if outlier filtering needed, this will remove from dataset
if (!gsgLS$allOK){
  # Print information about the removed genes:
  if (sum(!gsgLS$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExprLS[[1]]$data)[!gsgLS$goodGenes],
                                              collapse = ", ")))
  for (set in 1:exprSizeLS$nSets)
  {
    2
    if (sum(!gsgLS$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExprLS[[set]]$data)[!gsgLS$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExprLS[[set]]$data = multiExprLS[[set]]$data[gsgLS$goodSamples[[set]], gsgLS$goodGenes];
  }
  # Update exprSizeLS
  exprSizeLS = checkSets(multiExprLS)
}

#for backup
multiExprLS2 <- multiExprLS


#now cluster samples for each set
multiExprLS <- multiExprLS2#reload original dataset if filtering errors occur
sampleTrees = list()
for (set in 1:nSets){
  sampleTrees[[set]] = hclust(dist(multiExprLS[[set]]$data), method = "average")
}

#check cluster trees for outliers, if outliers identified check tutorial for sample trimming.
pdf(file = "Plots/SampleClusteringLS.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
dev.off();

# Choose the "base" cut height for the female data set
baseHeight = 100
# Adjust the cut height for the male data set for the number of samples
cutHeights = c(16, 16*exprSizeLS$nSamplesLS[2]/exprSizeLS$nSamplesLS[1]);
# Re-plot the dendrograms including the cut lines
pdf(file = "Plots/SampleClustering_recut.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets){
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
  abline(h=cutHeights[set], col = "red");
}
dev.off();

for (set in 1:nSets){
  # Find clusters cut by the line
  labels = cutreeStatic(sampleTrees[[set]], cutHeight = cutHeights[set])
  # Keep the largest one (labeled by the number 1)
  keep = (labels==1)
  multiExprLS[[set]]$data = multiExprLS[[set]]$data[keep, ]
}
collectGarbage();
# Check the size of the leftover data
exprSizeLS = checkSets(multiExprLS)
exprSizeLS


############################################
## Skip to here if filtering unnecessary ###
############################################

#Process with sample metadata
dim(colData_LS)
names(colData_LS)

# remove columns that hold information we do not need.
allTraitsLS = colData_LS
# See how big the traits are and what are the trait and sample names
dim(allTraitsLS)
names(allTraitsLS)
allTraitsLS$sampleID
# Form a multi-set structure that will hold the clinical traits.
TraitsLS = vector(mode="list", length = nSets);
for (set in 1:nSets){
  setSamples = rownames(multiExprLS[[set]]$data);
  traitRows = match(setSamples, allTraitsLS$sampleID);
  TraitsLS[[set]] = list(data = allTraitsLS[traitRows, -1]);
  rownames(TraitsLS[[set]]$data) = allTraitsLS[traitRows, 1];
}
collectGarbage();
# Define data set dimensions
nGenesLS = exprSizeLS$nGenesLS;
nSamplesLS = exprSizeLS$nSamplesLS;

save(multiExprLS, TraitsLS, nGenesLS, nSamplesLS, setLabels, shortLabels, exprSizeLS,
     file = "WGCNA-dataInput_LS.RData");


#####################################
# Time to actually do some analysis #
#####################################

# Load the data saved in the first part
lnames = load(file = "WGCNA-dataInput_LS.RData");

# The variable lnames contains the names of loaded variables.
lnames
# Get the number of sets in the multiExprLS structure.
nSets = checkSets(multiExprLS)$nSets


# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets) powerTables[[set]] = list(data = pickSoftThreshold(multiExprLS[[set]]$data, powerVector=powers,
                                                                        networkType = "signed", verbose = 2)[[2]]);
collectGarbage();
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets){
  for (col in 1:length(plotCols)){
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
pdf(file = "Plots/scaleFreeAnalysis_LS.pdf", wi = 8, he = 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets){
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off();



netLS = blockwiseConsensusModules(
  multiExprLS, power = 12, minModuleSize = 30, deepSplit = 0,
  pamRespectsDendro = FALSE, maxBlockSize = 26000, 
  mergeCutHeight = 0.25, numericLabels = TRUE,
  minKMEtoStay = 0, networkType = "signed",
  saveTOMs = TRUE, verbose = 5)

names(netLS)
netLS$individualTOMInfo

consMEsLS = netLS$multiMEs;
moduleLabelsLS = netLS$colors;
# Convert the numeric labels to color labels
moduleColorsLS = labels2colors(moduleLabelsLS)
consTreeLS = netLS$dendrograms[[1]];

cbind(moduleLabelsLS,data.frame(moduleColorsLS))[,1:2]

sizeGrWindow(8,6);
pdf(file = "Plots/ConsensusDendrogramSigned-LS_only.pdf", wi = 12, he = 6)
plotDendroAndColors(consTreeLS, moduleColorsLS,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()

write.csv(consMEsLS, file = "consensusAnalysis-MEs_LS.csv",
          row.names = TRUE, quote = FALSE);

write.csv(cbind(moduleLabelsLS,data.frame(moduleColorsLS)), file = "consensusAnalysis-Gene-Module_key_LS.csv",
          row.names = TRUE, quote = FALSE);

# Define numbers of genes and samples
nGenesLS = ncol(netLS);
nSamplesLS = nrow(netLS);

# Set up variables to contain the module-trait correlations
moduleTraitCorLS = list();
moduleTraitPvalueLS = list();
# Calculate the correlations
moduleTraitCorLS = cor(consMEsLS[[1]]$data, TraitsLS[[1]]$data, use = "p");
moduleTraitPvalueLS = corPvalueFisher(moduleTraitCorLS, exprSizeLS$nSamples);

# Convert numerical lables to colors for labeling of modules in the plot
MEColors = labels2colors(as.numeric(substring(names(consMEsLS[[1]]$data), 3)));
MEColorNames = paste("ME", MEColors, sep="");
# Open a suitably sized window (the user should change the window size if necessary)
sizeGrWindow(10,7)
pdf(file = "Plots/ModuleTraitRelationships_LS.pdf", wi = 10, he = 15);
# Plot the module-trait relationship table for set number 1
set = 1
textMatrix = paste(signif(moduleTraitCorLS, 2), "\n(",
                   signif(moduleTraitPvalueLS, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCorLS)
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCorLS,
               xLabels = names(TraitsLS[[1]]$data),
               yLabels = MEColorNames,
               ySymbols = row.names(moduleTraitCorLS),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module--trait relationships in", setLabels[set]))
dev.off()

write.csv(moduleTraitCorLS, file = "moduleTraitCorLS-MEs_LS.csv",
          row.names = TRUE, quote = FALSE);

write.csv(moduleTraitPvalueLS, file = "moduleTraitPvalueLS-MEs_LS.csv",
          row.names = TRUE, quote = FALSE);

### LRTs and Pairwise contrasts for LS ###
colData_LS_LRT <- read.csv("PHENO_DATA_filt_LS.csv")
              
consMEs_LS <- merge(data.frame(netLS$multiMEs), data.frame(colData_LS_LRT), 
                    by.x = 0, by.y = 1, all = TRUE, na.omit = TRUE)
rownames(consMEs_LS) <- consMEs_LS[,1]
consMEs_LS[,1] <- NULL

## Look at 2-way interaction LRTs
nModules <- ncol(consMEs_LS)-ncol(colData_LS_LRT)+1
LRTresults_2wayLS = vector(mode="list", length = nModules)
for (ME in 1:nModules){
  nested <- lm(consMEs_LS[[ME]]~1,data=consMEs_LS)
  complex <- lm(consMEs_LS[[ME]]~treat,data=consMEs_LS)
  LRTresults_2wayLS[[ME]] = lrtest(nested, complex)
  names(LRTresults_2wayLS)[[ME]] <- colnames(consMEs_LS)[[ME]]
}

# Fit regression model for 2-way interaction and look at 2-way contrasts
lm_2wayLS = vector(mode="list", length = nModules)
contrasts_2wayLS = vector(mode="list", length = nModules)
for (ME in 1:nModules){
  lm_2wayLS[[ME]] <- lm(consMEs_LS[[ME]]~treat,data=consMEs_LS)
  contrasts_2wayLS[[ME]] = contrast(emmeans(lm_2wayLS[[ME]], ~ treat), interaction = "pairwise")
  names(contrasts_2wayLS)[[ME]] <- colnames(consMEs_LS)[[ME]]
}

# Make boxplots
for (i in colnames(consMEs_LS[1:nModules])){
  pdf(paste("WGCNA_outlierFilt/plots/LS_signed_power12_", i, "_outFilt.pdf",sep=""), width = 6, height = 4) # Open a new pdf file
  print(ggplot(data = consMEs_LS, aes(x=treat, y=unlist(consMEs_LS[i]))) + geom_boxplot(aes(fill=treat)) + ylab(i) + xlab("time-point"))
  dev.off() # Close the file
}
