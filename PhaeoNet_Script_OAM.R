
##########################################################
# PhaeoNet creation script
# Ait-Mohamed Ouardia 
# This script describes major steps in the PhaeoNet creation 
# for the paper 
# PhaeoNet: A holistic RNA-seq-Based Portrait of Transcriptional Coordination in 
# the Model Diatom Phaeodactylum tricornutum
# frontiers in Plant Science
#doi: 10.3389/fpls.2020.590949
##########################################################

#########################
# Call required Libraries 
#########################
library("DESeq2")
library("ggplot2")
library("gplots")
library("dplyr")
library("ggrepel")
library("flashClust")
library("genefilter")
library("RColorBrewer")
library("BiocParallel")
library("data.table")
library("limma")
library("edgeR")
library("class") 
library("cluster")   
library("Hmisc")
library("impute")
library("WGCNA")
library("FactoMineR")
library("amap")
library("pamr")
library("bladderbatch")
library("ggfortify")
library("Biobase")
options(stringsAsFactors=FALSE)
disableWGCNAThreads()


#=====================================================================================
#
#  Code part 1
#
#=====================================================================================
#---------------------------------------------------------------------------------------#
# Create the counting matrix                                                            #
#---------------------------------------------------------------------------------------#

#### Read the PhaeoNet count files 
# The folder contains files, each file is a two column table
# The first column is gene IDs and the second column is gene counts

basedir <- ("/export/home/users/mpb/Dana_Scully")
#setwd(basedir)
cntdir <- paste(basedir, "Count_files", sep="/")
pat <- ".txt"
tophat.all <- list.files(path = cntdir,
                         pattern = pat,
                         all.files = TRUE,
                         recursive = FALSE,
                         ignore.case = FALSE,
                         include.dirs = FALSE)

# we take the 'all' series
myfiles <- tophat.all
DT <- list()

# I Read each file as an array element of DT and rename the last 2 colonnes
# I create a list of single sample tables
for (i in 1:length(myfiles) ) {
  infile = paste(cntdir, myfiles[i], sep = "/")
  DT[[myfiles[i]]] <- read.table(infile, header = F, stringsAsFactors = FALSE)
  cnts <- gsub("(.*)_all_counts.txt", "\\1", myfiles[i])
  colnames(DT[[myfiles[i]]]) <- c("ID", cnts)
}

# Merge all files based on the first ID columns
data <- DT[[myfiles[1]]]

# Check that I have all the tables
head(data)
for (i in 2:length(myfiles)) {
  y <- DT[[myfiles[i]]]
  z <- merge(data, y, by = c("ID"))
  data <- z
}

# ID column becomes rownames
rownames(data) <- data$ID
data <- data[,-1]

# Add the total counts per sample
data <- rbind(data, tot.counts=colSums(data))

# Inspect and look at the top row names of the table
head(data)

# Create a table organized  like in the colData file
rawCountTable <- data
names(rawCountTable) = gsub(pattern = ".txt", replacement = "", x = names(rawCountTable))
cts <- as.matrix(rawCountTable)

# remove the last row of the counting matrix that contains the column sum of reads
cts <- cts[-nrow(cts),] 

######This part ends the preparation of gene matrix counts

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================
#---------------------------------------------------------------------------------------#
# Create the DDS object from DESeq2 pckage                                              #
#---------------------------------------------------------------------------------------#

# The DDS object will help to create all quality check plots of data 
# Two variables were used in the design : Batch effect and the Treatment
############## DDS object ##################
# Read the design file table
# Sample names in the counting matrix table MUST be the same as in the coldata.csv file
coldata <- read.csv("/export/home/users/Dana_Scully/Design_file/design_file_Phaeonet.csv")
#DESeq2 object creation
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Batch + Treatment)


dds$Treatment <- relevel(dds$Treatment, ref = "Ctl")
dds$group <- factor(paste0(dds$Batch, dds$Treatment))
design(dds) <- ~ group

#### Remove samples with a median lower than 10 #######

idx <- dds[rowMedians(norm_counts)>=10,]
dds <- DESeq(idx)

#=====================================================================================
#
#  Code part 3
#
#=====================================================================================
#---------------------------------------------------------------------------------------#
# Check the quality of Data                                                             #
#---------------------------------------------------------------------------------------#

#### Choose the DESeq2 vst normalisation due to the high number of samples
vsd <- vst(dds, blind=FALSE)


#### Margins and window sizes are given all through the script as 
####indicative, please change them according to your code
# Plot PhaeoNet data distribution 

par(mar=c(9,5,2,2))
p1 <- plotDispEsts(dds)

vsdPhaeo <- assay(vsd)
# Plot PhaeoNet vst gene counts Spearman correlation 
p2 <- plot(hcluster(dist(t(vsdPhaeo)), method= "spearman", link="average"), 
           cex=0.5, cex.main=1, main= "PhaeoNet Hierarchical Clustering based on Spearman correlation")


#PCA plot of PhaeoNet data to visualise the Batch effect

p3 <- plotPCA(vsd, intgroup=c("Batch", "Treatment"), returnData=TRUE)

#Heatmap plot based on Euclidean distance of PhaeoNet data to visualise the Batch effect

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Batch,
                                    vsd$Treatment,
                                    vsd$Replicate, sep="-")
colnames(sampleDistMatrix) <- paste(vsd$Batch,
                                    vsd$Treatment,
                                    vsd$Replicate, sep="-")

colours = colorRampPalette(rev(brewer.pal(9, "OrRd")))(256)
p4 <- heatmap.2(sampleDistMatrix, trace="none", col=colours, margins = c(9, 9)) 


#=====================================================================================
#
#  Code part 4
#
#=====================================================================================

#---------------------------------------------------------------------------------------#
# Remove the batch effect and recheck the quality of Data                               #
#---------------------------------------------------------------------------------------#

# Call the function of removeBatchEffect from Limma
# Several functions from other packages have been tested, we kept the removeBatchEffect
# As it gave the best results

vsd.removed=removeBatchEffect(vsd, batch=coldata$Batch, design= design)

# Replace the VST normalised count table with the one without the batch effect
assay(vsd) <- vsd.removed

#PCA plot of PhaeoNet data to visualise the data after Batch effect removal
p5 <- plotPCA(vsd, intgroup=c("Batch", "Treatment"))

#Heatmap plot based on Euclidean distance of PhaeoNet after Batch effect removal
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Treatment,
                                    vsd$Batch,
                                    vsd$Replicate, sep="-")
colnames(sampleDistMatrix) <- paste(vsd$Treatment,
                                    vsd$Batch,
                                    vsd$Replicate, sep="-")


colours = colorRampPalette(rev(brewer.pal(9, "OrRd")))(256)
p6 <- heatmap.2(sampleDistMatrix, trace="none", col=colours, margins = c(9, 9)) 

#=====================================================================================
#
#  Code part 5
#
#=====================================================================================

#---------------------------------------------------------------------------------------#
# Prepare the count matrix relived from batch effect to WGCNA modules creation          #
#---------------------------------------------------------------------------------------#

#####Change the column names and put $Symbol for the gene column

Phaeo = read.table("/home/users/Dana_Scully/Combat_normalized_data_PheoNet.txt", header = TRUE, sep = "\t")
datExpr0 = as.data.frame(t(Phaeo[, -1]))
names(datExpr0) = Phaeo$Symbol
rownames(datExpr0) = names(Phaeo)[-1]

#####Check the quality of the total gene count table

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

#=====================================================================================
#
#  Code part 6
#
#=====================================================================================

#####Create sample tree
sampleTree = hclust(dist(datExpr0), method = "average")
#####Plot the sampleTree
# Window dimensions should be changed if the window is too large or too small.
sizeGrWindow(12,9)
par(mar = c(0,4,2,0))
p7 <- plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1, 
      cex.axis = 1, cex.main = 1.5)
# Plot a line to show the cut
abline(h = 150, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 150, minSize = 10)
# The first cluster contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

######This part ends the preparation of gene counts for WGCNA modules creation

#=====================================================================================
#
#  Code part 7
#
#=====================================================================================

#---------------------------------------------------------------------------------------#
# Soft threshold power choose                                                           #
#---------------------------------------------------------------------------------------#


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
p8 <- plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
      main = paste("Scale independence"));
      text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
p9 <- plot(sft$fitIndices[,1], sft$fitIndices[,5],
      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
      main = paste("Mean connectivity"))
      text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# PhaeoNet was built with a softPowerThreshold of 12
beta1=12
k.dataOne=softConnectivity(datExpr,power=beta1)-1 
par(mfrow=c(2,2))
p10 <- scaleFreePlot(k.dataOne, main=paste("Phaeo Data, power=",beta1), truncated=F) 

#=====================================================================================
#
#  Code part 9
#
#=====================================================================================

#---------------------------------------------------------------------------------------#
# Step by step PhaeoNet network construction                                            #
#---------------------------------------------------------------------------------------#

# Create adjacency matrix
softPower = 12;
adjacency = adjacency(
  datExpr, 
  power = softPower, 
  type = "signed", 
  corFnc = "bicor")

# From adjacency matrix to topological overlap matrix
TOM =TOMsimilarityFromExpr(
  datExpr,
  networkType = "signed",
  power = softPower,
  TOMType = "none",
  corType = "bicor",
  maxPOutliers = 0.05)

# Create a dissimilarity matrix
dissTOM = 1-TOM
# Plot the TOM matrix set to power 7 to get rid from lower corrleations 
# This step requires memory and time to be computed 
p11 <- plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;

# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), method = "average")

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
p12 <- plot(geneTree, xlab="", sub="", main = " PhaeoNet gene clustering on TOM-based dissimilarity",
      labels = FALSE, hang = 0.04);

############# Step by step PhaeoNet network construction ####################

# For PhaeoNet, we set the minimum module size relatively high (40)
minModuleSize = 40;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
module_table <- table(dynamicColors)
# Save a table that contains the modules with their corresponding gene numbers 
write.table(as.data.frame(module_table), row.names = FALSE, sep= "\t", file="/export/home/users/Dana_Scully/step_by_step_network_PhaeoNet_modules.txt")

# Save the modules with their corresponding gene
intModules <- names(module_table)
probes = names(datExpr)
for (module in intModules)
{
  # Select module probes
  modGenes = (dynamicColors== module)
  # Get their entrez ID codes
  modLLIDs = probes[modGenes]
  # Write them into a file
  fileName = paste("/export/home/users/Dana_Scully/step_by_step_Phaeonet_modules_Gene_IDs-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = "geneID")
}

# Plot the dendrogram that shows PhaeoNet modules 
sizeGrWindow(10,6)
p13 <- plotDendroAndColors(geneTree, dynamicColors, "PhaeoNet Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "PhaeoNet gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 15)
p14 <- plot(METree, main = " PhaeoNet clustering of module eigengenes",
      xlab = "", sub = "")

## MDS plot   
cmd1 = cmdscale(as.dist(dissTOM), 2)
par(mfrow = c(1,1))
p15 <- plot(cmd1, col= dynamicColors, main="", xlab ="PhaeoNet[,1]", ylab="PhaeoNet[,2]")


#=====================================================================================
#
#  Code part 10
#
#=====================================================================================

#---------------------------------------------------------------------------------------#
# Creation of PhaeoNet merged modules                                                   #
#---------------------------------------------------------------------------------------#
############# Creation of PhaeoNet merged modules ####################
# In this section the script shows the steps to create the PhaetoNet merged modules
# The threshold has been set to 0.25, other thresholds have been tested and we decided 
# to work with a threshold of 0.25 

################Merge modules with 0.25 threshold

# Set the merging threshold
MEDissThres = 0.25

# Call the merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The PhaeoNet merged module colors
mergedColors = merge$colors;

# Eigengenes of PhaeoNet merged modules:
mergedMEs = merge$newMEs;

###############Create a tree of modules after the clustering###################

nMEList = moduleEigengenes(datExpr, colors = mergedColors)
nMEs = nMEList$eigengenes
# Calculate dissimilarity of PhaeoNet merged module eigengenes
nMEDiss = 1-cor(nMEs);
# Cluster PhaeoNet merged module eigengenes
nMETree = hclust(as.dist(nMEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 15)
p16 <- plot(nMETree, main = " Phaeo clustering of merged module eigengenes after merging",
      xlab = "", sub = "")


################Create a color table of PhaeoNet merged module colors
p17 <- plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                          c("Dynamic Tree Cut", "PhaeoNet Merged modules"),
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05)

moduleColors = mergedColors
# Visualize PhaeoNet merged modules
table(moduleColors)
# Save a table that contains PhaeoNet merged modules with their corresponding gene numbers 
write.table(as.data.frame(table(moduleColors)), row.names = FALSE, sep= "\t", file="/export/home/users/Dana_Scully/step_by_step_PhaeoNet_merged_modules_after_merging_0.25.txt")
# Save thePhaeoNet merged modules with their corresponding gene IDs
intModules <- names(moduleColors)
probes = names(datExpr)

for (module in intModules)
{
  # Select module probes
  modGenes = (mergedColors== module)
  # Get their entrez ID codes
  modLLIDs = probes[modGenes]
  # Write them into a file
  fileName = paste("/export/home/users/Dana_Scully/merged_modules_0.25/step_by_step_Phaeonet_merged_modules_Gene_IDs-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = "geneID")
}

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(mergedColors, colorOrder)-1;
MEs = mergedMEs;

###############Create a tree of PhaeoNet merged modules ###################
nMEList = moduleEigengenes(datExpr, colors = mergedColors)
nMEs = nMEList$eigengenes
# Calculate dissimilarity of module eigengenes
nMEDiss = 1-cor(nMEs);
# Cluster module eigengenes
nMETree = hclust(as.dist(nMEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 15)
p18 <- plot(nMETree, main = " Phaeo clustering of module eigengenes after merging",
     xlab = "", sub = "")

###############Create heatmaps and Histograms of PhaeoNet merged modules ###################

  for (i in c(1:length(mergedColors)))
{
  whichmodule=colorlevels[[i]]
  ME=MEs[, paste("ME",which.module, sep="")]
  par(mfrow=c(2,1), mar=c(0.2, 5.5, 2.6, 2))
  plotMat(t(scale(datExpr[,mergedColors==which.module ]) ),
          nrgcols=30,rlabels=T,rcols=which.module,
          main= paste("Module",which.module, sep=" "), cex.main=2)
  par(mar=c(5, 4.2, 0, 0.7))
  barplot(ME, col=which.module, main="", cex.main=2,
          ylab="Eigengene expression",xlab="Samples")
}

###############Create PhaeoNet merged modules eigengenes pearson correlation heatmap and tree ###################

# Heatmap
par(cex = 1.0)
p19 <- plotEigengeneNetworks(MEs, "PhaeoNet Eigengene merged modules Pearson correlation heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
# Tree

p20 <- plotEigengeneNetworks(MEs, "PhaeoNet Eigengene merged modules Pearson correlation tree", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)

#=====================================================================================
#
#  Code part 11
#
#=====================================================================================

#---------------------------------------------------------------------------------------#
# Export pf the paleturquoise merged module to Cytoscape for visualization              #
#---------------------------------------------------------------------------------------#
# This section gives the script to export the paleturquoise PhaeoNet merged module 
# Cytoscape tables based on a threshold of 0.2 for connectivity

# Select the modue paleturquoise
modules = c("paleturquoise")
# Select module probes;
probes = names(datExpr)
inModule = is.finite(match(mergedColors, modules))
modProbes = probes[inModule]

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("/export/home/users/Dana_Scully/Cytoscape/CytoscapeInput-edges-02_dissTOM-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("/export/home/users/Dana_Scully/Cytoscape/CytoscapeInput-nodes-02_dissTOM-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.2,
                               nodeNames = modProbes,
                               altNodeNames = modProbes,
                               nodeAttr = dynamicColors[inModule])

