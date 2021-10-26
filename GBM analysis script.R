#At first load the library's required
library(WGCNA)
library(SNFtool)
library(RColorBrewer)
library(dhga)

# read in the R libraries 
library(MASS) # standard, no need to install 
library(class) # standard, no need to install 
library(cluster) 
library(impute)# install it for imputing missing value


#next we loaded the gbm data txt file using delim function
gbm_data<-read.delim("C:/project_5/raw_data/GLIO_Gene_Expression.txt ", sep="", header=T)

#lets check the dimensions of the gbm_data and we got 12042 entries with 215 columns.
dim(gbm_data) 
dimnames(gbm_data)

#lets keep the gbm_data as gene.exp
gene.exp<-gbm_data

#next we cut the TCGA ids  
colnames_gene.exp <- colnames(gene.exp)
colnames_gene.exp <- strsplit(as.character(colnames_gene.exp), "[.]")
colnames_gene.exp <- do.call("rbind", colnames_gene.exp)
colnames_gene.exp <- as.data.frame(colnames_gene.exp)
colnames_gene.exp_final <- paste(colnames_gene.exp[,1], "-", colnames_gene.exp[,2], "-", colnames_gene.exp[,3], "-", colnames_gene.exp[,4], sep="") 
colnames_gene.exp_final <- as.data.frame(colnames_gene.exp_final)
colnames_gene.exp_final[,2] <- ""
colnames_gene.exp_final [,2] <- colnames(gene.exp)
colnames(colnames_gene.exp_final) <- c("new_sample_ids", "old_sample_ids")
colnames(gene.exp) <- colnames_gene.exp_final$new_sample_ids

#next we boxplot the gene-expression data using boxplot function before normalization of gene.exp data 
boxplot(gene.exp , las=2,ylab= "gene expression of each sample before normalization")

#next we standard normalize the data 
Ngene.exp <- standardNormalization(gene.exp)
# Normalize each column of Ngene.exp to have mean 0 and standard deviation 1.
standardNormalization = function(Ngene.exp) {
  Ngene.exp= as.matrix(Ngene.exp); Ngene.exp <- t(Ngene.exp)
  mean = apply(Ngene.exp, 2, mean)
  sd = apply(Ngene.exp, 2, sd)
  sd[sd==0] = 1
  Ngene.expNorm = t((t(Ngene.exp) - mean) / sd)
  Ngene.expNorm = t(Ngene.expNorm)
  return(Ngene.expNorm)
}

#next we again plot the boxplot the expression data after normalization
boxplot(Ngene.exp , las=2, ylab= "gene expression of each sample after normalization")


#next we do quantile normalization of the normalised data for better outcome 
#quantile normalization, Input: rows where genes and columns should be samples 


quantile_normalisation<-function(df){
  df_rank<-apply(df,2,rank,ties.method="min")
  df_sorted<-data.frame(apply(df,2,sort))
  df_mean<-apply(df_sorted,1,mean)
  
  
  index_to_mean<-function(my_index,my_mean){
    return(my_mean[my_index])
  }
  
  
  df_final<-apply(df_rank,2,index_to_mean,my_mean=df_mean)
  rownames(df_final)<-rownames(df)
  return(df_final)
}


quant_normgene.exp<-quantile_normalisation(Ngene.exp)


#next we again plot the boxplot the expression data after quantile_normalization
boxplot(quant_normgene.exp, las=2, ylab= "gene expression of each sample after quantile_normalization")


#Next we plot heatmaps of the data to check the qauality of data 
#to plot a heatmap we convert data.frame to data.matrix , where we require a matrix to plot a heatmap
#At first we covert a data.frame in to numeric matrix to plot a boxplot 
df1<-data.frame(gene.exp)
df2<-data.matrix(df1)


#now we plot the heatmap using heatmap function with unnormalised data using scale as row 
heatmap(df2 , scale="row")

#now we plot the heatmap using heatmap function with normalized data
heatmap(Ngene.exp , scale="none")

#now we plot the heatmap using heatmap function with quantile normalized data
heatmap(quant_normgene.exp , scale="none")


#WGCNA ANALYSIS#

#transpose of qaunt normgene expression data 

datExpr = as.data.frame(t(gene.exp[ , ]))

#chosing soft-threshold power and checking mean connectivity of genes
## Choose a set of soft-thresholding powers as part of WGCNA 
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

#this line corresponds to using an R^2 cut-off of h
abline(h=0.86,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#ATER all identifications from above plots we decided softpower as 9
softPower<-9


#connectivity plot for chosing better soft-threshold power to get biologically meaningful modules 

#connectivity / scale free plot #######
####scalefreeplot######
#R^2 value must be greather than 0.8 , to statisfy sacale free topology 
Connectivity=softConnectivity(datExpr,power=softPower)-1
# Let's create a scale free topology plot. 
# The black curve corresponds to scale free topology and 
# the red curve corresponds to truncated scale free topology. 
par(mfrow=c(1,1)) 
scaleFreePlot(Connectivity, main=paste("soft threshold, power=",softPower), truncated=F); 



#converting similarity matrix into adjacency matrix
#calclute the adjacency matrix
adj= adjacency(datExpr,type = "signed", power = softPower);

#turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(datExpr,networkType = "signed", TOMType = "signed", power = softPower);

#calculation of dissTOM(dissimilarity Topological overlap matrix)
dissTOM=1-TOM

#next we cluster the data using hclust function
#hierarchial clustering(hierTOM)
hierTOM = hclust(as.dist(dissTOM),method="average"); 
par(mfrow=c(1,1)) 
plot(hierTOM,labels=F)



#The reason is: there is a conflict between WGCNA with other software packages. In another package of WGCNA run r studio has the same functionality as the other. WGCNA has its own function "cor", which namespace "cor" is associated.
#solution for: temporarily re-allocation before using this function
cor <- WGCNA::cor

#next we constructed modules 
net<-blockwiseModules(datExpr, power = 9,TOMType = "signed",
                      minModuleSize = 30,reassignThreshold = 0, 
                      mergeCutHeight = 0.25,numericLabels = TRUE, 
                      pamRespectsDendro = FALSE,
                      saveTOMs = TRUE,
                      saveTOMFileBase = "GBM_dataTOM",verbose = 3)


#by using table(net$colors) function we get number of modules we detected , we detected 3 modules 
table(net$colors)

#end
cor<-stats::cor

#and further plotted the modules which where constructed above 
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]],
                    mergedColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)



# According to our definition, modules correspond to branches of the tree. 
# The question is what height cut-off should be used? This depends on the 
# biology. Large heigth values lead to big modules, small values lead to small 
# but tight modules. 
# In reality, the user should use different thresholds to see how robust the findings are.
# The function cutreeStatistColor colors each gene by the branches that 
# result from choosing a particular height cut-off. 
# GREY IS RESERVED to color genes that are not part of any module. 
# We only consider modules that contain at least 30 genes. 
colorh1= cutreeStaticColor(hierTOM,cutHeight = 0.94, minSize = 30) 


#cluster dendrogram with module(branch) color #
par(mfrow=c(2,1),mar=c(2,4,1,1)) 
plot(hierTOM, main="Cluster Dendrogram", labels=F, xlab="", sub=""); 
plotColorUnderTree(hierTOM,colors=data.frame(module=colorh1)) 
title("Module (branch) color")

#To know how many genes allocated for each module after height cut 
table(colorh1)


# An alternative view of this is the so called TOM plot that is generated by the 
# function TOMplot 
# Inputs: TOM distance measure, hierarchical (hclust) object, color 
# Warning: for large gene sets, say more than 2000 genes 
#this will take a while. I recommend you skip this. 
TOMplot(dissTOM , hierTOM, colorh1)


#caution much time #
# We also propose to use classical multi-dimensional scaling plots 
# for visualizing the network. Here we chose 2 scaling dimensions 
# This also takes about 30 minutes... 
cmd1=cmdscale(as.dist(dissTOM),2) 
par(mfrow=c(1,1)) 
plot(cmd1, col=as.character(colorh1), main="MDS plot",xlab="Scaling Dimension 
1",ylab="Scaling Dimension 2")


# The following produces heatmap plots for each module. 
# Here the rows are genes and the columns are samples. 
# Well defined modules results in characteristic band structures since the corresponding genes are 
# highly correlated. 
#####
sizeGrWindow(8,9)
par(mfrow=c(3,1), mar=c(1, 2, 4, 1))
which.module="turquoise";
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
# for the second (blue) module we use
which.module="blue";
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
which.module="brown";
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
which.module="lightgreen";
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
which.module="red";
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )


which.module="yellow";
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
which.module="grey";
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
which.module="pink";
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
which.module="black";
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
which.module="purple";
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )


# To get a sense of how related the modules are one can summarize each module 
# by its first eigengene (referred to as principal components). 
# and then correlate these module eigengenes with each other. 
datME=moduleEigengenes(datExpr,colorh1)

# We define a dissimilarity measure between the module eigengenes that keeps track of the sign of 
#the correlation between the module eigengenes. 
dissimME=1-(t(cor(datME, method="p")))/2 
hclustdatME=hclust(as.dist(dissimME), method="average" ) 
par(mfrow=c(1,1)) 
plot(hclustdatME, main="Clustering tree based on the module eigengenes of modules")




#Relation between module eigengenes 
#which depicts how modules related with each other
# Now we create scatter plots of the samples (arrays) along the module eigengenes. 

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
#pairs(z,upper.panel=panel.cor)

panel.hist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks
  nB <- length(breaks)
  y <- h$counts
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "white", ...)
}

datMEordered=datME[,hclustdatME$order] 
pairs( datMEordered, upper.panel = panel.smooth, lower.panel = panel.cor , 
       diag.panel=panel.hist ,main="Relation between module eigengenes") 





#network visualization for each module we detected after height cut .
#network viualization 

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 9);
# Read in the annotation file
annot = read.delim("C:/project_5/raw_data/Gene_annotation.csv", sep=",", header=T)
# Select modules
modules = c("turquoise");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(mergedColors, modules));
modProbes = probes[inModule];
modGenes = annot$gene_symbol
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = mergedColors[inModule]);



#network viualization 

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 9);
# Read in the annotation file
annot = read.delim("C:/project_5/raw_data/Gene_annotation.csv", sep=",", header=T)
# Select modules
modules = c("blue");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(mergedColors, modules));
modProbes = probes[inModule];
modGenes = annot$gene_symbol
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = mergedColors[inModule]);


#network viualization 

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 9);
# Read in the annotation file
annot = read.delim("C:/project_5/raw_data/Gene_annotation.csv", sep=",", header=T)
# Select modules
modules = c("brown");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(mergedColors, modules));
modProbes = probes[inModule];
modGenes = annot$gene_symbol
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = mergedColors[inModule]);





#network viualization 

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 9);
# Read in the annotation file
annot = read.delim("C:/project_5/raw_data/Gene_annotation.csv", sep=",", header=T)
# Select modules
modules = c("lightgreen");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(mergedColors, modules));
modProbes = probes[inModule];
modGenes = annot$gene_symbol
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = mergedColors[inModule]);



#network viualization 

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 9);
# Read in the annotation file
annot = read.delim("C:/project_5/raw_data/Gene_annotation.csv", sep=",", header=T)
# Select modules
modules = c("red");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(mergedColors, modules));
modProbes = probes[inModule];
modGenes = annot$gene_symbol
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = mergedColors[inModule]);


#network viualization 

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 9);
# Read in the annotation file
annot = read.delim("C:/project_5/raw_data/Gene_annotation.csv", sep=",", header=T)
# Select modules
modules = c("yellow");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(mergedColors, modules));
modProbes = probes[inModule];
modGenes = annot$gene_symbol
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = mergedColors[inModule]);


#network viualization 

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 9);
# Read in the annotation file
annot = read.delim("C:/project_5/raw_data/Gene_annotation.csv", sep=",", header=T)
# Select modules
modules = c("grey");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(mergedColors, modules));
modProbes = probes[inModule];
modGenes = annot$gene_symbol
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = mergedColors[inModule]);


#network viualization 

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 9);
# Read in the annotation file
annot = read.delim("C:/project_5/raw_data/Gene_annotation.csv", sep=",", header=T)
# Select modules
modules = c("pink");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(mergedColors, modules));
modProbes = probes[inModule];
modGenes = annot$gene_symbol
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = mergedColors[inModule]);



#network viualization 

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 9);
# Read in the annotation file
annot = read.delim("C:/project_5/raw_data/Gene_annotation.csv", sep=",", header=T)
# Select modules
modules = c("green");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(mergedColors, modules));
modProbes = probes[inModule];
modGenes = annot$gene_symbol
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = mergedColors[inModule]);




#network viualization 

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 9);
# Read in the annotation file
annot = read.delim("C:/project_5/raw_data/Gene_annotation.csv", sep=",", header=T)
# Select modules
modules = c("black");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(mergedColors, modules));
modProbes = probes[inModule];
modGenes = annot$gene_symbol
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = mergedColors[inModule]);





#network viualization 

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 9);
# Read in the annotation file
annot = read.delim("C:/project_5/raw_data/Gene_annotation.csv", sep=",", header=T)
# Select modules
modules = c("purple");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(mergedColors, modules));
modProbes = probes[inModule];
modGenes = annot$gene_symbol
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = mergedColors[inModule]);





#Based on relation between module eigengenes
#network viualization 
#At first we plotted the cytoscape network for red and turquoise because thier relation
#of module eigengenes was high( it was 0.86)
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 9);
# Read in the annotation file
annot = read.delim("C:/project_5/raw_data/Gene_annotation.csv", sep=",", header=T)
# Select modules
modules = c("red","turquoise");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(mergedColors, modules));
modProbes = probes[inModule];
modGenes = annot$gene_symbol
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = mergedColors[inModule]);                               






#Based on relation between module eigengenes
#network viualization 
#At first we plotted the cytoscape network for grey  and yellow because thier relation
#of module eigengenes was secound highest( it was 0.80)
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 9);
# Read in the annotation file
annot = read.delim("C:/project_5/raw_data/Gene_annotation.csv", sep=",", header=T)
# Select modules
modules = c("grey","yellow");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(mergedColors, modules));
modProbes = probes[inModule];
modGenes = annot$gene_symbol
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = mergedColors[inModule]);                               



#Based on relation between module eigengenes
#network viualization 
#At first we plotted the cytoscape network for red and turquoise because thier relation
#of module eigengenes was high( it was 0.86)
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 9);
# Read in the annotation file
annot = read.delim("C:/project_5/raw_data/Gene_annotation.csv", sep=",", header=T)
# Select modules
modules = c("green","turquoise");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(mergedColors, modules));
modProbes = probes[inModule];
modGenes = annot$gene_symbol
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = mergedColors[inModule]);                               

#################################################################################################

