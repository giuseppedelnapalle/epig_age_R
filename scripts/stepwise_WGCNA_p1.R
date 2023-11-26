#!/usr/bin/env Rscript
# $title stepwide WGCNA part_1
# $description
#   stepwise weighted gene correlation network analysis
#   adapted from "Corrected R code from chapter 12 of the book" of Horvath
# $input
#   rna.hg38
# $main_output
#   m_eigen, gene.w, moduleColorsManual2, phenoDT.rna
# $author giuseppe
# $date Apr 2020

library(DESeq2)
library(flashClust)
options(stringsAsFactors = F)

# section_1 preprocess RNA data -------------------------------------------

# 1 gene filtering --------------------------------------------------------

# removing all features that have a count of less than n in more than p (proportion) of the samples
dat_rna0 <- assay(rna.hg38)
dim(dat_rna0)
# obtain data of tumour sample whose DNAm data available
smp_r <- sample_name(colnames(dat_rna0), n.field = 4)
smp_m <- sample_name(phenoDT.met$sample, n.field = 4)
dat_sub <- dat_rna0[,smp_r %in% smp_m]
dim(dat_sub)

# subset top 1:n_sub genes by variance ------------------------------------

n_sub <- round(nrow(dat_sub)*top_p)
variance <- apply(dat_sub, 1, var)
#summary(variance)
dat_filt <- dat_sub[order(variance, decreasing = T)[1:n_sub],]

# define background genes -------------------------------------------------

top_p2 <- .75
n_sub2 <- round(nrow(dat_sub)*top_p2)
var_p2 <- variance[order(variance, decreasing = T)][1:n_sub2]
var_thres <- var_p2[length(var_p2)]
dat_filt2 <- dat_sub[order(variance, decreasing = T)[1:n_sub2],]
gene_bg.w <- rownames(dat_filt2)
save(gene_bg.w, file = paste0(dir_rdata, "gene_bg.w.RData"))
rm(variance, var_p2)

# alternative filtering method --------------------------------------------

## removing all features that have a count of less than n in more than p (proportion) of the samples
# paramters for filtering
#filt_n.w & filt_p.w
#filt_n.w <- apply(dat_sub, 2, median) %>% median()
#n.smp <- ncol(dat_sub)
#del <- apply(dat_sub, 1, function(x) sum(x < filt_n.w) > floor(filt_p.w*n.smp))
#dat_filt <- dat_sub[!del,]
#dim(dat_filt); rm(del)

# 2 nomalisation ----------------------------------------------------------

# variance-stabilizing transformation

# generate DESeqDataSet
smp_m2 <- smp_m[smp_m %in% smp_r]
smp_r2 <- sample_name(colnames(dat_filt), n.field = 4)
countData <- dat_filt[,match(smp_m2, smp_r2)]
phenoDT.rna <- phenoDT.met[phenoDT.met$sample %in% smp_m2,]
# label NA as missing in specified column
for (var in var_vst) {
  phenoDT.rna[,var] <- as.character(phenoDT.rna[,var])
  phenoDT.rna[,var][is.na(phenoDT.rna[,var])] <- "missing"
}
# replace " ", "," with "_"
for (var in var_vst) {
  phenoDT.rna[,var] <- gsub(" ","_", phenoDT.rna[,var], fixed = T)
  phenoDT.rna[,var] <- gsub(",","", phenoDT.rna[,var], fixed = T)
}
dds <- DESeqDataSetFromMatrix(countData, phenoDT.rna, fml_vst)
rm(var)

# perform variance-stabilizing transformation -----------------------------

# blind=F
#system.time({vsd <- varianceStabilizingTransformation(dds, blind = F)})
system.time({vsd <- vst(dds, blind = F, nsub = 1e4)})
dat_rna.w <- assay(vsd)
# genes for WGCNA
gene.w <- rownames(dat_rna.w)

# save the data
save(phenoDT.rna, file = paste0(dir_rdata, "phenoDT.rna.RData"))
save(dat_rna.w, file = paste0(dir_rdata, "dat_rna.w.RData"))
save(gene.w, file = paste0(dir_rdata, "gene.w.RData"))


# section_2 WGCNA procedure -----------------------------------------------

# 1 binarize specified columns of clinical data ---------------------------

traitDt <- binarizeCategoricalColumns.forPlots(phenoDT.rna[,c(base_col, binar_col)], 
                                                    convertColumns = binar_col)
colnames(traitDt) <- gsub("^paper_", "", colnames(traitDt), fixed = F)

# 2 constructing a sample network for outlier detection -------------------

datExpr <- t(dat_rna.w)
# show that row names agree 
table(rownames(datExpr)==rownames(traitDt)) 

# sample network based on squared Euclidean distance 
# note that we transpose the data 
A=adjacency(t(datExpr),type="distance") 
# this calculates the whole network connectivity 
k=as.numeric(apply(A,2,sum))-1 
# standardized connectivity 
Z.k=scale(k) 

# Designate samples as outlying 
# if their Z.k value is below the threshold 
thresholdZ.k=-3.5 # often -2.5 

# the color vector indicates outlyingness (red) 
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black") 

# calculate the cluster tree using flahsClust or hclust 
sampleTree = flashClust(as.dist(1-A), method = "average") 

# Convert traits to a color representation: 
# where red indicates high values 
traitColors=data.frame(numbers2colors(traitDt,signed=FALSE)) 
dimnames(traitColors)[[2]]=paste(colnames(traitDt),"C",sep="") 
datColors=data.frame(outlierC=outlierColor,traitColors) 

# Plot the sample dendrogram and the colors underneath. 
pdf(file = paste0(dir_plot, "Sample dendrogram and trait heatmap.pdf"),
    width = 10, height = 20, pointsize = p.size_plt)
plotDendroAndColors(sampleTree,groupLabels=names(datColors), 
                    colors=datColors,main="Sample dendrogram and trait heatmap") 
dev.off()

# Remove outlying samples from expression and trait data 
remove.samples= Z.k<thresholdZ.k | is.na(Z.k) 
sum(remove.samples)

rownames(datExpr)[remove.samples]  ## confirm outlier
datExpr=datExpr[!remove.samples,] 
traitDt=traitDt[!remove.samples,] 
# Recompute the sample network among the remaining samples 
A=adjacency(t(datExpr),type="distance") 
# Let's recompute the Z.k values of outlyingness 
k=as.numeric(apply(A,2,sum))-1 
Z.k=scale(k) 

#save(datExpr, file = paste0(dir_rdata, "datExpr_temp.RData"))
#save(traitDt, file = paste0(dir_rdata, "traitDt_temp.RData"))

# 3 Choosing the soft threshold beta via scale free topology --------------

# Choose a set of soft thresholding powers 
powers=c(1:30) # in practice this should include powers up to 20. 
# choose power based on SFT criterion 
system.time({
  sft=pickSoftThreshold(datExpr,powerVector=powers, networkType = "signed")
})

#Digression: if you want to pick a soft threshold for a signed network write 
#sft=pickSoftThreshold(datExprFemale,powerVector=powers, networkType = "signed") 
# but then you should consider higher powers. Default beta=12. 

# Plot the results: 
par(mfrow=c(1,2)) 
# SFT index as a function of different powers 
pdf(file = paste0(dir_plot, "Scale independence.pdf"),
    width = 8, height = 8, pointsize = p.size_plt)
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, 
     signed R^2",type="n",main=paste("Scale independence")) 
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     labels=powers,col="red") 
# this line corresponds to using an R^2 cut-off of h 
abline(h=0.90,col="red") 
dev.off()

# Mean connectivity as a function of different powers 
pdf(file = paste0(dir_plot, "Mean connectivity.pdf"),
    width = 8, height = 8, pointsize = p.size_plt)
plot(sft$fitIndices[,1],sft$fitIndices[,5],type="n", 
     xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity")) 
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,col="red") 
dev.off()

# remove temp data
rm(smp_r, smp_r2, smp_m, smp_m2, rna.hg38, dat_rna0, dat_sub, dat_filt, dat_filt2, 
   countData, dds, vsd, A, dat_rna.w)

# to be cont'd ------------------------------------------------------------

