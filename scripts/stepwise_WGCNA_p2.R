#!/usr/bin/env Rscript
# $title stepwide WGCNA part_2
# $description
#   stepwise weighted gene correlation network analysis
#   adapted from "Corrected R code from chapter 12 of the book" of Horvath
# $input
#   rna.hg38
# $main_output
#   m_eigen, gene.w, moduleColorsManual2, phenoDT.rna
# $author giuseppe
# $date Apr 2020


# cont'd ------------------------------------------------------------------

# 4 Manual, stepwise module detection -------------------------------------

# 4.1 cut tree dynamic ----------------------------------------------------

# We now calculate the weighted adjacency matrix, using the power 13: 
# pwr <- 16
system.time({
  A.gene = adjacency(datExpr, power = pwr, type = "signed") 
})

#define a dissimilarity based on the topological overlap 
system.time({
  dissTOM = TOMdist(A.gene) 
})
#save(dissTOM, file = paste0(dir_rdata, "dissTOM_temp.RData"))

#hierarchical clustering 
geneTree = flashClust(as.dist(dissTOM),method="average") 
save(geneTree, file = paste0(dir_rdata, "geneTree.RData"))

# here we define the modules by cutting branches 
# deepSplit (range 0 to 4) provides a rough control over sensitivity to cluster splitting. 
# The higher the value (or if TRUE), the more and smaller clusters will be produced.
# deep.split <- 2
system.time({
  moduleLabelsManual1=cutreeDynamic(dendro=geneTree,distM=dissTOM, 
                                    method="hybrid",deepSplit=deep.split,
                                    pamRespectsDendro=F,minClusterSize=30)
}) 
moduleColorsManual1=labels2colors(moduleLabelsManual1)

# 4.2 plot module detection by dynamic tree cut & merged dynamic ----------

#mergingThresh <- .25
merge=mergeCloseModules(datExpr,moduleColorsManual1, 
                        cutHeight=mergingThresh) 
# resulting merged module colors 
moduleColorsManual2 = merge$colors 
save(moduleColorsManual2, file = paste0(dir_rdata, "moduleColorsManual2.RData"))

# eigengenes of the newly merged modules: 
dim(merge$newMEs) 

# Show the effect of module merging by plotting the 
# original and merged module colors below the tree 
datColors_2=data.frame(moduleColorsManual1,moduleColorsManual2) 
pdf(file = paste0(dir_plot, "hierarchical_cluster_tree_manual_unmerged_merged_",mergingThresh,".pdf"),
    width = 12, height = 8, pointsize = 20)
plotDendroAndColors(geneTree,colors=datColors_2, 
                    groupLabels=c("Dynamic tree cut","Merged dynamic"), 
                    dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05) 
dev.off()

# check the agreement between manual and automatic module labels 
mean_mod.mtch <- mean(moduleColorsManual2==moduleColorsManual1)
mean_mod.mtch

## merged dynamic column only
pdf(file = paste0(dir_plot, "hierarchical_cluster_tree_manual_merged_",mergingThresh,".pdf"),
    width = 12, height = 8, pointsize = 20)
plotDendroAndColors(geneTree,colors=moduleColorsManual2, 
                    groupLabels="Merged dynamic", 
                    dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05) 
dev.off()

# 4.3 Calculate eigengenes ------------------------------------------------

MEList=moduleEigengenes(datExpr,colors=moduleColorsManual2) 
m_eigen = MEList$eigengenes 
## check number of modules
dim(m_eigen)
save(MEList, file = paste0(dir_rdata, "MEList.RData"))
save(m_eigen, file = paste0(dir_rdata, "m_eigen.RData"))

# Plot the relationships among the eigengenes
ME_cor <- orderMEs(m_eigen)
pdf(file = paste0(dir_plot, "eigengene_cluster_tree_heatmap_merged_",mergingThresh,".pdf"),
    width = 8, height = 10, pointsize = 16)
plotEigengeneNetworks(ME_cor,"",marDendro=c(0,4,1,2), 
                      marHeatmap=c(3,4,1,2),cex.lab=0.8,xLabelsAngle=90)
dev.off()

# 5 Relating modules to clinical traits -----------------------------------

# Define numbers of genes and samples 
nGenes = ncol(datExpr) 
nSamples = nrow(datExpr) 

## subset trait data from traitDt
keep <- sapply(v_cor.heatmap, function(v) grep(v, colnames(traitDt), fixed = T))
traitDt_sub = traitDt[,unlist(keep)]
rm(keep)

# define plot function
heatmap_mod_trait_cor <- function(ME_cor = ME_cor,
                                  trait = trait, 
                                  width = 8, 
                                  height = 8, 
                                  point.size_plt = 14,
                                  directory = directory,
                                  mergingThresh = mergingThresh,
                                  name = name) {
  modTraitCor = cor(ME_cor, trait, use = "p") 
  modTraitP = corPvalueStudent(modTraitCor, nSamples) 
  #Since we have a moderately large number of modules and traits, 
  #a suitable graphical representation will help in reading 
  #the table. We color code each association by the correlation value: 
  # Will display correlations and their p-values 
  textMatrix = paste(signif(modTraitCor, 2), "\n(", 
                     signif(modTraitP, 1), ")", sep = "") 
  dim(textMatrix) = dim(modTraitCor) 
  par(mar = c(6, 8.5, 3, 3)) 
  # Display the correlation values within a heatmap plot 
  pdf(file = paste0(directory,"heatmap_module-trait_correlations_",mergingThresh,"_",name,".pdf"),
      width = width, height = height, pointsize = point.size_plt)
  labeledHeatmap(Matrix = modTraitCor, xLabels = names(trait), 
                 yLabels = names(ME_cor), ySymbols = names(ME_cor),  
                 colorLabels =FALSE,colors=greenWhiteRed(50),textMatrix=textMatrix, 
                 setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), 
                 main = paste("Module-trait relationships")) 
  dev.off()
}

# plot 1 all columns
heatmap_mod_trait_cor(ME_cor = ME_cor, trait = traitDt, width = 22, height = 12, point.size_plt = 14,
                      directory = dir_plot, mergingThresh = mergingThresh, name = "all")

# plot 2 selected column
heatmap_mod_trait_cor(ME_cor = ME_cor, trait = traitDt_sub,  width = 18, height = 12, point.size_plt = 14, 
                      directory = dir_plot, mergingThresh = mergingThresh, name = "subset")

# 6 Visualizing the network -----------------------------------------------

# 6.1 topological overlap matrix (TOM) plot -------------------------------

# Set the diagonal of the TOM disscimilarity to NA
diag(dissTOM) = NA
# Transform dissTOM with a power to enhance visibility
# pdf file too large. save in png & jpeg.
# png
png(file = paste0(dir_plot, "topological_overlap_matrix_plot_",mergingThresh,".png"),
    width = 12, height = 12, units = "in", res = 450, pointsize = 20)
TOMplot(dissim=dissTOM^7,dendro=geneTree,Colors=moduleColorsManual2, 
        main = "Network heatmap plot, all genes")
dev.off()
# jpeg
jpeg(file = paste0(dir_plot, "topological_overlap_matrix_plot_",mergingThresh,".jpeg"),
     width = 12, height = 12, units = "in", res = 450, pointsize = 20)
TOMplot(dissim=dissTOM^7,dendro=geneTree,Colors=moduleColorsManual2, 
        main = "Network heatmap plot, all genes")
dev.off()

# 6.2 Multidimensional scaling (MDS) plot ---------------------------------

cmd1=cmdscale(as.dist(dissTOM),2)
pdf(file = paste0(dir_plot, "Multidimensional_scaling_plot_",mergingThresh,".pdf"),
    width = 10, height = 10, pointsize = 20)
par(mfrow=c(1,1))
plot(cmd1,col=moduleColorsManual2, main="MDS plot", 
     xlab="Scaling Dimension 1",ylab="Scaling Dimension 2")
dev.off()

# section_3 functional annotation of gene modules -------------------------

# set directory
dir_fun_annot <- "~/Project_Data/R_data/script_function/gene_set_annot/"
dir_go <- paste0(dir_plot, "enrichment_analysis/GO/")
dir_kegg <- paste0(dir_plot, "enrichment_analysis/KEGG/")
dir.create(dir_go, recursive = T)
dir.create(dir_kegg, recursive = T)

# load function
source(file = paste0(dir_fun_annot, "gene_set_annot.R"))

# enrich GO
system.time({res_eGO <- enrichGO_multi(gene_list = gene.w, gene_label = moduleColorsManual2, 
                                       drop_label = "grey",keyType = "ENSEMBL", 
                                       n_showCategory = 8, width = 10, height = 8, directory = dir_go)})

# enrich KEGG
system.time({res_eKEGG <- enrichKEGG_multi(gene_list = gene.w, gene_label = moduleColorsManual2, 
                                      drop_label = "grey", fromType = "ENSEMBL", 
                                      n_showCategory = 8, directory = dir_kegg)})

# remove temp data
rm(datExpr, A.gene, dissTOM, res_eGO, res_eKEGG, gene_bg.w)

