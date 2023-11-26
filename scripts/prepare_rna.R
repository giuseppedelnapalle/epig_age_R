#!/usr/bin/env Rscript
# $title RNA data preparation
# $description
#   filter & normalise RNA profiles
# $input
#   rna.hg38
# $output
#   dat_rna, phenoDT.rna
# $author giuseppe
# $date Apr 2020

library(DESeq2)

# 1 gene filtering --------------------------------------------------------

# removing all features that have a count of less than n in more than p (proportion) of the samples
dat_rna0 <- assay(rna.hg38)
dim(dat_rna0)
# obtain data of tumour sample whose DNAm data available
smp_r <- sample_name(colnames(dat_rna0), n.field = 4)
smp_m <- sample_name(phenoDT.met$sample, n.field = 4)
dat_sub <- dat_rna0[,smp_r %in% smp_m]
dim(dat_sub)

# paramters for filtering
#filt_n & filt_p
filt_n <- apply(dat_sub, 2, median) %>% median()
n.smp <- ncol(dat_sub)
del <- apply(dat_sub, 1, function(x) sum(x < filt_n) > floor(filt_p*n.smp))
dat_filt <- dat_sub[!del,]
dim(dat_filt)
rm(del)

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
system.time({
  vsd <- varianceStabilizingTransformation(dds, blind = F)
})
dat_rna <- assay(vsd)

# blind=T, for QC
system.time({
  vsd2 <- varianceStabilizingTransformation(dds, blind = T)
})
# plot PCA
pdf(file = paste0(dir_plot, "PCA_rna_dt_tumour.pdf"),
    width = 8, height = 8, pointsize = p.size_plt)
plotPCA(vsd2, intgroup=group_PCA)
dev.off()

# save the data
save(dds, file = paste0(dir_rdata, "dds_temp.RData"))
save(vsd, file = paste0(dir_rdata, "vsd_temp.RData"))
save(vsd2, file = paste0(dir_rdata, "vsd2_blind_temp.RData"))
save(phenoDT.rna, file = paste0(dir_rdata, "phenoDT.rna.RData"))
save(dat_rna, file = paste0(dir_rdata, "dat_rna.RData"))

# remove temp data
rm(smp_r, smp_r2, smp_m, smp_m2, dat_rna0, dat_sub, dat_filt, countData, dds, vsd, vsd2)

