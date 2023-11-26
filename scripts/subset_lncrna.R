#!/usr/bin/env Rscript
# $title subset lncRNA data
# $input
#   dat_rna
# $output
#   dat_lnc
# $author giuseppe
# $date Apr 2020


# set directory of gene information
dir_gene <- "~/Project_Data/R_data/cancer_research/suppl_TCGA/lncrna_biomart/"

# load ensbl_lnc
load(paste0(dir_gene, "ensbl_lnc.RData"))

# subset lncRNA data
dat_lnc <- dat_rna[rownames(dat_rna) %in% ensbl_lnc,]

# remove sample not in met data
name_lnc <- sample_name(colnames(dat_lnc),4)
name_met <- sample_name(rownames(phenoDT.met),4)
dat_lnc <- dat_lnc[,name_lnc %in% name_met]

# filter out low variance genes
#p_var <- .5
variance <- apply(dat_lnc, 1, var)
keep.idx <- order(variance, decreasing = T)[1:(p_var*nrow(dat_lnc))]
dat_lnc.f <- dat_lnc[keep.idx,]
dat_lnc.f <- t(dat_lnc.f)
dat_lnc.f <- as.matrix(dat_lnc.f)

save(dat_lnc, file = paste0(dir_rdata, "dat_lnc.RData"))
save(dat_lnc.f, file = paste0(dir_rdata, "dat_lnc.f.RData"))
rm(dat_rna, p_var, variance, keep.idx, name_lnc, name_met)

