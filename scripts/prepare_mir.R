#!/usr/bin/env Rscript
# $title mirna data preparation
# $description
#   subset reads_per_million data
# $input
#   mir.hg38, min.var
# $output
#   dat_mir
# $author giuseppe
# $date Apr 2020

# subset reads_per_million data
dat_mir <- mir.hg38
row_name <- dat_mir[,1]
dat_mir <- mir.hg38[,grep("^reads_per_million", colnames(mir.hg38), fixed = F)]
colnames(dat_mir) <- gsub("^reads_per_million_miRNA_mapped_", "", colnames(dat_mir), fixed = F)
rownames(dat_mir) <- row_name

# remove sample not in met data
name_mir <- sample_name(colnames(dat_mir),4)
name_met <- sample_name(rownames(phenoDT.met),4)
dat_mir <- dat_mir[,name_mir %in% name_met]

# filter out genes whose variance across obs is smaller than min.var
#min.var=.1
variance <- apply(dat_mir, 1, var)
dat_mir.f <- dat_mir[!variance < min.var,]
dat_mir.f <- t(dat_mir.f)
dat_mir.f <- as.matrix(dat_mir.f)

save(dat_mir, file = paste0(dir_rdata, "dat_mir.RData"))
save(dat_mir.f, file = paste0(dir_rdata, "dat_mir.f.RData"))
rm(mir.hg38, row_name, variance, min.var, name_mir, name_met)


