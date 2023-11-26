#!/usr/bin/env Rscript
# $title transform cnv data
# $input
#   cnv.hg38, p_smp
# $output
#   dat_cnv
# $author giuseppe
# $date Apr 2020

# transform tibble to data frame
dat_cnv <- as.data.frame(cnv.hg38)
# rename row as ensembl_id
rownames(dat_cnv) <- dat_cnv[,1]

# remove gene info columns
dat_cnv <- dat_cnv[,-c(1:3)]
dat_cnv <- as.matrix(dat_cnv)

# remove sample not in met data
name_cnv <- sample_name(colnames(dat_cnv),4)
name_met <- sample_name(rownames(phenoDT.met),4)
dat_cnv <- dat_cnv[,name_cnv %in% name_met]

# filter out with zero value in more than prop of sample
keep.idx <- apply(dat_cnv, 1, function(x) sum(x==0)) < ceiling(p_smp*ncol(dat_cnv))
dat_cnv.f <- dat_cnv[keep.idx,]
# transpose dat_cnv.f
dat_cnv.f <- t(dat_cnv.f)

# save the data
save(dat_cnv, file = paste0(dir_rdata, "dat_cnv.RData"))
save(dat_cnv.f, file = paste0(dir_rdata, "dat_cnv.f.RData"))
rm(cnv.hg38, name_cnv, name_met, keep.idx, p_smp)
