#!/usr/bin/env Rscript
# $title somatic mutation data preparation
# $description
#   subset & arrange maf data for DNA methy. analysis
# $input
#   maf.hg38, maf.hg38_2, prep_maf2, phenoDT.met, p_smp
# $output
#   dat_mut, dat_mut2, dat_mut.f, dat_mut2.f, maf_sub, maf_sub2
# $author giuseppe
# $date Apr 2020

#prep_maf2 <- T
#p_smp <- 1/40
#p_smp2 <- 1/12
#min_count <- 5

library(reshape2)

# retrieve required columns of maf.hg38
col_mut <- c("Hugo_Symbol", "Entrez_Gene_Id", "Start_Position", "End_Position",
             "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", 
             "Tumor_Seq_Allele2", "dbSNP_RS", "Tumor_Sample_Barcode", 
             "Matched_Norm_Sample_Barcode", "Mutation_Status", "HGVSc", "HGVSp", 
             "HGVSp_Short", "Transcript_ID", "t_depth", "t_ref_count", "t_alt_count", 
             "n_depth", "all_effects", "Allele", "Gene", "Feature", "Feature_type", 
             "One_Consequence","Consequence","Amino_acids", "BIOTYPE", "CANONICAL", 
             "ENSP", "SWISSPROT", "RefSeq", "SIFT", "PolyPhen", "EXON", "INTRON", 
             "DOMAINS", "IMPACT", "COSMIC")

if (prep_maf2) {
  # mutect2
  maf_sub <- maf.hg38[,col_mut]
  maf_sub <- transform(maf_sub, 
                       mutation_ID = paste(maf_sub$Hugo_Symbol, maf_sub$Start_Position, 
                                           maf_sub$End_Position,maf_sub$Variant_Classification, sep = "."))
  save(maf_sub, file = paste0(dir_rdata, "maf_sub_mutect2.RData"))
  maf_shrt <- maf_sub[,c("Tumor_Sample_Barcode", "Hugo_Symbol", "mutation_ID")]
  
  # varscan2
  maf_sub2 <- maf.hg38_2[,col_mut]
  maf_sub2 <- transform(maf_sub2, 
                       mutation_ID = paste(maf_sub2$Hugo_Symbol, maf_sub2$Start_Position, 
                                           maf_sub2$End_Position,maf_sub2$Variant_Classification, sep = "."))
  save(maf_sub2, file = paste0(dir_rdata, "maf_sub2_varscan2.RData"))
  maf_shrt2 <- maf_sub2[,c("Tumor_Sample_Barcode", "Hugo_Symbol", "mutation_ID")]
} else {
  maf_sub <- maf.hg38[,col_mut]
  maf_sub <- transform(maf_sub, 
                       mutation_ID = paste(maf_sub$Hugo_Symbol, maf_sub$Start_Position, 
                                           maf_sub$End_Position,maf_sub$Variant_Classification, sep = "."))
  save(maf_sub, file = paste0(dir_rdata, "maf_sub_mutect2.RData"))
  maf_shrt <- maf_sub[,c("Tumor_Sample_Barcode", "Hugo_Symbol", "mutation_ID")]
}

# rearrange data frame & filter out rare mutation
if (prep_maf2) {
  # mutect2
  dat_mut <- melt(maf_shrt, id = c("Tumor_Sample_Barcode", "Hugo_Symbol"))
  dat_mut <- dat_mut[,-3]
  colnames(dat_mut)[3] <- "variable"
  dat_mut <- transform(dat_mut, value = rep(1, nrow(dat_mut)))
  dat_mut <- dcast(dat_mut, Tumor_Sample_Barcode~Hugo_Symbol, sum)
  rownames(dat_mut) <- dat_mut[,1]
  dat_mut <- dat_mut[,-1]
  # remove sample  not in met data
  name_mut <- sample_name(rownames(dat_mut),4)
  name_met <- sample_name(rownames(phenoDT.met),4)
  dat_mut <- dat_mut[name_mut %in% name_met,]
  # update name
  name_mut2 <- sample_name(rownames(dat_mut),4)
  name_met2 <- name_met[name_met %in% name_mut2]
  dat_mut <- dat_mut[match(name_met2, name_mut2),]
  save(dat_mut, file = paste0(dir_rdata, "dat_mut_mutect2.RData"))
  # filter out rare mutation
  keep_idx <- apply(dat_mut, 2, sum) > max(min_count, round(nrow(dat_mut)*p_smp))
  dat_mut.f <- dat_mut[,keep_idx]
  dat_mut.f <- as.matrix(dat_mut.f)
  save(dat_mut.f, file = paste0(dir_rdata, "dat_mut.f_mutect2.RData"))
  # filter 2
  keep_idx2 <- apply(dat_mut, 2, sum) > max(min_count, round(nrow(dat_mut)*p_smp2))
  dat_mut.f.s <- dat_mut[,keep_idx2]
  dat_mut.f.s <- as.matrix(dat_mut.f.s)
  save(dat_mut.f.s, file = paste0(dir_rdata, "dat_mut.f.s_mutect2.RData"))
  rm(name_mut, name_met, name_met2, name_mut2)
  
  # varscan2
  dat_mut2 <- melt(maf_shrt2, id = c("Tumor_Sample_Barcode", "Hugo_Symbol"))
  dat_mut2 <- dat_mut2[,-3]
  colnames(dat_mut2)[3] <- "variable"
  dat_mut2 <- transform(dat_mut2, value = rep(1, nrow(dat_mut2)))
  dat_mut2 <- dcast(dat_mut2, Tumor_Sample_Barcode~Hugo_Symbol, sum)
  rownames(dat_mut2) <- dat_mut2[,1]
  dat_mut2 <- dat_mut2[,-1]
  # remove sample  not in met data
  name_mut <- sample_name(rownames(dat_mut2),4)
  name_met <- sample_name(rownames(phenoDT.met),4)
  dat_mut2 <- dat_mut2[name_mut %in% name_met,]
  # update name
  name_mut2 <- sample_name(rownames(dat_mut2),4)
  name_met2 <- name_met[name_met %in% name_mut2]
  dat_mut2 <- dat_mut2[match(name_met2, name_mut2),]
  save(dat_mut2, file = paste0(dir_rdata, "dat_mut2_varscan2.RData"))
  # filter out rare mutation
  keep_idx <- apply(dat_mut2, 2, sum) > max(min_count, round(nrow(dat_mut2)*p_smp))
  dat_mut2.f <- dat_mut2[,keep_idx]
  dat_mut2.f <- as.matrix(dat_mut2.f)
  save(dat_mut2.f, file = paste0(dir_rdata, "dat_mut2.f_varscan2.RData"))
  # filter 2
  keep_idx2 <- apply(dat_mut2, 2, sum) > max(min_count, round(nrow(dat_mut2)*p_smp2))
  dat_mut2.f.s <- dat_mut2[,keep_idx2]
  dat_mut2.f.s <- as.matrix(dat_mut2.f.s)
  save(dat_mut2.f.s, file = paste0(dir_rdata, "dat_mut2.f.s_varscan2.RData"))
  rm(name_mut, name_met, name_met2, name_mut2)
} else {
  # mutect2
  dat_mut <- melt(maf_shrt, id = c("Tumor_Sample_Barcode", "Hugo_Symbol"))
  dat_mut <- dat_mut[,-3]
  colnames(dat_mut)[3] <- "variable"
  dat_mut <- transform(dat_mut, value = rep(1, nrow(dat_mut)))
  dat_mut <- dcast(dat_mut, Tumor_Sample_Barcode~Hugo_Symbol, sum)
  rownames(dat_mut) <- dat_mut[,1]
  dat_mut <- dat_mut[,-1]
  # remove sample not in met data
  name_mut <- sample_name(rownames(dat_mut),4)
  name_met <- sample_name(rownames(phenoDT.met),4)
  dat_mut <- dat_mut[name_mut %in% name_met,]
  # update name
  name_mut2 <- sample_name(rownames(dat_mut),4)
  name_met2 <- name_met[name_met %in% name_mut2]
  dat_mut <- dat_mut[match(name_met2, name_mut2),]
  save(dat_mut, file = paste0(dir_rdata, "dat_mut_mutect2.RData"))
  # filter out rare mutation
  keep_idx <- apply(dat_mut, 2, sum) > max(min_count, round(nrow(dat_mut)*p_smp))
  dat_mut.f <- dat_mut[,keep_idx]
  dat_mut.f <- as.matrix(dat_mut.f)
  save(dat_mut.f, file = paste0(dir_rdata, "dat_mut.f_mutect2.RData"))
  # filter 2
  keep_idx2 <- apply(dat_mut, 2, sum) > max(min_count, round(nrow(dat_mut)*p_smp2))
  dat_mut.f.s <- dat_mut[,keep_idx2]
  dat_mut.f.s <- as.matrix(dat_mut.f.s)
  save(dat_mut.f.s, file = paste0(dir_rdata, "dat_mut.f.s_mutect2.RData"))
  rm(name_mut, name_met, name_met2, name_mut2)
}

rm(maf.hg38, maf.hg38_2, maf_sub, maf_sub2, maf_shrt ,maf_shrt2, p_smp, keep_idx, keep_idx2)

