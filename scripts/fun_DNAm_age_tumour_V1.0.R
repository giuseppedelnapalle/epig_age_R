#!/usr/bin/env Rscript
# $title functions for epigenetic age analysis for tumour projects
# $author giuseppe
# $date Apr 2020

library(MASS)
library(multtest)

# SECTION_1 DNA methy. model ----------------------------------------------

# 1 Horvath model ---------------------------------------------------------

# define age transformation functions
trafo= function(x,adult.age=20){ 
  x=(x+1)/(1+adult.age); y=ifelse(x<=1, log(x),x-1);y 
}

anti.trafo= function(x,adult.age=20){ 
  ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) 
}

# define sample filtering function
# note rows are CpGs, colums are samples
norm_filter <- function(x){
  # calculate mean inter-array correlation
  cor.mat = cor(x, use = "complete.obs")
  mean_cor = apply(cor.mat, 1, mean)
  
  # calculate max beta for each sample
  max_b = apply(x, 2, max, na.rm = T)
  
  # filter sample by mean inter-array correlation & max beta
  x = x[, mean_cor >= .9 & max_b >= .96]
  return(x)
}

# median absolute deviation
medianAbsDev <- function(x,y) median(abs(x-y),na.rm=T)

# 2 epiTOC model ----------------------------------------------------------

# probe filtering function
filter_probe <- function(data = data, 
                         n_sample = n_sample, 
                         keep_probe = keep_probe){
  keep = apply(as.matrix(is.na(data[,-1])),1,sum) < n_sample | data[,1] %in% keep_probe
  data = data[keep,]
}

# 3 general function ------------------------------------------------------

# Age Acceleration defined as estimated_age - age
ageAcc.Diff <- function(datout = datout, 
                        phenoDT = phenoDT, 
                        epiAge = epiAge, 
                        age = age, 
                        is.select = is.select,
                        sample.name = sample.name){
  restNonMissing = !is.na(datout[,epiAge]) & !is.na(phenoDT[,age])
  if (sum(restNonMissing,na.rm=T) >3) {
    age.Acc = rep(NA, nrow(phenoDT))
    aae.Acc = as.numeric(age.Acc)
    age.Acc[restNonMissing & is.select] = datout[,epiAge][restNonMissing & is.select] - phenoDT[,age][restNonMissing & is.select]
    names(age.Acc) = gsub("^X", "", phenoDT[,sample.name], fixed = F)
    return(age.Acc)
  } else print("# of restNonMissing <= 3")
} # end of function


# SECTION_2 plot ----------------------------------------------------------

# scatter plot function 1
scatter_plt_1 <- function(dat_DNAm = dat_DNAm,
                          pheno_dt = pheno_dt,
                          is.select = is.select,
                          factor.v= factor.v,
                          epiAge = epiAge,
                          age = "age",
                          model = model,
                          shape = "tumour_Yes_No",
                          point.size = 2,
                          height = 10, width = 16,
                          ref_line = T,
                          fit_line = F,
                          directory = directory){
  restNonMissing = !is.na(dat_DNAm[,epiAge]) & !is.na(pheno_dt[,age])
  data0 = cbind(pheno_dt, dat_DNAm)
  data_plt0 = data0[restNonMissing,]
  for (fct in factor.v) {
    nm_fct = gsub("^paper_","", fct, fixed = F)
    # remove sample with NA for specified factor
    del = is.na(data_plt0[,fct]) | data_plt0[,fct] == "NA"
    data_plt = data_plt0[!del,]
    plot = ggplot(data_plt, aes_string(x = epiAge, y = age, shape = shape, color = fct)) +
      geom_point(size = point.size) + 
      labs(title=paste0("Scatterplot of Chronological Age vs ",epiAge), 
                                           x=epiAge, y="Chronological Age (years)")
    if (ref_line) {
      if (fit_line) {
        # fit robust linear model on normal tissue data
        reg = rlm(as.formula(paste(age,"~",epiAge,sep = " ")), data = data0[!is.select,], na.action = "na.omit")
        coeff = coefficients(reg)
        plot = plot + 
          geom_abline(intercept = round(coeff[1],1), slope = round(coeff[2],1), color="red", size = .8) + 
          geom_abline(intercept = 0, slope = 1, color= "black", size = .6, linetype="dashed")
        ggsave(filename = paste0(directory, "Scatterplot_",epiAge,"_Age_",nm_fct,"_",model,".pdf"), 
               plot = plot, height = height, width = width, units = "cm")
      } else {
        plot = plot + 
          geom_abline(intercept = 0, slope = 1, color= "black", size = .6, linetype="dashed")
        ggsave(filename = paste0(directory, "Scatterplot_",epiAge,"_Age_",nm_fct,"_",model,".pdf"), 
               plot = plot, height = height, width = width, units = "cm")
      } # end of else
    } else {
      if (fit_line) {
        reg = rlm(as.formula(paste(age,"~",epiAge,sep = " ")), data = data0[!is.select,], na.action = "na.omit")
        coeff = coefficients(reg)
        plot = plot + 
          geom_abline(intercept = round(coeff[1],1), slope = round(coeff[2],1), color="red", size = .8, linetype="solid")
        ggsave(filename = paste0(directory, "Scatterplot_",epiAge,"_Age_",nm_fct,"_",model,".pdf"), 
               plot = plot, height = height, width = width, units = "cm")
      } else {
        ggsave(filename = paste0(directory, "Scatterplot_",epiAge,"_Age_",nm_fct,"_",model,".pdf"), 
               plot = plot, height = height, width = width, units = "cm")
      } # end of else
    } # end of else
  } # end of for
} # end of function

# scatter plot function 2
scatter_plt_2 <- function(dat_DNAm = dat_DNAm,
                          pheno_dt = pheno_dt,
                          is.select = is.select,
                          epiAge = epiAge,
                          age = "age",
                          height = 8, width = 8, point.size_plt = 14,
                          model = model,
                          directory = directory,
                          name = name){
  restNonMissing = !is.na(dat_DNAm[,epiAge]) & !is.na(pheno_dt[,age])
  medianAbsDev1=signif(medianAbsDev(dat_DNAm[,epiAge][restNonMissing][is.select], 
                                    pheno_dt[,age][restNonMissing][is.select]),2)
  pdf(file = paste0(directory, "Scatterplot_",epiAge,"_Age_",model,"_",name,".pdf"),
      width = width, height = height, pointsize = point.size_plt)
  par(mfrow=c(1,1))
  verboseScatterplot(dat_DNAm[,epiAge][restNonMissing & is.select], 
                     pheno_dt[,age][restNonMissing & is.select],
                     xlab=epiAge, ylab="Chronological Age",
                     main=paste("All, err=",medianAbsDev1))
  abline(0,1)
  dev.off()
}

# scatter plot function 3
scatter_plt_3 <- function(X = X,
                          y = y,
                          ylab = ylab,
                          factor.v= factor.v,
                          point.size = 2,
                          height = 10, width = 12,
                          ref_line = F,
                          fit_line = T,
                          fit_type = "rlm",  # c("rlm", "lm")
                          cor_method = "pearson",
                          alternative = "two.sided",
                          theme_ggplot = theme_ggplot,
                          directory = directory,
                          name_plt = name_plt,
                          model = model){
  # match obs of X & y
  rownames(X) = sample_name(rownames(X), 4)
  names(y) = sample_name(names(y), 4)
  X = X[rownames(X) %in% names(y),]
  y = y[names(y) %in% rownames(X)]
  X = X[match(names(y), rownames(X)),]
  # remove obs with missing response value
  del.idx = is.na(y)
  y = y[!del.idx]
  X = X[!del.idx,]
  # rename column of X
  colnames(X) = gsub("-", "_", colnames(X), fixed = T)
  # combine X & y
  data0 = cbind(as.data.frame(y), X)
  colnames(data0)[1] = ylab
  # rename factor.v
  factor.v = gsub("-", "_", factor.v, fixed = T)
  for (fct in factor.v) {
    cor = cor(X[,fct], y, method = cor_method)
    cor_test = cor.test(X[,fct], y, alternative = alternative, method = cor_method)
    nm_fct = gsub("^paper_","", fct, fixed = F)
    # remove sample with NA in specified factor & response variable
    del = is.na(data0[,fct]) | data0[,fct] == "NA" | is.na(data0[,ylab])
    data = data0[!del,]
    
    plot = ggplot(data, aes_string(x = fct, y = ylab)) +
      geom_point(size = point.size) + 
      labs(title=paste0(nm_fct,", ","cor=", round(cor, digits = 3), ", ", "p=", round(unlist(cor_test[3]),digits = 3)), 
           x=nm_fct, y=ylab) + 
      theme_ggplot
    if (ref_line) {
      if (fit_line) {
        # fit linear model on normal tissue data
        if (fit_type=="lm") reg = lm(as.formula(paste(ylab,"~",fct,sep = " ")), data = data, na.action = "na.omit") else
          reg = rlm(as.formula(paste(ylab,"~",fct,sep = " ")), data = data, na.action = "na.omit")
        coeff = coefficients(reg)
        plot = plot + 
          geom_abline(intercept = round(coeff[1],1), slope = round(coeff[2],1), color="red", size = .8, linetype="solid") + 
          geom_abline(intercept = 0, slope = 1, color= "black", size = .6, linetype="dashed")
        ggsave(filename = paste0(directory, "Scatterplot_",ylab,"_",nm_fct,"_",name_plt,"_",model,".pdf"), 
               plot = plot, height = height, width = width, units = "cm")
      } else {
        plot = plot + 
          geom_abline(intercept = 0, slope = 1, color= "black", size = .6, linetype="dashed")
        ggsave(filename = paste0(directory, "Scatterplot_",ylab,"_",nm_fct,"_",name_plt,"_",model,".pdf"), 
               plot = plot, height = height, width = width, units = "cm")
      } # end of else
    } else {
      if (fit_line) {
        if (fit_type=="lm") reg = lm(as.formula(paste(ylab,"~",fct,sep = " ")), data = data, na.action = "na.omit") else
          reg = rlm(as.formula(paste(ylab,"~",fct,sep = " ")), data = data, na.action = "na.omit")
        coeff = coefficients(reg)
        plot = plot + 
          geom_abline(intercept = round(coeff[1],1), slope = round(coeff[2],1), color="red", size = .8, linetype="solid")
        ggsave(filename = paste0(directory, "Scatterplot_",ylab,"_",nm_fct,"_",name_plt,"_",model,".pdf"), 
               plot = plot, height = height, width = width, units = "cm")
      } else {
        ggsave(filename = paste0(directory, "Scatterplot_",ylab,"_",nm_fct,"_",name_plt,"_",model,".pdf"), 
               plot = plot, height = height, width = width, units = "cm")
      } # end of else
    } # end of else
  } # end of for
} # end of function

# box & scatter plot function 1
box_scatter_plt_1 <- function(age.Acc = age.Acc, 
                              pheno_dt = pheno_dt, 
                              age = "age", 
                              factor.v = factor.v,
                              ylab = ylab, 
                              model = model, 
                              point.size = 2,
                              height = 10, 
                              width = 16,
                              directory = directory,
                              theme_ggplot = theme_ggplot,
                              name_epiAge = name_epiAge){
  restNonMissing = !is.na(age.Acc) & !is.na(pheno_dt[,age])
  data_plt0 = cbind(pheno_dt, as.data.frame(age.Acc))[restNonMissing,]
  for (fct in factor.v) {
    nm_fct = gsub("^paper_","", fct, fixed = F)
    # remove sample with NA for specified factor & epig.age
    del = is.na(data_plt0[,fct]) | data_plt0[,fct] == "NA"
    data_plt = data_plt0[!del,]
    plot = ggplot(data_plt, aes_string(x = fct, y = "age.Acc", color = fct, fill = fct)) +
      # argument "outlier.shape = NA" remove outlier in boxplot
      geom_boxplot(color = "black", alpha = .3, outlier.shape = NA) +
      geom_point(position = "jitter", size = point.size, alpha = .8) +
      labs(title=paste0(ylab ," among ",nm_fct," groups"), x=nm_fct, y=ylab) +
      theme_ggplot
    # wide: width = 20, (point) size = 2.5
    ggsave(filename = paste0(directory,"Box_plot_",name_epiAge,"_",nm_fct,"_",model,"_width",width,".pdf"), 
           plot = plot, height = height, width = width, units = "cm")
  }
}

# box & scatter plot function 2
box_scatter_plt_2 <- function(age.Acc = age.Acc, 
                              pheno_dt = pheno_dt, 
                              factor.v = factor.v,
                              ylab = ylab, 
                              model = model, 
                              point.size = 2,
                              height = 8, 
                              width = 8,
                              directory = directory,
                              theme_ggplot = theme_ggplot,
                              name_epiAge = name_epiAge,
                              name_plt = name_plt){
  # match sample
  names(age.Acc) = sample_name(names(age.Acc), 4)
  rownames(pheno_dt) = sample_name(rownames(pheno_dt), 4)
  age.Acc = age.Acc[names(age.Acc) %in% rownames(pheno_dt)]
  pheno_dt = pheno_dt[rownames(pheno_dt) %in% names(age.Acc),]
  pheno_dt = pheno_dt[match(names(age.Acc), rownames(pheno_dt)),]
  # remove obs with NA age.Acc & combine data
  restNonMissing = !is.na(age.Acc)
  data_plt0 = cbind(pheno_dt, as.data.frame(age.Acc))[restNonMissing,]
  for (fct in factor.v) {
    nm_fct = gsub("^paper_","", fct, fixed = F)
    # remove sample with NA for specified factor & epig.age
    del = is.na(data_plt0[,fct]) | data_plt0[,fct] == "NA"
    data_plt = data_plt0[!del,]
    plot = ggplot(data_plt, aes_string(x = fct, y = "age.Acc", color = fct, fill = fct)) +
      # argument "outlier.shape = NA" remove outlier in boxplot
      geom_boxplot(color = "black", alpha = .3, outlier.shape = NA) +
      geom_point(position = "jitter", size = point.size, alpha = .8) +
      labs(title=paste0(ylab ," among ",nm_fct," groups"), x=nm_fct, y=ylab) +
      theme_ggplot
    # wide: width = 20, (point) size = 2.5
    ggsave(filename = paste0(directory,"Box_plot_",name_epiAge,"_",nm_fct,"_",model,"_width",width,"_",name_plt,".pdf"), 
           plot = plot, height = height, width = width, units = "cm")
  }
}

# ROC curve
# ROCR package
# from An Introduction to Statistical Learning, 2017
rocplot=function(pred = pred, 
                 truth = truth, 
                 age.Acc = age.Acc,
                 model = model, 
                 height = 8, 
                 width = 8,
                 point.size_plt = 14,
                 directory = directory, ...){
  restNonMissing = !is.na(pred) & !is.na(truth)
  predob = prediction(pred[restNonMissing], truth[restNonMissing])
  perf = performance(predob, "tpr", "fpr")
  pdf(file = paste0(directory, "ROC_curve_",age.Acc,"_",model,".pdf"),
      width = width, height = height, pointsize = point.size_plt)
  plot(perf,...)
  dev.off()
  }

# PRROC package
plot_PRROC <- function(scores.class0 = scores.class0,
                       class.label = class.label,
                       class.ref = class.ref,
                       curve = T,
                       auc.main = T,
                       legend = ifelse(is.logical(color) & color==T,4,NA),
                       main = NULL,
                       color = F,
                       add = F,
                       age.Acc = age.Acc,
                       model = model, 
                       height = 8, 
                       width = 8,
                       point.size_plt = 14,
                       directory = directory) {
  restNonMissing = !is.na(scores.class0) & !is.na(class.label)
  scores.class0 = scores.class0[restNonMissing]
  weights.class0 = ifelse(class.label == class.ref, 0, 1)[restNonMissing]
  PRROC_obj = roc.curve(scores.class0 = scores.class0, weights.class0 = weights.class0,
                         curve = curve)
  pdf(file = paste0(directory, "ROC_curve_",age.Acc,"_",model,".pdf"),
      width = width, height = height, pointsize = point.size_plt)
  plot(PRROC_obj, auc.main = auc.main, legend = legend, main = main, color = color, add = add)
  dev.off()
}


# SECTION_3 statistical analysis ------------------------------------------

# define F-P permutation test function for k sample (k>2)
perm_test_mult <- function(data_perm = data_perm, 
                           factor_perm = factor_perm, 
                           response_var = response_var, 
                           n_resample = n_resample){
  tmp_dat0 = data_perm
  # remove observation with NA ageAcc.diff
  tmp_dat0 = tmp_dat0[!is.na(tmp_dat0[,response_var]),]
  
  # loop over factor
  Res_FP_perm = lapply(factor_perm, function(factor){
    fmla = as.formula(paste0(response_var, " ~ ",factor))
    tmp_dat = tmp_dat0
    # remove NA observation
    tmp_dat = tmp_dat[!is.na(tmp_dat[,factor]),]
    tmp_dat = tmp_dat[!tmp_dat[,factor] == "NA",]
    #In the coin package, categorical variables and ordinal variables must
    #be coded as factors and ordered factors, respectively
    tmp_dat[,factor] = as.factor(tmp_dat[,factor])
    Res = oneway_test(fmla, data = tmp_dat, 
                                   distribution = approximate(nresample=n_resample))
    return(Res)
  })
  names(Res_FP_perm) = factor_perm
  return(Res_FP_perm)
}

# define F-P permutation test function for two sample
perm_test_bi <- function(data_perm = data_perm, 
                         factor_perm = factor_perm, 
                         response_var = response_var, 
                         n_resample = 1e5, 
                         name_normal = "Normal"){
  tmp_dat0 = data_perm
  # remove observation with NA ageAcc.diff
  tmp_dat0 = tmp_dat0[!is.na(tmp_dat0[,response_var]),]
  # remove factor with fewer than 3 levels
  num.level = sapply(factor_perm, function(x) length(levels(as.factor(tmp_dat0[,x]))))
  factor_perm = factor_perm[num.level > 2]
  
  # loop over factor
  Res_FP_perm = lapply(factor_perm, function(factor){
    fmla = as.formula(paste0(response_var, " ~ ",factor))
    level.fct = levels(as.factor(tmp_dat0[,factor]))
    # skip level normal & NA
    del1 = grep(name_normal, level.fct, fixed = T)
    del2 = grep("NA", level.fct, fixed = T)
    level.fct = level.fct[-c(del1, del2)]
    
    # loop over level
    Res.lvl = lapply(level.fct, function(level){
      # keep only sample of level fct or normal 
      tmp_dat = tmp_dat0
      tmp_dat = tmp_dat[tmp_dat[,factor] %in% c(name_normal, level),]
      # remove NA observation
      tmp_dat = tmp_dat[!is.na(tmp_dat[,factor]),]
      #tmp_dat = tmp_dat[!tmp_dat[,factor] == "NA",]
      #In the coin package, categorical variables and ordinal variables must be coded as factors and ordered factors, respectively
      tmp_dat[,factor] = as.factor(tmp_dat[,factor])
      res_tmp = oneway_test(fmla, data = tmp_dat, 
                            distribution = approximate(nresample=n_resample))
      pval_tmp = pvalue(res_tmp)
      res = list(res_tmp, pval_tmp, level)
      return(res) # result of level
    }) # end of function of level)
    return(Res.lvl) # result of factor
  }) # end of function of factor
  names(Res_FP_perm) = factor_perm
  return(Res_FP_perm)
} # end of function

# define F-P permutation test function for two sample with two level factor
perm_test_bi_2 <- function(data_perm = data_perm, 
                         factor_perm = factor_perm, 
                         response_var = response_var, 
                         n_resample = 1e5){
  tmp_dat0 = data_perm
  # remove observation with NA ageAcc.diff
  tmp_dat0 = tmp_dat0[!is.na(tmp_dat0[,response_var]),]
  # loop over factor
  Res_FP_perm = lapply(factor_perm, function(factor){
    fmla = as.formula(paste0(response_var, " ~ ",factor))
    # remove NA observation
    tmp_dat = tmp_dat0[!is.na(tmp_dat0[,factor]),]
    #tmp_dat = tmp_dat[!tmp_dat[,factor] == "NA",]
    #In the coin package, categorical variables and ordinal variables must be coded as factors and ordered factors, respectively
    tmp_dat[,factor] = as.factor(tmp_dat[,factor])
    res_tmp = oneway_test(fmla, data = tmp_dat, 
                          distribution = approximate(nresample=n_resample))
    pval_tmp = pvalue(res_tmp)
    res = list(res_tmp, pval_tmp)
    return(res) # result of factor
  }) # end of function of factor
  names(Res_FP_perm) = factor_perm
  return(Res_FP_perm)
} # end of function

# multtest package
# adjust p-value
adjust_p <- function(Res_perm_bi = Res_perm_bi, 
                    proc = "BH", 
                    alpha = .05){
  res_adj.p = list()
  pval.l = list()
  lvl.l = list()
  for (i in 1:length(Res_perm_bi)) {
    pval.l[[i]] = list()
    lvl.l[[i]] = list()
    for (j in 1:length(Res_perm_bi[[i]])) {
      pval.l[[i]] = list(pval.l[[i]], Res_perm_bi[[i]][[j]][2])
      lvl.l[[i]] = list(lvl.l[[i]], Res_perm_bi[[i]][[j]][3])
    } # end of for of j
  pval_tmp = unlist(pval.l[[i]])
  lvl_tmp = unlist(lvl.l[[i]])
  res_adj.p[[i]] = mt.rawp2adjp(pval_tmp, proc = proc, alpha = alpha, na.rm = T)
  # to_be_modified
  res_adj.p[[i]][1]$adjp = as.data.frame(res_adj.p[[i]][1]$adjp)
  res_adj.p[[i]][1]$adjp = transform(res_adj.p[[i]][1]$adjp, 
                                       contrast = lvl_tmp[match(res_adj.p[[i]][1]$adjp[,1], pval_tmp)])
  } # end of for of i
  names(res_adj.p) = names(Res_perm_bi)
  return(res_adj.p)
}

# adjust p-value for results of perm_test_bi_2
adjust_p_2 <- function(Res_perm_bi = Res_perm_bi, 
                     proc = "BH", 
                     alpha = .05){
  pval.l = list()
  for (i in 1:length(Res_perm_bi)) {
    pval.l[[i]] = as.numeric(Res_perm_bi[[i]][[2]])
  } # end of for
  pval_tmp = unlist(pval.l)
  res_adj.p = mt.rawp2adjp(pval_tmp, proc = "BH", alpha = .05, na.rm = T)
  # to_be_modified
  res_adj.p[[1]] = as.data.frame(res_adj.p[[1]])
  res_adj.p[[1]] = transform(res_adj.p[[1]], 
                             factor = names(Res_perm_bi)[match(res_adj.p[[1]][,1], pval_tmp)])
  return(res_adj.p)
} # end of function

# define function converting string to formula for Surv() function
surv_fml <- function(x, 
                     surv_time = surv_time, 
                     surv_status = surv_status) {
  fml = as.formula(paste(paste0("Surv(",as.name(surv_time),",",as.name(surv_status),")"),x,sep="~"))
  return(fml)
}

# define function of stratification for Cox model
cox_strat <- function(surv_data = surv_data,
                      surv_time = "surv_time",
                      surv_status = "surv_status",
                      strata_fct = strata_fct,
                      n_sample = 10,
                      cox_var = cox_var,
                      age.Acc = age.Acc,
                      model = model,
                      directory = directory){
  # remove level with sample fewer than 10
  strt_lvl = as.character(levels(as.factor(surv_data[,strata_fct])))
  lvl_keep = sapply(strt_lvl, function(level){
    keep = nrow(surv_data[surv_data[,strata_fct] == level,]) >= n_sample
    return(keep)
  })
  strt_lvl = strt_lvl[lvl_keep]
  # function for fitting Cox model at specified level
  fit_strata = lapply(strt_lvl, function(level = level){
    # find sample of specified level
    id_keep = surv_data[,strata_fct] == level
    id_keep[is.na(id_keep)] = F
    # retrieve subset for Cox model
    subdt = surv_data[id_keep, c(surv_time, surv_status, cox_var)]
    # fit Cox model
    fit = coxph(surv_fml(do.call(paste, c(as.list(cox_var), sep = "+")), surv_time, surv_status), data = subdt)
    # output the results
    name_strat = gsub("^paper_", "", strata_fct, fixed = F)
    fileName = paste0(directory,"summary_Cox_mod_", age.Acc, "_stratified_", name_strat, "_",model,".txt")
    capture.output(paste0("level: ", level), append = T, file = fileName)
    cat("\n", append = T, file = fileName)
    capture.output(summary(fit), append = T, file = fileName)
    #res = list(fit, subdt)
    return(fit)
  }) # end of function of level
  return(fit_strata)
} # end of function


# SECTION_4 miscellaneous -------------------------------------------------

# get sample_name from sample_ID
sample_name <- function(name.v = name.v, n.field = n.field) {
  str.mtx = sapply(strsplit(name.v, "-", fixed = T), "[", c(1:n.field))
  name = c()
  for (i in 1:length(name.v)) {
    name[i] = do.call(paste, c(as.list(str.mtx[,i]), sep = "-"))
  }
  return(name)
} 

# Format Number as Percentage
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(x * 100, format = format, digits = digits, ...), "%")
}

