# adapted from Horvath's script
# find optimal # of neighbour for min noMissingPerSample

#"
# Steve Horvath: Estimating DNAm age.
# This file assumes a data frame exists called dat1 whose rows correspond to CpGs
# and whose first column reports the CpG identifier
# and whose remaining columns corresponds to samples (e.g. Illumina arrays).
#"

fastImputation0=FALSE


#STEP 1: DEFINE QUALITY METRICS

meanMethBySample0 = as.numeric(apply(as.matrix(dat2[,-1]),2,mean,na.rm=TRUE))
minMethBySample0 = as.numeric(apply(as.matrix(dat2[,-1]),2,min,na.rm=TRUE))
maxMethBySample0 = as.numeric(apply(as.matrix(dat2[,-1]),2,max,na.rm=TRUE))

datMethUsed0 = t(dat2[,-1])
colnames(datMethUsed0) = as.character(dat2[,1])

# noMissingPerSample
noMissingPerSample0 = apply(as.matrix(is.na(datMethUsed0)),1,sum)
print(table(noMissingPerSample0))
print(paste("mean # of missing value per sample:",mean(noMissingPerSample0,na.rm=TRUE)))
print(paste("median # of missing value per sample:",median(noMissingPerSample0,na.rm=TRUE)))
print(paste("max # of missing value per sample:",max(noMissingPerSample0,na.rm=TRUE)))
print(paste("min # of missing value per sample:",min(noMissingPerSample0,na.rm=TRUE)))
cat("\n")

# noMissingPerProbe
noMissingPerProbe0 = apply(as.matrix(is.na(datMethUsed0)),2,sum)
print(table(noMissingPerProbe0))
print(paste("mean # of missing value per probe:",mean(noMissingPerProbe0,na.rm=TRUE)))
print(paste("median # of missing value per probe:",median(noMissingPerProbe0,na.rm=TRUE)))
print(paste("max # of missing value per probe:",max(noMissingPerProbe0,na.rm=TRUE)))
print(paste("min # of missing value per probe:",min(noMissingPerProbe0,na.rm=TRUE)))
cat("\n")


#STEP 2: Imputing 
if (max(noMissingPerSample0,na.rm=TRUE)>=3000) print("max # of missing value per sample over 3000, imputing not implemented.")

if (! fastImputation0 & nSamples>1 & max(noMissingPerSample0,na.rm=TRUE)< 3000 ){

# run the following code if there is at least one missing
  mean.na = c(); median.na = c(); max.na = c(); min.na = c();
if ( max(noMissingPerSample0,na.rm=TRUE)>0 ){
#dimnames0 = dimnames(datMethUsed0)
for (i in c(1:20)) {
  datMethUsed0 = data.frame(t(impute.knn(data=t(datMethUsed0), k=i+1, rowmax=.6)$data))
  noMissingPerSample0_2 = apply(as.matrix(is.na(datMethUsed0)),1,sum)
  mean.na[i] = mean(noMissingPerSample0_2,na.rm=TRUE)
  median.na[i] = median(noMissingPerSample0_2,na.rm=TRUE)
  max.na[i] = max(noMissingPerSample0_2,na.rm=TRUE)
  min.na[i] = min(noMissingPerSample0_2,na.rm=TRUE)
  print(paste("k:", i+1,";","statistics:","mean:",mean.na[i],"median:",median.na[i],"max:",max.na[i],"min:",min.na[i]))
}
#dimnames(datMethUsed0)=dimnames0
  stat.na = data.frame(k = c(2:21), mean.na, median.na, max.na, min.na)
  print("Program finieshed.")
} # end of if
} # end of if (! fastImputation0 )


