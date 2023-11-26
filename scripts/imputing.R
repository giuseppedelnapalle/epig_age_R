# adapted from Horvath's script

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


noMissingPerSample0 = apply(as.matrix(is.na(datMethUsed0)),1,sum)
table(noMissingPerSample0)

#STEP 2: Imputing 
if (max(noMissingPerSample0,na.rm=TRUE)>=6000) print("max # of missing value per sample over 6000")

if (! fastImputation0 & nSamples>1 & max(noMissingPerSample0,na.rm=TRUE)< 6000 ){

# run the following code if there is at least one missing
if ( max(noMissingPerSample0,na.rm=TRUE)>0 ){
dimnames0 = dimnames(datMethUsed0)
datMethUsed0 = data.frame(t(impute.knn(data=t(datMethUsed0), k=10, rowmax=.6)$data))
dimnames(datMethUsed0)=dimnames0
} # end of if
} # end of if (! fastImputation0 )

# update noMissingPerSample
noMissingPerSample0_2 = apply(as.matrix(is.na(datMethUsed0)),1,sum)

