### DoBMIQ.R for epiTOC analysis
# adpated from script DoBMIQ.R of Teschendorff lab


# 1 load required data & function -----------------------------------------

dir_450k <- "~/Project_Data/R_data/Illum_450k_manifest/"
rm(BMIQ)
source(file = paste0(dir_fun, "BMIQ_1.4_mdf.R"))
load(file = paste0(dir_450k, "probe_type_450k.RData"))
#load("./data.m.Rd")
#load("./anno450k.Rd")
print("Read Success.")


# 2 obtain parameter design.v ---------------------------------------------

# sorted index of probe corresponding to DNAm matrix
data.m <- t(datMethUsed0)
index <- which(probe_type$Name %in% rownames(data.m))
index <- index[match(rownames(data.m),probe_type$Name[index])]

# index of type_I & type_II
type1.idx <- which(probe_type$Type[index] == "I")
type2.idx <- which(probe_type$Type[index] == "II")

# obtain probe type of probes in DNAm matrix
design.v <- probe_type$Type[index]
design.v <- ifelse(design.v == "I", 1, 2)

# plot density
dir.create(paste0(dir_plot2, "BMIQ"))
pdf(file = paste0(dir_plot2,"BMIQ/","Profiles_probe_type_epiTOC.pdf"),
    width = 8, height = 6, pointsize = p.size_plt)
for(s in 1:ncol(data.m)){
   plot(density(data.m[type1.idx,s]))
   d.o <- density(data.m[type2.idx,s])
   points(d.o$x,d.o$y,type="l",col="red")
print(s)
}
dev.off()


# 3 perform BMIQ ----------------------------------------------------------

dir.create(paste0(dir_rdata2, "BMIQ"))
for(s in 1:ncol(data.m)){
  beta.v <- data.m[,s]
  bmiq.o <- BMIQ(beta.v,design.v,sampleID=s, dir_plot=paste0(dir_plot2, "BMIQ"))
  tmp.v <- bmiq.o$nbeta
  save(tmp.v,file=paste0(dir_rdata2,"BMIQ/", "bmiq",s,".Rd"))
  print(paste("Done BMIQ for sample ",s,sep=""))
}

bmiq.m <- data.m
rm(data.m)
for(s in 1:ncol(bmiq.m)){
  load(paste0(dir_rdata2,"BMIQ/", "bmiq",s,".Rd"))
  bmiq.m[,s] <- tmp.v
  print(s)
}

save(bmiq.m,file=paste0(dir_rdata2,"BMIQ/","bmiq.Rd"))
rm(tmp.v)
