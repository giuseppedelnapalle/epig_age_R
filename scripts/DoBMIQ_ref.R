### DoBMIQ.R
load("./data.m.Rd")
load("./probeInfoALL.lv.Rd")
load("./anno450k.Rd")
source("./BMIQ_1.4.R")

print("Read Success.")

index <- which(probeInfoALL.lv$probeID %in% rownames(data.m))
index <- index[match(rownames(data.m),probeInfoALL.lv$probeID[index])]


###
type1.idx <- which(probeInfoALL.lv[[2]][index]==1);
type2.idx <- which(probeInfoALL.lv[[2]][index]==2);
design.v <- probeInfoALL.lv[[2]][index];


pdf("Profiles.pdf",width=4,height=3);
for(s in 1:ncol(data.m)){
   plot(density(data.m[type1.idx,s]));
   d.o <- density(data.m[type2.idx,s]);
   points(d.o$x,d.o$y,type="l",col="red");
print(s);
}
dev.off();


for(s in 1:ncol(data.m)){
  beta.v <- data.m[,s];
  bmiq.o <- BMIQ(beta.v,design.v,sampleID=s);
  tmp.v <- bmiq.o$nbeta;
  save(tmp.v,file=paste("bmiq",s,".Rd",sep=""));
  print(paste("Done BMIQ for sample ",s,sep=""));
}

bmiq.m <- data.m;
rm(data.m);
for(s in 1:ncol(bmiq.m)){
  load(paste("bmiq",s,".Rd",sep=""));
  bmiq.m[,s] <- tmp.v;
  print(s);
}

save(bmiq.m,file="bmiq.Rd");
### select CpGs (remove CHs)
annoindex <- which(anno450k.m[,1] %in% rownames(bmiq.m))
annoindex <- annoindex[match(rownames(bmiq.m),anno450k.m[annoindex])]
anno.m <- anno450k.m[annoindex,]

cg.idx <- grep("cg",anno.m[,1]);
anno2.m <- anno.m[cg.idx,];
bmiq2.m <- bmiq.m[cg.idx,];
rm(bmiq.m);
save(bmiq2.m,file="BMIQ_RESULT.Rd")
save(anno2.m,file="anno2.m.Rd")
