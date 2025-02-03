firstrun=2
lastrun=9
steps=25
sysname="27A"
analysis="density_rda"
dir="data/"
ncores = 40

library(bio3d)
library(sna)
library(geometry)
library(doParallel)

ncores = 40
pdbname=paste(dir,sysname,".pdb",sep="")
mypdb<-read.pdb(pdbname)
dcdname=paste(dir,sysname,"_sum", firstrun, "to", lastrun,"_every_", steps, ".dcd",sep="")
mydcd<-read.dcd(dcdname)
nframes<-nrow(mydcd)
nframes<-52
load("data/component_resid_CSH_CSSC_3to7_cyl.rda")
out<-data.frame(matrix(ncol = 4, nrow = 0,byrow = T))
cluresid<-c()
k<-1
for(j in 1:nframes){
  myresids<-resids[[j]][[k]]
  cluresid[[j]]<-t(myresids)
  sel<-atom.select(mypdb,resno = myresids)
  particle<-matrix(mydcd[j,sel$xyz],ncol = 3,byrow = T)
  ab <- sum(myresids < 412 ) *2 +sum(myresids > 411 )
  tvol<-convhulln(particle, "FA")$vol
  rho<-ab/tvol
  out[nrow(out) + 1,]<-c(j,tvol, ab, rho)

}
outname=paste(dir,sysname,"_",analysis, "_",firstrun, "to", lastrun,".dat",sep="")

write.table(out, file=outname,col.names = F,row.names=F,quote = F)
