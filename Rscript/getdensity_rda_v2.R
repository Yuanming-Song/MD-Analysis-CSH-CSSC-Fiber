firstrun=2
lastrun=11
steps=100
sysname="27A"
prefix<-"27A_3to7_cyl"
analysis="density_rda"
dir="data/"
CSSCmax=411

library(bio3d)
library(sna)
library(geometry)
library(doParallel)
pdbname=paste(dir,sysname,".pdb",sep="")
dcdname=paste(dir,sysname,"_sum", firstrun, "to", lastrun,"_every_", steps, ".dcd",sep="")
mypdb<-read.pdb(pdbname)
mydcd<-read.dcd(dcdname)
nframes<-nrow(mydcd)
nframes<-52
registerDoParallel(cores = detectCores())
rdaname=paste0("data/", "component_resid_", prefix, "_", firstrun, "to", lastrun,"_every",steps,".rda")
load(rdaname)
#cluresid<-c()
k<-1
getdensity<-function(j) {
  myresids<-resids[[j]][[k]]
  #cluresid[[j]]<-t(myresids)
  sel<-atom.select(mypdb,resno = myresids)
  particle<-matrix(mydcd[j,sel$xyz],ncol = 3,byrow = T)
  ab <- sum(myresids <= CSSCmax ) *2 +sum(myresids > CSSCmax )
  tvol<-convhulln(particle, "FA")$vol
  rho<-ab/tvol
  c(j,tvol, ab, rho)
}
out<-foreach(j=1:nframes, .combine=rbind) %dopar% (
  getdensity(j)
)
outname=paste(dir,sysname,"_",analysis, "_",firstrun, "to", lastrun,"_every",steps,".dat",sep="")

write.table(out, file=outname,col.names = F,row.names=F,quote = F)
