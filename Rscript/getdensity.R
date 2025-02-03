sysname="27A"
analysis="density"
dir="data/"
steps=25
firstrun=2
lastrun=7
#set directory 
#setwd("/dfs2/tw/yuanmis1/mrsec/cg/analysis/csh_cssc_1to1_50mM_cg_big")
library(bio3d)
library(sna)
library(geometry)
library(doParallel)
ncores = 40
#pbc wrapped residue/protein only pdb and dcd files  
pdbname=paste(dir,sysname,".pdb",sep="")  
mypdb<-read.pdb(pdbname)
pdbname=paste(sysname,".pdb",sep="")
dcdname=paste(dir,sysname,"_sum", firstrun, "to", lastrun,"_every_", steps, ".dcd",sep="")
print(dcdname)
mydcd<-read.dcd(dcdname)
nframes<-nrow(mydcd)
colnames<-c("steps","Volume","NUnit","Density")
out<-data.frame(matrix(ncol = length(colnames), nrow = 0,byrow = T))
#names(out)<-a
  for(j in 1:nframes){
	sel<-atom.select(mypdb,elety = c("N1",  "C1",  "O1",  "C2",  "C3",  "S1",   "N2",  "C4",  "O2",  "C5",  "C6",  "C7",  "C8",  "C9",  "C10",  "S2",  "C11",  "C12",  "C13",  "C14",  "N3",  "C15",  "C16",  "C17",  "C18",  "C19",  "C20",  "O3",  "O4",  "N4"))
    	particle<-matrix(mydcd[j,sel$xyz],ncol = 3,byrow = T)
	units<- 1174
	tvol<-convhulln(particle, "FA")$vol
	rho<-units/tvol
	out[nrow(out) + 1,]<-c(j,tvol, units, rho)
  }
outname=paste(dir,sysname,"_",analysis, "_",firstrun, "to", lastrun,".dat",sep="")
write.table(out, file=outname,col.names = T,row.names=F,quote = F)
