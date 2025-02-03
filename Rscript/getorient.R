startframe<-1
sysname="mix_27A"
analysis="density"
dir="data/"
steps=25
firstrun=6
lastrun=10
mincssc<-1
maxcssc<-411 #411
mincsh<-maxcssc+1
maxcsh<-763
totalres<-maxcsh
### for pdb & traj file names
pdbname=paste0("data/CSH_CSSC_mix_27A_cyl_water_ions_redo.pdb")
#pdbname=paste(dir,sysname,"_resonly.pdb",sep="")
dcdname=paste(dir,sysname,"_sum", firstrun, "to", lastrun,"_all_every_", steps, ".dcd",sep="")
binsize<-5
binmax<-200
binmin<--200
nbins <- ceiling((binmax - binmin) / binsize)
histx<-seq(binmin + binsize / 2, binmax - binsize / 2, by = binsize)
moietyname <- c("N1", "C1", "O1", "C2", "C3", "S", "N2", "C4", "O2", "C5", "C6", "C7", "C8", "C9", "C10")
moietyname2 <- c("S2", "C11", "C12", "C13", "C14", "N3", "C15", "C16", "C17", "C18", "C19", "C20", "O3", "O4", "N4")
phemoietyname <- c("C5", "C6", "C7", "C8", "C9", "C10")
phemoietyname2 <- c("C12", "C13", "C14", "C16", "C18", "C19")

library(dplyr)
library(bio3d)
library(sna)
library(geometry)
library(pracma)
library(doParallel)
library(plyr)
registerDoParallel(cores = detectCores())
getcostheta <- function(v1, v2, v3,vref) {
	va <- v2 - v1
	vb <- v3 - v1
	vc<-cross(va,vb)
	vref<-as.vector(vref)
	sum(vc * vref)/(sqrt(sum(vref^2))*sqrt(sum(vc^2)))

}
getunit<-	function(v1, v2, v3) {
  va <- v2-v1
  va/sqrt(sum(va^2))
  
}
cbind.fill <- function(...){
    nm <- list(...) 
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}
vref<-c(0,0,1)
mycosthetalist <- c()
mycoorlist<-c()
update_orient_table <- function(table, zcoor,zvec) {
  if (zcoor < binmin) {
    table[1, 2] <- table[1, 2] + zvec
  } else if (zcoor > binmax) {
    table[nrow(table), 2] <- table[nrow(table), 2] + zvec
  } else {
    bin <- floor((zcoor - binmin) / binsize) + 1
    table[bin, 2] <- table[bin, 2] + zvec
  }
  table
}  
getphezdens<-function(framenum) {
  temphistogram <- matrix(0, nrow = nbins, ncol = 2)
  temphistogram[, 1] <- histx
  for(i in 1:totalres) {
    phezsel<-atom.select(mypdb,elety = phemoietyname,resno =i)
    phexyz<-mydcd[framenum,phezsel$xyz]
    phecom<-as.vector(com.xyz(phexyz))-com_t[framenum]
    phecomz<-phecom[3]
    sel1<-atom.select(mypdb,elety = c("C5"),resno =i)
    sel2<-atom.select(mypdb,elety = c("C8"),resno =i)
    sel3<-atom.select(mypdb,elety = c("C6"),resno =i)
    v1<-as.vector(mydcd[framenum,sel1$xyz])
    v2<-as.vector(mydcd[framenum,sel2$xyz])
    v3<-as.vector(mydcd[framenum,sel3$xyz])
    temphistogram<-update_orient_table(temphistogram,phecomz,getunit(v1, v2, v3)[3])
    if (i<=maxcssc) {
      phezsel<-atom.select(mypdb,elety = phemoietyname,resno =i)
      phexyz<-mydcd[framenum,phezsel$xyz]
      phecom<-as.vector(com.xyz(phexyz))-com_t[framenum]
      phecomz<-phecom[3]
      sel1<-atom.select(mypdb,elety = c("C16"),resno =i)
      sel2<-atom.select(mypdb,elety = c("C14"),resno =i)
      sel3<-atom.select(mypdb,elety = c("C13"),resno =i)	
      v1<-as.vector(mydcd[framenum,sel1$xyz])
      v2<-as.vector(mydcd[framenum,sel2$xyz])
      v3<-as.vector(mydcd[framenum,sel3$xyz])
      temphistogram<-update_orient_table(temphistogram,phecomz,getunit(v1, v2, v3)[3])
    }
  }
  temphistogram[,2]
}
#####   read dcd files
mydcd<-read.dcd(dcdname)
mypdb<-read.pdb(pdbname,hex=TRUE)
nframes<-nrow(mydcd)
resnoh<-atom.select(mypdb,elety=c(moietyname,moietyname2))
resnoh<-atom.select(mypdb,elety=c(moietyname,moietyname2))

for (framenum in startframe:nframes) {
  if (framenum==startframe) {
    com_t<-matrix(nrow=nframes*3,ncol = 3)
    counter<-1
    coorlab<-c("x","y","z")
  }
  #####   get coordinate
  resnohxyz<-mydcd[framenum,resnoh$xyz]
  #####   calculate COM of res i
  comi<-as.vector(com.xyz(resnohxyz))
  ##### for tracking com shift
  for (i in 1:3) {
    com_t[framenum,]<-comi
    counter<-counter+1
  }
}

phez <- foreach (framenum = startframe:nframes, .combine = "+") %dopar% {
  getphezdens(framenum)
}
phez <- cbind(histx, phez)
colnames(phez)<-c("z","Sum")
write.table(phez, file=paste0(sysname,"_phe_orient", ".dat"), row.names=FALSE,quote = F)

if (0) {
phez<-read.table(file="~/Documents/Research/HPC/dfs2/mrsec/cyl1/mix_27A/analysis/mix_27A_phe_orient.dat",header=TRUE)
ggplot()+geom_line(data = phez,aes(x=z,y=Sum/200))+ylab("Averge z component of unit vector per Frame")
}