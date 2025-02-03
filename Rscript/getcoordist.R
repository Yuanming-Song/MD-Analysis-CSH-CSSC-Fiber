sysname="mix_27A"
analysis="density"
dir="data/"
steps=1
firstrunlist=c(1)
lastrunlist=c(8)

#firstrun=43
#lastrun=46
library(dplyr)
library(bio3d)
library(sna)
library(geometry)
#setwd("/dfs2/tw/yuanmis1/mrsec/cg/analysis/csh_cssc_1to1_50mM_cg_big")
source("/pub/jfreites/mrsec/analysis/networks/pbcs.R")
source("/pub/jfreites/mrsec/analysis/networks/shape/shape.R")
library(doParallel)
ncores = 40
pdbname=paste(dir,sysname,".pdb",sep="")
mypdb<-read.pdb(pdbname)
ipdbname=paste(sysname,".pdb",sep="")
a<-1
for (firstrun in firstrunlist) {
lastrun<-lastrunlist[a]
a<-a+1
dcdname=paste(dir,sysname,"_sum", firstrun, "to", lastrun,"_every_", steps, ".dcd",sep="")
mydcd<-read.dcd(dcdname)
nframes<-nrow(mydcd)
print(nframes)  

for(j in 1:nframes){
	                sel<-atom.select(mypdb,elena = c("N1",  "C1",  "O1",  "C2",  "C3",  "S1",   "N2",  "C4",  "O2",  "C5",  "C6",  "C7",  "C8",  "C9",  "C10",  "S2",  "C11",  "C12",  "C13",  "C14",  "N3",  "C15",  "C16",  "C17",  "C18",  "C19",  "C20",  "O3",  "O4",  "N4"))
	if (j==1) {	

		particle<-matrix(mydcd[j,sel$xyz],ncol = 3,byrow = T)
		print(dim(particle))
	} else {
		particle2<-matrix(mydcd[j,sel$xyz],ncol = 3,byrow = T)
		particle<-rbind(particle,particle2)
  	}
}
particle<-as.data.frame(particle)
binsizez <- 1
binnumz <- (0-min(particle$V3)+max(particle$V3))/binsizez - 1
zdist <- as.data.frame(table(cut(particle$V3,breaks =seq(-200,200,binsizez))))
zdistnorm <- zdist %>% mutate(x=seq(-200,200-binsizez,binsizez))
zdistnorm <- zdistnorm %>% mutate(y=Freq/sum(zdistnorm$Freq))


particle<-particle %>% mutate(r = sqrt(V1^2+V2^2))
binsizer <- 1
binnumr <- (0-min(particle$r)+max(particle$r))/binsizer - 1
rdist <- as.data.frame(table(cut(particle$r,breaks =seq(min(particle$r),max(particle$r),binsizer))))
rdistnorm <- rdist %>% mutate(x=seq(min(particle$r),max(particle$r)-binsizer,binsizer))
rdistnorm <- rdistnorm %>% mutate(y=Freq/(sum(rdistnorm$Freq)*((x+binsizer)^2-x^2)*3.14))
filename <- paste(sysname, "_zdist_all_",firstrun, "to", lastrun,".dat",sep="")
write.table(zdistnorm, file=filename,col.names = T,row.names=F,quote = F)
filename <- paste(sysname, "_rdist_all_",firstrun, "to", lastrun,".dat",sep="")
write.table(rdistnorm, file=filename,col.names = T,row.names=F,quote = F)

reslist=c("CSH", "CSSC")
for (i in reslist) {
	for(j in 1:nframes){
		sel<-atom.select(mypdb,elena = c("N1",  "C1",  "O1",  "C2",  "C3",  "S1",   "N2",  "C4",  "O2",  "C5",  "C6",  "C7",  "C8",  "C9",  "C10",  "S2",  "C11",  "C12",  "C13",  "C14",  "N3",  "C15",  "C16",  "C17",  "C18",  "C19",  "C20",  "O3",  "O4",  "N4"), resid = i)
		if (j==1) {

        	        particle1<-matrix(mydcd[j,sel$xyz],ncol = 3,byrow = T)
        	                                print(dim(particle1))


		} else {
        	        particle2<-matrix(mydcd[j,sel$xyz],ncol = 3,byrow = T)
        	        particle1<-rbind(particle1,particle2)
        	}
	}
	particle1<-as.data.frame(particle1)
	#zdist <- as.data.frame(table(cut(particle1$V3,breaks =seq(min(particle$V3),max(particle$V3),binsizez))))
	#zdistnorm1 <- zdist %>% mutate(x=seq(min(particle$V3),max(particle$V3)-binsizez,binsizez))
	zdist <- as.data.frame(table(cut(particle1$V3,breaks =seq(-200,200,binsizez))))
	zdistnorm1 <- zdist %>% mutate(x=seq(-200,200-binsizez,binsizez))
	zdistnorm1 <- zdistnorm1 %>% mutate(y=Freq/sum(zdistnorm1$Freq))


	particle1<-particle1 %>% mutate(r = sqrt(V1^2+V2^2))
	#rdist <- as.data.frame(table(cut(particle1$r,breaks =seq(min(particle$r),max(particle$r),binsizer))))
	#rdistnorm1 <- rdist %>% mutate(x=seq(min(particle$r),max(particle$r)-binsizer,binsizer))
	rdist <- as.data.frame(table(cut(particle1$r,breaks =seq(0,16,binsizer))))
	rdistnorm1 <- rdist %>% mutate(x=seq(0,16-binsizer,binsizer))
	rdistnorm1 <- rdistnorm1 %>% mutate(y=Freq/(sum(rdistnorm1$Freq)*((x+binsizer)^2-x^2)*3.14))
	filenamez=paste("run2_zdist_",i,"_",firstrun, "to", lastrun,".dat",sep="")
	filenamer=paste("run2_rdist_",i,"_",firstrun, "to", lastrun,".dat",sep="")
	write.table(zdistnorm1, file=filenamez,col.names = T,row.names=F,quote = F)
	write.table(rdistnorm1, file=filenamer,col.names = T,row.names=F,quote = F)

}
}
