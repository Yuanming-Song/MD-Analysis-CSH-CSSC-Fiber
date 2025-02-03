firstrun=18
lastrun=43
totalres=269
steps=25
sysname="run2"
analysis="phe_com"
dir="data/"
ncores = 40

outname=paste("test_",sysname,"_",analysis, "_",firstrun, "to", lastrun,".dat",sep="")


#####   PHE list for each residue
reslist=c("CSH", "CSSC")
minCSSC <- 1
maxCSSC <- 145
moiety1<-c("C5","C6","C7","C8","C9","C10")
moiety2<-c("C12","C13","C16","C19","C14","C18")

library(dplyr)
library(bio3d)
library(sna)
library(geometry)
library(pracma)
library(doParallel)
library(plyr)

#####   read pdb and dcd files
pdbname=paste(dir,sysname,".pdb",sep="")
mypdb<-read.pdb(pdbname)
dcdname=paste(dir,sysname,"_sum", firstrun, "to", lastrun,"_every_", steps, ".dcd",sep="")
mydcd<-read.dcd(dcdname)
nframes<-nrow(mydcd)
#####	for testing...
#nframes<-1

#####   function for calculating distance
euclidean <- function(a, b) sqrt(sum((a - b)^2))

#####   function for getting COM and update list
getdist <- function() {
  selphei<-atom.select(mypdb,elety = moietylist,resno =i)
  selpheixyz<-mydcd[j,selphei$xyz]
  comi<-as.vector(com.xyz(selpheixyz))
  for (k in i:totalres) {
    #####   avoid same residue
    if (i==k) {
      next
    }
    selphek<-atom.select(mypdb,elety = moiety1,resno =k)
    selphekxyz<-mydcd[j,selphek$xyz]
    comk<-as.vector(com.xyz(selphekxyz))    
    dis<-euclidean(comi,comk)
    dislist<<-c(dislist,dis)
    if (k>=minCSSC && k<=maxCSSC) {
      selphek<-atom.select(mypdb,elety = moiety2,resno =k)
      selphekxyz<-mydcd[j,selphek$xyz]
      comk<-as.vector(com.xyz(selphekxyz))    
      dis<-euclidean(comi,comk)
      dislist<<-c(dislist,dis)    
    }
  }
}




dislist<-c()
for (j in 1:nframes) {
  for (i in 1:totalres) {
    moietylist<-moiety1
    getdist()
    if (i>=minCSSC && i<=maxCSSC) {
      moietylist<-moiety2
      getdist()
    }
  }
}   
write.table(as.data.frame(dislist), file=outname,col.names = F,row.names=F,quote = F)

analysis="rdf"

binsize<-1
rdf<-as.data.frame(table(cut(distlist,breaks =seq(0,200,binsize))))
rdfnorm<-rdf %>% mutate(x=seq(0,200-binsize,binsize))
rdfnorm<-rdfnorm %>% mutate(y=Freq/(sum(Freq)*((x+binsize)^2-x^2)*3.14))

outname=paste("test_",sysname,"_",analysis, "_",firstrun, "to", lastrun,".dat",sep="")
write.table(rdfnorm, file=outname,col.names = T,row.names=F,quote = F)


# binsize <- (0-min(distlist)+max(distlist))/(binnum+1)
# disthis <- as.data.frame(table(cut(distlist,breaks =seq(min(distlist),max(distlist),binsize))))
# library(dplyr)
# disthisnorm <- disthis %>% mutate(x=seq(min(disthis$V1),max(disthis$V1)-1,binw))
# #rden<-rden %>% mutate(y=Freq/(sum(Freq)*((x+binw)^2-x^2)*3.14))
# disthisnorm <- disthisnorm %>% mutate(y=Freq/(sum(disthis$Freq)*((x+binw)^3-x^3)*4*3.14/3))
# filename <- paste(prefix, "all_his",".dat",sep="")
# write.table(as.data.frame(disthisnorm), file=filename,col.names = T,row.names=F,quote = F)

