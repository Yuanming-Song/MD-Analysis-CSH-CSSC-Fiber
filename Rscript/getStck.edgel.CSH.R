#
# convert single-edgelist files (.edges) to a sna edgelist snd->rcv stack
# this is the multiple file version
# 
#
library(doParallel)
registerDoParallel(cores = detectCores())
nnodes<-4696
firstrun=6
lastrun=10
myDir<-"edge/"
steps=1
thePrefixes<-"mix_27A_3to7_cyl_CSH_every_10"
begframes<-1
endframes<-1775
outname<-paste("data/", thePrefixes, sep="")

####################################### end of UI
nframes<- sum(endframes-begframes)+length(begframes)
treatedgelist<-function(j) {
  if(j<10) {
    framenum<-paste0("00",j)
  } else {
    if(j<100) {
      framenum<-paste0("0",j)
    } else {
      framenum<-j
    }
  }
  edgefile.name<-paste0(myDir,thePrefixes[i],"_",framenum,".edges")
  if (file.exists(edgefile.name)) {        
    edgelist<-matrix(scan(edgefile.name),ncol=3,byrow=T)
    edgelist[,3]<-1
    edgelist<-rbind(edgelist,edgelist[,c(2,1,3)])
    edgelist<-edgelist[!duplicated(edgelist),]
    attr(edgelist,"n")<-nnodes
    edgelist
  }
}
for(i in 1:length(begframes)) {
  gs<-foreach(j=begframes:endframes) %dopar% (
    treatedgelist(j)
  )
  gs[sapply(gs, is.null)] <- NULL
  save(gs,file=paste0(outname,".edgel.stack.rda"))
}


