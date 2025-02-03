#
# convert single-edgelist files (.edges) to a sna edgelist snd->rcv stack
# this is the multiple file version
# 
#
library(doParallel)
registerDoParallel(cores = detectCores())
nnodes<-3522
load("mix_27A_edges_COM.rda")


treatedgelist<-function(framenum) {
  edgelist<-edges[[framenum]][[1]]
  edgelist<-rbind(edgelist,edgelist[,c(2,1,3)])
  edgelist<-edgelist[!duplicated(edgelist),]
  attr(edgelist,"n")<-nnodes
  edgelist
}
gs<-foreach(j=1:length(edges)) %dopar% (
  treatedgelist(j)
)
gs[sapply(gs, is.null)] <- NULL
save(gs,file=paste0("mix_27A_COM.edgel.stack.rda"))



