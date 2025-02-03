#
# Moiety degree distributions
#
# input: a SNA edgelist stack
#        nodes file without the VMD indices (only cols 1-6)
# output: multi-column file with the degree frequencies for each moiety 
#         first column is degree
#

library(sna)
library(doParallel)
ncores = 40
#this nodes file doesn't have the VMD indices
nodes.file<-"data/CSH_CSSC_3to7_cyl_water_ions.psf.data.nodes"

infile<-"data/CSH_CSSC_3to7_cyl_every50.edgel.stack.rda"

outfile<-"data/degree-clargest-CSH_CSSC_3to7_cyl_every50.dat"

#initial frame
begframe<-1
#final frame (leave it as -1 if you want to go to the end of the stack)
endframe= -1

#### end of UI ##########

nodes<-read.table(nodes.file,col.names  = c("id","resname","segid","resid","name","type"))
load(infile)

if(endframe==-1)
  endframe<-length(gs)

dd<-degree(gs,g=begframe:endframe,gmode = "graph")
maxdeg<-max(dd)

mynames<-levels(as.factor(nodes$name))
mymat<-matrix(0,nrow=length(mynames),ncol=maxdeg+1)
row.names(mymat)<-mynames
colnames(mymat)<-0:maxdeg
for(i in begframe:endframe) {
  z<-gs[[i]]
  mycc<-which(component.largest(z))
  g<-get.inducedSubgraph(as.network(z),mycc)
  mycc.nodes<-nodes[mycc,]
  mycc.nodes$deg<-factor(degree(g,gmode = "graph"),levels=0:maxdeg)
  myt<-table(mycc.nodes$name,mycc.nodes$deg)
  mymat<-mymat+myt
}

mymat<-mymat/(endframe-begframe+1)

write.table(t(mymat),file=outfile,quote = F)

