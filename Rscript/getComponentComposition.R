#
# largest connected component composition (a contingecy table)
#
#
# input: a SNA edgelist stack
#        nodes file without the VMD indices (only cols 1-6)
# output: multi-column file with  quantities proportional to 
#         the molar fractions of the largest connected component 
#	  the first column is a running index
#

library(sna)
library(doParallel)
ncores = 40
#this nodes file doesn't have the VMD indices
nodes.file<-"data/run2_redo.psf.data.nodes"


infile<-"data/CSH_CSSC_3to7_cyl_every50.edgel.stack.rda"

outfile<-"data/composition-CSH_CSSC_3to7_cyl_every50.dat"

#multiply the largest CC proportion of the graph size by this number
scaleFactor<-414

#### end of UI ##########

nodes<-read.table(nodes.file)

load(infile)
gsize<-attr(gs[[1]],"n")

compo<-sapply(gs, function(z){
  cd<-component.dist(z);
  mymax<-ifelse(cd$membership %in% which(cd$csize==max(cd$csize)),cd$membership,NA);
  table(nodes$V2,mymax)})

compo.evol<-matrix( unlist(compo),ncol =2,byrow = T)*(scaleFactor/gsize)
write.table(compo.evol,file=outfile,col.names = F,quote = F)



