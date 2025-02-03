library(sna)
library(doParallel)
registerDoParallel(cores = detectCores())
#setwd("/dfs2/tw/yuanmis1/mrsec/cg/analysis/csh_cssc_1to1_50mM_cg_big2/")
#this nodes file doesn't have the VMD indices
firstrun=2
lastrun=11
steps=100
#nodes.file<-"data/CSH_CSSC_3to7_cyl_water_ions.psf.data.nodes"
nodes.file<-"data/27A_redo.psf.data.nodes"
prefix<-"27A_3to7_cyl"

infile<-paste("data/", prefix, "_", firstrun, "to", lastrun,"_every",steps,".edgel.stack.rda",sep="")

outPre<-paste("data/resids/", "component_resid_", prefix, "_", firstrun, "to", lastrun,"_every",steps,sep="")

#infile<-"../csh_cssc_1to1_50mM_cg_big_every10_pt8and9.edgel.stack.rda"
#outPre<-"component_resid_csh_cssc_1to1_50mM_cg_5all_19to22_every10"

#--------------------- END of UI

#frame ranges
#330 to 470
#120 to 260
#550 to 670
#find out if any have more than one largest connected component 
#buba<-lapply(gs, function(z){components(component.largest(z,result = "graph"))})
#bubas<-unlist(buba)
#which(bubas>1)
# 1   9  12  16  17  32  35  44  46  48  49 265

nodes<-read.table(nodes.file)
colnames(nodes)<-c("id","resname","resid","molid","name","type")

load(infile)
gsize<-attr(gs[[1]],"n")


getResid<-function(g,nodes,howmany=1){
  cd<-component.dist(g)
#  myclus<-order(cd$csize,decreasing = T)[1:howmany]
#  myresids<-lapply(myclus,function(z,cd,nodes){
#    mynodes=which(cd$membership==z)
#    return(unique(nodes[mynodes,"resid"]))
#  },nodes=nodes,cd=cd)
#  rm(cd)
#  rm(myclus)
#  myresids
cd
}
#load(paste0(outPre,".rda"))
#a<-c()
#for (i in 1:length(resids)) {
#	if (length(resids[[i]][[1]])==0) {
#		a<-c(a,i)
#	}
#}
resids<-foreach(j=1:length(gs)) %dopar% (
#residstemp<-foreach(j=a) %dopar% (
	getResid(gs[[j]],nodes=nodes)
)
#c<-1
#for (i in a ) {
#	resids[[i]]<-residstemp[[c]]
#	c<-c+1
#}
#resids<-mclapply(gs,getResid,nodes=nodes,mc.cores=detectCores())
save(resids,file=paste0(outPre,".rda"))
#for (i in 1:length(gs)) {
#	outname<-getResid(gs[[i]],nodes=nodes)
#	save(outname, file=paste0(outPre,i,".rda"))
#	rm(outname)
#	print(paste0("frame",i))
#}



