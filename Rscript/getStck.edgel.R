#
# convert single-edgelist files (.edges) to a sna edgelist snd->rcv stack
#
# 
#
firstrun=2
lastrun=7
steps=25
begframe<-1
endframe<-295
myDir<-"edge/"
prefix<-"27A_3to7_cyl"
outname<-paste("data/", prefix, "_", firstrun, "to", lastrun,"_every",steps,sep="")
nnodes<-1656

####################################### end of UI
nframes<- endframe-begframe+1

gs<-vector("list",nframes)

for(j in 1:nframes) {
        k<-j+begframe-1
    if(k<10)
      framenum<-paste0("00",k)
    else
      if(k<100)
        framenum<-paste0("0",k)
      else
        framenum<-k
    edgefile.name<-paste0(myDir,prefix,"_",framenum,".edges")
    edgelist<-matrix(scan(edgefile.name),ncol=3,byrow=T)
    edgelist[,3]<-1
    edgelist<-rbind(edgelist,edgelist[,c(2,1,3)])
    edgelist<-edgelist[!duplicated(edgelist),]
    attr(edgelist,"n")<-nnodes
    gs[[j]]<-edgelist
}
save(gs,file=paste0(outname,".edgel.stack.rda"))
