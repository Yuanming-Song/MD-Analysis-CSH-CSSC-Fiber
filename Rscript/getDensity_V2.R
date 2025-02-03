source("Rscript/package.R")
source("Rscript/base_fun.R")
step<-10
pdbname<-"data/mix_27A.pdb"
mypdb<-read.pdb(pdbname)
dcdname<-paste0("mix_27A_sum_every",step,".dcd")
mydcd<-read.dcd(dcdname,verbose=FALSE,big=TRUE)
mydcdcell<-read.dcd(dcdname,cell=T,verbose=FALSE)
### how many residues in total
mincssc<-1
maxcssc<-411 #411
mincsh<-maxcssc+1
maxcsh<-763
totcssc<-maxcssc
totcsh<-maxcsh-maxcssc
perbatch<-40
resname<-c("CSH","CSSC")
source("Rscript/getDensity_base.R")
source("Rscript/loadCOM_base.R")
totframe<-nrow(mydcd)
rholist<-c()
for (batch in 1:ceiling(totframe/perbatch)) {
  iniframe<-(batch-1)*perbatch+1
  finiframe<-ifelse(batch<ceiling(totframe/perbatch),batch*perbatch,totframe)
  rholisttemp<-foreach(frame=iniframe:finiframe,.combine=rbind) %dopar% (
    getrho_fibre_network(frame)
  )
  rholist<-rbind(rholist,rholisttemp)
  save(rholist,file="data/mix_27A_density_COM_net.rda")
  print(batch)
}
