source("Rscript/package.R")
source("Rscript/base_fun.R")
step<-100


source("Rscript/getDensity_base.R")
firstrun<-1
outputdir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/mix_27A/analysis/data/"
midfix<-""
if (FALSE) {
  if (TRUE) {
    ### how many residues in total
    mincssc<-1
    maxcssc<-411 #411
    mincsh<-maxcssc+1
    maxcsh<-763
    totcssc<-maxcssc
    totcsh<-maxcsh-maxcssc
    dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/mix_27A/"
    # mypdb<-read.pdb(file.path(dcddir, "setup/CSH_CSSC_mix_27A_cyl_water_ions.pdb"),hex = TRUE)
    pdbname<-"data/mix_27A.pdb"
    mypdb<-read.pdb(pdbname)
    outnamepre<-"mix_27A_fibre_rho_network_"
    uniseltextlist<-c("CSSC","CSH","TIP3") #unique resname in sim
    if (TRUE) {
      midfix<-"cylcssc"
      outnamepre<-paste0(outnamepre,midfix,"_")
      midfix<-paste0(midfix,".")
      firstrun<-12
      
    }
  } else {
    dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/cssc_27A/"
    mypdb<-read.pdb(file.path(dcddir, "setup/CSSC_pure_27A_cyl_water_ions.pdb"),hex = TRUE)
    mypdb$atom<-mypdb$atom[which(mypdb$atom$resid=="CSSC"),]
    maxcssc<-max(mypdb$atom$resno)
    uniseltextlist<-c("CSSC","TIP3") #unique resname in sim
    outnamepre<-"pure_27A_fibre_rho_network_"
  }
} else {
  dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/CSH_27A/"
  mypdb<-read.pdb("/dfs9/tw/yuanmis1/mrsec/cyl1/run8/CSH_3to7_cyl_water_ions.pdb",hex = TRUE)
  mypdb$atom<-mypdb$atom[which(mypdb$atom$resid=="CSH"),]
  maxcssc<-0
  uniseltextlist<-c("CSH","TIP3") #unique resname in sim
  outnamepre<-"CSH_27A_fibre_rho_network_"
  firstrun<-2
}
source("Rscript/getCOM_moiety_base.R")

print(outnamepre)
dcdprefix<-"nvt"
lastrun<-100
dcdlist<-getdcdlist(firstrun,midfix = midfix)

rhistot<-0 #empty matrix
zhistot<-0 #empty matrix
fibdist<-list()
processedframe<-0
rholist<-c()
processedframe<-0
oldoutname<-"randomplaceholder"
for (dcdname in dcdlist) {
  print(c(dcdname,format(Sys.time(), "%Y-%m-%d %H:%M:%S"))) #print time
  mydcd<-read.dcd(file.path(dcddir,dcdname),verbose = FALSE,big = TRUE) #read dcd file
  mydcdcell<-read.dcd(file.path(dcddir,dcdname),verbose = FALSE,cell = TRUE) #read box size
  totframe<-nrow(mydcdcell) #how many frames
  print(c(totframe,format(Sys.time(), "%Y-%m-%d %H:%M:%S"))) #print time again
  #parallel compute with multiple frames
  
  if (length(seq(1,totframe,step))>1) { 
    rholisttemp<-foreach(frame=seq(1,totframe,step),.combine=rbind) %dopar% (
      getrho_fibre_network(frame)
    )
    rholisttemp[,1]<-rholisttemp[,1]+processedframe
  } else {
    rholisttemp<-getrho_fibre_network(1)
    rholisttemp[1]<-rholisttemp[1]+processedframe
    
  }
  print(c(length(seq(1,totframe,step)),format(Sys.time(), "%Y-%m-%d %H:%M:%S"))) #print time again 
  processedframe<-processedframe+totframe
  rholist<-rbind(rholist,rholisttemp)
  outname<-paste0(outputdir,outnamepre,firstrun,"to",regmatches(dcdname, regexpr("\\d+", dcdname)),".rda")
  save(rholist,file=outname)
  print(outname)
  if (file.exists(oldoutname)) {
    file.remove(oldoutname)
  }
  oldoutname<-outname
}
colnames(rholist)<-c("t","rho","totres","nCSH")
save(rholist,file=outname)