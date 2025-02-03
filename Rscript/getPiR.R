source("Rscript/package.R")
source("Rscript/base_fun.R")
abs<-TRUE
outputdir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/mix_27A/analysis/data/"
if (TRUE) {
  outnamepre<-"mix_27A_cylcssc_fibre_Pi"
  dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/mix_27A/"
  mypdb<-read.pdb(file.path(dcddir, "setup/CSH_CSSC_mix_27A_cyl_water_ions.pdb"),hex = TRUE)
  UniqResnameList<-c("CSH","CSSC")
  maxCSH<-max(mypdb$atom[which(mypdb$atom$resid=="CSH"),]$resno)
  if (FALSE) {
    outnamepre<-"mix_27A_cylcssc_fibre_Pi"
    midfix<-"cylcssc."
  } else {
    outnamepre<-"mix_27A_fibre_Pi"
    midfix<-""
  }
} else {
  dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/cssc_27A/"
  mypdb<-read.pdb(file.path(dcddir, "setup/CSSC_pure_27A_cyl_water_ions.pdb"),hex = TRUE)
  UniqResnameList<-c("CSSC") #unique resname in sim
  outnamepre<-"pure_27A_fibre_Pi"
  midfix<-""
}
maxCSSC<-max(mypdb$atom[which(mypdb$atom$resid=="CSSC"),]$resno)
fibresel<-atom.select(mypdb,resid=UniqResnameList) #select CSH/CSSC atoms only
dcdprefix<-"nvt"
firstrun<-12
lastrun<-100
step<-10
dcdlist<-getdcdlist(firstrun,midfix=midfix)
pdis<-4.9 #parallel distance maximum
pang<-0.9 #parallel angle cutoff
tdisl<-5.5 #t shape distance minimum
tdish<-5.8 #t shape distance maximum
tang<-0.23 #t shape angle
PHEname<-list()
PHEname[[1]]<-c("C5","C8","C6")
PHEname[[2]]<-c("C16","C14","C13")
PHElistC<-list()
for (i in 1:2) {
  PHElistC[[i]]<-list()
  for (j in 1:3) {
    PHEsel<-atom.select(mypdb,elety=PHEname[[i]][[j]])
    PHElistC[[i]][[j]]<-PHEsel$atom
  }
}
processedframe<-0
for (dcdname in dcdlist) {
  print(c(dcdname,Sys.time()))
  mydcd<-read.dcd(file.path(dcddir,dcdname),verbose = FALSE,big=TRUE)
  mydcdcell<-read.dcd(file.path(dcddir,dcdname),verbose = FALSE,cell = TRUE,big=TRUE)
  totframe<-nrow(mydcdcell)
  print(c(totframe,Sys.time()))
  PiIntTemp<-foreach(frame=seq(1,totframe,step)) %dopar% (
    getPi_spatialcutoff(frame)
    )
  if (exists("outname")) {
    load(outname)
    file.remove(outname)
  } else {
    PiAnaOut<-list()
    for (name in names(PiIntTemp[[1]])) {
      PiAnaOut[[name]]<-0
    }
  }
  for (name in names(PiIntTemp[[1]])) {
  PiAnaOut[[name]]<-foreach(frame=1:length(PiIntTemp),.combine="+") %do% (
    PiIntTemp[[frame]][[name]]
  ) +PiAnaOut[[name]]
  }
  processedframe<-processedframe+length(PiIntTemp)
  PiAnaOut[["Frame"]]<-processedframe
  outname<-paste0(outputdir,outnamepre,"_PiInt_",firstrun,"to",regmatches(dcdname, regexpr("\\d+", dcdname)),".rda")
  save(PiAnaOut,file = outname)
  print(c(dcdname,processedframe))
}