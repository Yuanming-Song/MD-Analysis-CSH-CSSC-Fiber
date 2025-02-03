source("Rscript/package.R")
source("Rscript/base_fun.R")
abs<-TRUE
if (FALSE) {
  outnamepre<-"mix_27A_cylcssc_fibre_PHEdis_in_out"
  dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/mix_27A/"
  #dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/cssc_27A/"
  outputdir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/mix_27A/analysis/data/"
  mypdb<-read.pdb(file.path(dcddir, "setup/CSH_CSSC_mix_27A_cyl_water_ions.pdb"),hex = TRUE)
  UniqResnameList<-c("CSH","CSSC")
  maxCSH<-max(mypdb$atom[which(mypdb$atom$resid=="CSH"),]$resno)
} else {
  dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/cssc_27A/"
  mypdb<-read.pdb(file.path(dcddir, "setup/CSSC_pure_27A_cyl_water_ions.pdb"),hex = TRUE)
  UniqResnameList<-c("CSSC") #unique resname in sim
  outputpre<-"pure_27A_fibre_hbond_spati"
}
maxCSSC<-max(mypdb$atom[which(mypdb$atom$resid=="CSSC"),]$resno)

dcdprefix<-"nvt"
firstrun<-5
lastrun<-100
step<-1
dcdlist<-getdcdlist(firstrun,midfix="")

#Determine outer or inner cyl
FibreCoorPercentCutZ<-0.9
FibreRangeR<-17
PHElistC<-c("C5", "C16")
PHElistSelC<-atom.select(mypdb,resid="CSSC",elety=PHElistC)
PHEClist<-matrix(PHElistSelC$atom,ncol = 2,byrow = TRUE)
PHEdisBreak<-seq(1,29,1)
PHEdisBreakBin<-seq(0.5,29.5,1)
PHEdis<-list()
for (dcdname in dcdlist) {
  print(c(dcdname,Sys.time()))
  mydcd<-read.dcd(file.path(dcddir,dcdname),verbose = FALSE,big=TRUE)
  mydcdcell<-read.dcd(file.path(dcddir,dcdname),verbose = FALSE,cell = TRUE,big=TRUE)
  totframe<-nrow(mydcdcell)
  print(c(totframe,Sys.time()))
  if (length(seq(1,totframe,step))>1) {
    PHEdisTemp<-foreach(frame=1:totframe) %dopar% (getConfor(frame))
  } else {
    if (dcdname!=dcdlist[length(dcdlist)]) {
      PHEdisTemp<-list()
      PHEdisTemp[[1]]<-getConfor(1)
    } else {
      break
    }
  }
  PHEdis<-c(PHEdis,PHEdisTemp)
  PHEdisout<-c()
  PHEdis_sum<-list()
  for (Resname in unique(names(PHEdis[[1]]))) {
    for (Position in unique(names(PHEdis[[1]][[Resname]]))) {
      PHEdis_sum[[Resname]][[Position]]<-0
      for (frame in 1:length(PHEdis)) {
        PHEdis_sum[[Resname]][[Position]]<-PHEdis_sum[[Resname]][[Position]]+PHEdis[[frame]][[Resname]][[Position]]
      }
      Freq<-as.numeric(PHEdis_sum[[Resname]][[Position]])
      Freq_norm<-Freq/sum(Freq)
      Type<-paste(Resname,Position,sep="_")
      tempout<-as.data.frame(cbind(cbind(cbind(PHEdisBreakBin,Freq),Freq_norm),Type))
      PHEdisout<-rbind(PHEdisout,tempout)
      
    }
  }
  if (exists("outname")) {
    if (file.exists(outname)) {
      file.remove(outname)
    }
  }
  outname<-paste0(outputdir,outnamepre,firstrun,"to",regmatches(dcdname, regexpr("\\d+", dcdname)),".rda")
  PHEdisout<-as.data.frame(PHEdisout)
  for (i in 1:3) {
    PHEdisout[,i]<-as.numeric(PHEdisout[,i])
  }
  save(PHEdisout,file = outname)
  rm(PHEdisout)
  rm(PHEdisTemp)
}