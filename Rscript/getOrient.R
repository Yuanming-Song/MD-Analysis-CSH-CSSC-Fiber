source("Rscript/package.R")
source("Rscript/base_fun.R")
abs<-FALSE

outputdir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/mix_27A/analysis/data/"
midfix<-""
if (FALSE) {
  dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/mix_27A/"
  mypdb<-read.pdb(file.path(dcddir, "setup/CSH_CSSC_mix_27A_cyl_water_ions.pdb"),hex = TRUE)
  outnamepre<-"mix_27A_fibre_orient_"
  uniseltextlist<-c("CSSC","CSH","TIP3") #unique resname in sim
  if (FALSE) {
    midfix<-"cylcssc"
    outnamepre<-paste0(outnamepre,midfix,"_")
    midfix<-paste0(midfix,".")
  }
  maxCSH<-max(mypdb$atom[which(mypdb$atom$resid=="CSH"),]$resno)
} else {
  dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/cssc_27A/"
  mypdb<-read.pdb(file.path(dcddir, "setup/CSSC_pure_27A_cyl_water_ions.pdb"),hex = TRUE)
  uniseltextlist<-c("CSSC","TIP3") #unique resname in sim
  outnamepre<-"pure_27A_fibre_orient_"
}

dcdprefix<-"nvt"
firstrun<-12
lastrun<-100
step<-1
dcdlist<-getdcdlist(firstrun,midfix = midfix)
mypdb$atom[which(mypdb$atom$resid=="CSSC"),]$resno <- with(mypdb$atom[which(mypdb$atom$resid=="CSSC"),], rep(1:ceiling(max(eleno)/26), each = 26)[1:length(eleno)])
maxCSSC<-max(mypdb$atom[which(mypdb$atom$resid=="CSSC"),]$resno)

#Determine outer or inner cyl
FibreCoorPercentCutZ<-0.9
FibreRangeR<-17
OrientSelText<-list()
OrientSelText[["CSH"]][[1]]<-c("C2","C8")
OrientSelText[["CSH"]][[2]]<-c("C2","C8")
OrientSelText[["CSSC"]][[1]]<-c("C2","C8")
OrientSelText[["CSSC"]][[2]]<-c("C15","C14")
UniqResnameList<-c("CSH","CSSC")
namelist<-c("C2","C8","C15","C14")
NameIndexList<-list()
for (i in 1:4) {
  name<-namelist[i]
  tempsel<-atom.select(mypdb,elety=name)
  NameIndexList[[i]]<-tempsel$atom
}

OrientBreaks<-seq(-0.95,1.05,0.1)
AbsOrientBreaks<-seq(0.05,1.05,0.1)

ProcessedFrame<-0
for (dcdname in dcdlist) {
  print(c(dcdname,Sys.time()))
  mydcd<-read.dcd(file.path(dcddir,dcdname),verbose = FALSE,big=TRUE)
  mydcdcell<-read.dcd(file.path(dcddir,dcdname),verbose = FALSE,cell = TRUE,big=TRUE)
  totframe<-nrow(mydcdcell)
  print(c(totframe,Sys.time()))
  if (length(seq(1,totframe,step))>1) {
    Orienttemp<-foreach(frame=seq(1,totframe,step)) %dopar% (getOrient(frame))
  } else {
    if (dcdname!=dcdlist[length(dcdlist)]) {
      Orienttemp<-list()
      Orienttemp[[1]]<-getOrient(1)
    } else {
      break
    }
  }
  ProcessedFrame<-length(seq(1,totframe,step))+ProcessedFrame
  if (exists("outname")) {
    if (file.exists(outname)) {
      load(outname)
      file.remove(outname)
    }
  } else {
    OrientDist<-list()
    for (name in names(Orienttemp[[1]])) {
      OrientDist[[name]]<-0
    }
    OrientDist[["OrientBreaks"]]<-OrientBreaks
    OrientDist[["ZBreaks"]]<-seq(-159.5, 159.5,1)
    OrientDist[["RBreaks"]]<-seq(0.5, 29.5,1)
    OrientDist[["AbsOrientBreaks"]]<-AbsOrientBreaks
  }
  
  for (name in names(Orienttemp[[1]])) {
    OrientDist[[name]]<-0
    for (frame in 1:length(Orienttemp)) {
      OrientDist[[name]]<-OrientDist[[name]]+Orienttemp[[frame]][[name]]
    }
  }
  OrientDist[["frame"]]<-ProcessedFrame
  if (abs==FALSE) {
    outname<-paste0(outputdir,outnamepre,firstrun,"to",regmatches(dcdname, regexpr("\\d+", dcdname)),".rda")
  } else {
    outname<-paste0(outputdir,outnamepre,"_abs_",firstrun,"to",regmatches(dcdname, regexpr("\\d+", dcdname)),".rda")
  }
  save(OrientDist,file = outname)
  rm(OrientDist)
  rm(Orienttemp)
}
if (0) {
  Orient<-c(Orient,Orienttemp)
  Orientout<-c()
  Orient_sum<-list()
  for (Resname in unique(names(Orient[[1]]))) {
    for (Position in unique(names(Orient[[1]][[Resname]]))) {
      Orient_sum[[Resname]][[Position]]<-0
      for (frame in 1:length(Orient)) {
        Orient_sum[[Resname]][[Position]]<-Orient_sum[[Resname]][[Position]]+Orient[[frame]][[Resname]][[Position]]
      }
      Freq<-as.numeric(Orient_sum[[Resname]][[Position]])
      Freq_norm<-Freq/sum(Freq)
      Type<-paste(Resname,Position,sep="_")
      tempout<-as.data.frame(cbind(cbind(cbind(OrientBreaks,Freq),Freq_norm),Type))
      Orientout<-rbind(Orientout,tempout)
    }
  }
  
}