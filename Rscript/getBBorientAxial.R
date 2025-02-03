source("Rscript/package.R")
source("Rscript/base_fun.R")
abs<-TRUE
Bin<-0.1
OrientBreaks<-seq(-0.95,1.05,0.1)
AbsOrientBreaks<-seq(0.05,1.05,0.1)
dcdprefix<-"nvt"
firstrun<-12
lastrun<-100
step<-1
outputdir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/mix_27A/analysis/data/"
midfix<-""
if (TRUE) {
  dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/mix_27A/"
  #dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/cssc_27A/"
  mypdb<-read.pdb(file.path(dcddir, "setup/CSH_CSSC_mix_27A_cyl_water_ions.pdb"),hex = TRUE)
  UniqResnameList<-c("CSH","CSSC")
  maxCSH<-max(mypdb$atom[which(mypdb$atom$resid=="CSH"),]$resno)
  if (FALSE) {
    outnamepre<-"mix_27A_cylcssc_fibre_BB_orient_axial"
    midfix<-"cylcssc."
  } else {
    outnamepre<-"mix_27A_fibre_BB_orient_axial"
  }
} else {
  outnamepre<-"pure_27A_fibre_BB_orient_axial"
  dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/cssc_27A/"
  mypdb<-read.pdb(file.path(dcddir, "setup/CSSC_pure_27A_cyl_water_ions.pdb"),hex = TRUE)
  UniqResnameList<-c("CSSC") #unique resname in sim
}
maxCSSC<-max(mypdb$atom[which(mypdb$atom$resid=="CSSC"),]$resno)
fibresel<-atom.select(mypdb,resid=c("CSSC","CSH")) 
dcdlist<-getdcdlist(firstrun,midfix=midfix)
moietylist<- c("AM", "PHE", "THI")
load("data/COMcutoff.rda")
resname<-"CSSC"
mincssc<-1
maxcssc<-maxCSSC
moietyxyz<-list()
for (moiety in c("am","phe","thi")) {
  moietyxyz[[moiety]]<-list()
  moietyxyz[[moiety]][[1]]<-list()
  for (i in 1:maxCSSC) {
    #####   select input res as input for getdist()
    moietysel<-atom.select(mypdb,resno=i,resid=resname,elety=get(paste0(moiety,"moietyname")),operator = "AND")
    #####   store to a list
    moietyxyz[[moiety]][[1]][[i]]<-moietysel$xyz
  }
  ## repeat
  moietyxyz[[moiety]][[2]]<-list()
  for (i in 1:maxCSSC) {
    #####   select input res as input for getdist()
    moietysel<-atom.select(mypdb,resno=i,resid=resname,elety=get(paste0(moiety,"moietyname2")),operator = "AND")
    #####   store to a list
    moietyxyz[[moiety]][[2]][[i]]<-moietysel$xyz
  }
}

BBlist<-c("C1","C20")
BBxyz<-list()
for (BB in BBlist) {
  BBxyz[[BB]]<-list()
  for (i in 1:maxCSSC) {
    BBsel<-atom.select(mypdb,resno=i,resid=resname,elety=BB,operator = "AND")
    BBxyz[[BB]][[i]]<-BBsel$xyz
  }
}
processedframe<-0
HistMatrix<-list()
HistMatrix[["OrientBreaks"]]<-OrientBreaks
HistMatrix[["ZBreaks"]]<-seq(-159.5, 159.5,1)
HistMatrix[["RBreaks"]]<-seq(0.5, 29.5,1)
HistMatrix[["AbsOrientBreaks"]]<-AbsOrientBreaks
for (dcdname in dcdlist) {
  print(c(dcdname,Sys.time()))
  mydcd<-read.dcd(file.path(dcddir,dcdname),verbose = FALSE,big=TRUE)
  mydcdcell<-read.dcd(file.path(dcddir,dcdname),verbose = FALSE,cell = TRUE,big=TRUE)
  totframe<-nrow(mydcdcell)
  print(c(totframe,Sys.time()))
  if (length(seq(1,totframe,step))>1) {
    HistMatrixTemp<-foreach(frame=seq(1,totframe,step)) %dopar% (
      getBBaxialAngHist(frame)
    )
    
    processedframe<-processedframe+length(seq(1,totframe,step))
    
  } else {
    if (dcdname!=dcdlist[length(dcdlist)]) {
      HistMatrixTemp<-list()
      HistMatrixTemp[[1]]<-getBBaxialAngHist(1)#run getConforR for first frame 
      processedframe<-processedframe+1
    } else {
      
      break
    }
  }
  HistMatrixTempsum<-list()
  for (HisName in names(HistMatrixTemp[[1]])) {
    HistMatrixTempsum[[HisName]]<-0
    for (i in length(HistMatrixTemp)) {
      HistMatrixTempsum[[HisName]]<-HistMatrixTempsum[[HisName]]+HistMatrixTemp[[i]][[HisName]]
    }
    if (HisName == "OrdPar") {
      if (!exists("outname")) {
        HistMatrix[[HisName]]<-list()
        HistMatrix[[HisName]][[dcdname]]<-HistMatrixTempsum[[HisName]]
        
      } else {
        HistMatrix[[HisName]][[dcdname]]<-HistMatrixTempsum[[HisName]]
      }
    } else {
      if (!exists("outname")) {
        HistMatrix[[HisName]]<-HistMatrixTempsum[[HisName]]
      } else {
        HistMatrix[[HisName]]<-HistMatrix[[HisName]]+HistMatrixTempsum[[HisName]]
      }
    }
  }
  #  HistMatrixOut<-cbind(seq(-1+Bin/2,1+Bin/2,Bin),HistMatrix,HistMatrix/sum(HistMatrix))
  HistMatrix[["frame"]]<-processedframe
  if (exists("outname")) {
    file.remove(outname)
  }
  outname<-paste0(outputdir,outnamepre,firstrun,"to",regmatches(dcdname, regexpr("\\d+", dcdname)),".rda")
  save(HistMatrix,file=outname)  
  print(c(dcdname,processedframe))
  
}


