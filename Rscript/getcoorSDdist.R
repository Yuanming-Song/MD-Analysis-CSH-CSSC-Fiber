source("Rscript/package.R")
source("Rscript/base_fun.R")
tauList<-round(10^(seq(2,4,0.5)))
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
if (FALSE) {
  dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/mix_27A/"
  #dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/cssc_27A/"
  mypdb<-read.pdb(file.path(dcddir, "setup/CSH_CSSC_mix_27A_cyl_water_ions.pdb"),hex = TRUE)
  maxCSH<-max(mypdb$atom[which(mypdb$atom$resid=="CSH"),]$resno)
  UniqResnameList<-c("CSH","CSSC")
  maxCSH<-max(mypdb$atom[which(mypdb$atom$resid=="CSH"),]$resno)
  if (FALSE) {
    outnamepre<-"mix_27A_cylcssc_fibre_coor_sd_dist"
    midfix<-"cylcssc."
  } else {
    outnamepre<-"mix_27A_fibre_coor_sd_dist"
  }
} else {
  outnamepre<-"pure_27A_fibre_coor_sd_dist"
  dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/cssc_27A/"
  mypdb<-read.pdb(file.path(dcddir, "setup/CSSC_pure_27A_cyl_water_ions.pdb"),hex = TRUE)
  UniqResnameList<-c("CSSC") #unique resname in sim
  
}
maxCSSC<-max(mypdb$atom[which(mypdb$atom$resid=="CSSC"),]$resno)
fibresel<-atom.select(mypdb,resid=c("CSSC","CSH")) 
dcdlist<-getdcdlist(firstrun,midfix=midfix)
moietylist<- c("AM", "PHE", "THI")
resname<-"CSSC"
mincssc<-1
maxcssc<-maxCSSC
moietyxyz<-list()
for (moiety in UniqResnameList) {
  moietyxyz[[moiety]]<-list()
  for (i in 1:get(paste0("max",moiety))) {
    #####   select input res as input for getdist()
    moietysel<-atom.select(mypdb,"noh",resno=i,resid=resname,operator = "AND")
    #####   store to a list
    moietyxyz[[moiety]][[i]]<-moietysel$xyz
  }
}

processedframe<-0
HistMatrix<-list()
HistMatrix[["OrientBreaks"]]<-OrientBreaks
HistMatrix[["ZBreaks"]]<-seq(-159.5, 159.5,1)
HistMatrix[["RBreaks"]]<-seq(0.5, 29.5,1)
HistMatrix[["AbsOrientBreaks"]]<-AbsOrientBreaks
COM<-list()
Xcoor<-list()
Ycoor<-list()
Zcoor<-list()
for (resname in UniqResnameList) {
  Xcoor[[resname]]<-list()
  Ycoor[[resname]]<-list()
  Zcoor[[resname]]<-list()
  for (resid in 1:get(paste0("max",resname))) {
    Xcoor[[resname]][[resid]]<-list()
    Ycoor[[resname]][[resid]]<-list()
    Zcoor[[resname]][[resid]]<-list()
  }
}
outname<-paste0(outputdir,outnamepre,firstrun,"to",firstrun,".rda")
for (dcdname in dcdlist) {
  print(c(dcdname,Sys.time()))
  mydcd<-read.dcd(file.path(dcddir,dcdname),verbose = FALSE,big=TRUE)
  mydcdcell<-read.dcd(file.path(dcddir,dcdname),verbose = FALSE,cell = TRUE,big=TRUE)
  totframe<-nrow(mydcdcell)
  print(c(totframe,Sys.time()))
  if (length(seq(1,totframe,step))>1) {
    TempCOM<-foreach(frame=seq(1,totframe,step)) %dopar% (
      getCOM_perResID(frame,FALSE)
    )
    COM<-c(COM,TempCOM)
    processedframe<-processedframe+length(seq(1,totframe,step))
    
  } else {
    if (dcdname!=dcdlist[length(dcdlist)]) {
      COM[[length(COM)+1]]<-getCOM_perResID(1,FALSE)
      processedframe<-processedframe+1
      
    } else {
      
      break
    }
  }
}
outnameraw<-paste0(outputdir,outnamepre,firstrun,"to",regmatches(dcdname, regexpr("\\d+", dcdname)),"_raw.rda")
save(COM,file = outnameraw)
#load("data/pure_27A_fibre_coor_sd_dist12to31_raw.rda")
cartlist<-c("x","y","z")
frame_chunks <- split(1:length(COM), sort(rep(1:corenum, length = length(COM))))
for (resname in UniqResnameList) {
  for (cooind in 1:3) {
    out <- foreach(chunk = frame_chunks, .combine = 'rbind') %dopar% {
      temp_out <- NULL
      for (frame in chunk) {
      frame<-chunk[1]  
      temp_out <- rbind(temp_out, t(COM[[frame]][[resname]])[cooind,])

      }
      temp_out
    }
    # out<-c()
    # for (frame in 1:length(COM)) {
    #   out<-rbind(out,t(COM[[frame]][[resname]])[cooind,])
    #   if (frame%%1000==0){
    #     print(frame)
    #   }
    # }
    #save(out,file = "MSDtestCOM.rda")
    write.table(out, 
                file = paste0("MSD_analysis/",
                              outnamepre,firstrun,
                              "to",
                              regmatches(dcdname, regexpr("\\d+", dcdname)),
                              "_",
                              resname,
                              "_",
                              cartlist[cooind],
                              ".dat")
                , col.names = FALSE, row.names = FALSE)
  }
}
#paste0(outputdir,outnamepre,firstrun,"to",regmatches(dcdname, regexpr("\\d+", dcdname)),"_raw.rda")

if (0) {
  MSDraw<-c()
  for (resname in UniqResnameList) {
    Xcoor[[resname]]<-list()
    Ycoor[[resname]]<-list()
    Zcoor[[resname]]<-list()
    for (tau in c(tauList,length(COM))) {
      Xcoor[[resname]][[tau]]<-list()
      Ycoor[[resname]][[tau]]<-list()
      Zcoor[[resname]][[tau]]<-list()
      for (ntau in 1:(floor(length(COM)/tau)+1)) {
        Xcoor[[resname]][[tau]][[ntau]]<-list()
        Ycoor[[resname]][[tau]][[ntau]]<-list()
        Zcoor[[resname]][[tau]][[ntau]]<-list()
        for (resid in 1:get(paste0("max",resname))) {
          Ycoor[[resname]][[tau]][[ntau]][[resid]]<-0
          Xcoor[[resname]][[tau]][[ntau]][[resid]]<-0
          Zcoor[[resname]][[tau]][[ntau]][[resid]]<-0
          for (frame in (ntau-1)*tau+2:(ntau*tau)) {
            if (frame > length(COM)) {
              break
            }
            Xcoor[[resname]][[tau]][[ntau]][[resid]]<-Xcoor[[resname]][[tau]][[ntau]][[resid]]+(COM[[frame]][[resname]][resid,1]-COM[[(ntau-1)*tau+1]][[resname]][resid,1])^2
            Ycoor[[resname]][[tau]][[ntau]][[resid]]<-Ycoor[[resname]][[tau]][[ntau]][[resid]]+(COM[[frame]][[resname]][resid,2]-COM[[(ntau-1)*tau+1]][[resname]][resid,2])^2
            Zcoor[[resname]][[tau]][[ntau]][[resid]]<-Zcoor[[resname]][[tau]][[ntau]][[resid]]+(COM[[frame]][[resname]][resid,3]-COM[[(ntau-1)*tau+1]][[resname]][resid,3])^2
          }
          Xcoor[[resname]][[tau]][[ntau]][[resid]]<-Xcoor[[resname]][[tau]][[ntau]][[resid]]/(tau)
          Ycoor[[resname]][[tau]][[ntau]][[resid]]<-Ycoor[[resname]][[tau]][[ntau]][[resid]]/(tau)
          Zcoor[[resname]][[tau]][[ntau]][[resid]]<-Zcoor[[resname]][[tau]][[ntau]][[resid]]/(tau)
        }
        Xcoor[[resname]][[tau]][[ntau]]<-mean(unlist(Xcoor[[resname]][[tau]][[ntau]]))
        Ycoor[[resname]][[tau]][[ntau]]<-mean(unlist(Ycoor[[resname]][[tau]][[ntau]]))
        Zcoor[[resname]][[tau]][[ntau]]<-mean(unlist(Zcoor[[resname]][[tau]][[ntau]]))
      }
      if (length(1:(floor(length(COM)/tau)))>1) {
        TempMSD<-foreach(ntau=1:(floor(length(COM)/tau))) %dopar% (
          getMSD(ntau,resname)
        ) 
      }else {
        TempMSD<-list()
        TempMSD[[1]]<-getMSD(ntau,resname)
        
      }
      TempMSDx <- sapply(TempMSD, function(x) x[[1]])
      TempMSDy <- sapply(TempMSD, function(x) x[[2]])
      TempMSDz <- sapply(TempMSD, function(x) x[[3]])
      Xcoor[[resname]][[tau]]<-mean(TempMSDx)
      Ycoor[[resname]][[tau]]<-mean(TempMSDy)
      Zcoor[[resname]][[tau]]<-mean(TempMSDz)
      for (coor in c("X","Y","Z")) {
        MSDraw<-rbind(MSDraw,c(resname,tau,coor,get(paste0(coor,"coor"))[[resname]][[tau]]))
      }
    }
  }
  outname<-paste0(outputdir,outnamepre,firstrun,"to",regmatches(dcdname, regexpr("\\d+", dcdname)),".rda")
  save(MSDraw,file=outname)  
  
  
  OutHist<-list()
  for (resname in UniqResnameList) {
    SDx<-c()
    SDy<-c()
    SDz<-c()
    for (resid in 1:get(paste0("max",resname))) {
      SDx<-c(SDx,sd(unlist(Xcoor[[resname]][[resid]])))
      SDy<-c(SDy,sd(unlist(Ycoor[[resname]][[resid]])))
      SDz<-c(SDz,sd(unlist(Zcoor[[resname]][[resid]])))
    }
    OutHist[[resname]]<-list()
    OutHist[[resname]][["x"]]<-cbind(seq(0,10-0.1,0.1),table(cut(SDx,breaks=seq(0,10,0.1))))
    OutHist[[resname]][["y"]]<-cbind(seq(0,10-0.1,0.1),table(cut(SDy,breaks=seq(0,10,0.1))))
    OutHist[[resname]][["z"]]<-cbind(seq(0,70-0.1,0.1),table(cut(SDz,breaks=seq(0,70,0.1))))
    
  }
  if (file.exists(outname)) {
    file.remove(outname)
  }
  outname<-paste0(outputdir,outnamepre,firstrun,"to",regmatches(dcdname, regexpr("\\d+", dcdname)),".rda")
  save(OutHist,file=outname)  
  print(c(dcdname,processedframe))
}
