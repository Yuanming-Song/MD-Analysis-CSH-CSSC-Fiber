source("Rscript/package.R")
source("Rscript/base_fun.R")
abs<-TRUE
BinR<-1
BinMaxR<-30
BinD<-0.5
BinMaxD<-30
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
  UniqResnameList<-c("CSH","CSSC")
  maxCSH<-max(mypdb$atom[which(mypdb$atom$resid=="CSH"),]$resno)
  if (FALSE) {
    outnamepre<-"mix_27A_cylcssc_fibre_PHEdis_R+_"
    midfix<-"cylcssc."
  } else {
    outnamepre<-"mix_27A_fibre_PHEdis_SS"
  }
} else {
  outnamepre<-"pure_27A_fibre_PHEdis_SS"
  dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/cssc_27A/"
  mypdb<-read.pdb(file.path(dcddir, "setup/CSSC_pure_27A_cyl_water_ions.pdb"),hex = TRUE)
  UniqResnameList<-c("CSSC") #unique resname in sim
}
maxCSSC<-max(mypdb$atom[which(mypdb$atom$resid=="CSSC"),]$resno)
fibresel<-atom.select(mypdb,resid=c("CSSC","CSH")) 
Ssel<-atom.select(mypdb,elety=c("S1","S2")) 
Ssel<-matrix(Ssel$atom,ncol = 2,byrow = TRUE) #convert coord to a 2 column matrix 
dcdlist<-getdcdlist(firstrun,midfix=midfix)
#Determine outer or inner cyl
FibreCoorPercentCutZ<-0.9
FibreRangeR<-17
PHElistC<-c("C5", "C16")
PHElistSelC<-atom.select(mypdb,resid="CSSC",elety=PHElistC)
PHEClist<-matrix(PHElistSelC$atom,ncol = 2,byrow = TRUE)
HistMatrix<-0
processedframe<-0
for (dcdname in dcdlist) {
  print(c(dcdname,Sys.time()))
  mydcd<-read.dcd(file.path(dcddir,dcdname),verbose = FALSE,big=TRUE)
  mydcdcell<-read.dcd(file.path(dcddir,dcdname),verbose = FALSE,cell = TRUE,big=TRUE)
  totframe<-nrow(mydcdcell)
  print(c(totframe,Sys.time()))
  if (length(seq(1,totframe,step))>1) {
    HistMatrixTemp<-foreach(frame=1:totframe,.combine = "+") %dopar% (
      getConforR(frame)
    )
    processedframe<-processedframe+length(seq(1,totframe,step))
    
  } else {
    if (dcdname!=dcdlist[length(dcdlist)]) {
      HistMatrixTemp<-getConforR(1)#run getConforR for first frame 
      processedframe<-processedframe+1
    } else {
      break
    }
  }
  
  HistMatrix<-HistMatrix+HistMatrixTemp
  HistMatrixNormR <- apply(HistMatrix, 2, function(col) {
    if (sum(col) != 0) {
      col / sum(col)
    } else {
      col
    }
  })
  colnames(HistMatrixNormR)<-seq(0+BinR/2,BinMaxR-BinR/2,BinR)
  rownames(HistMatrixNormR)<-seq(0+BinD/2,BinMaxD-BinD/2,BinD)
  HistMatrixNormROut <- as.data.frame(HistMatrixNormR)
  HistMatrixNormROut$rowname <- as.numeric(rownames(HistMatrixNormROut))
  HistMatrixNormROut_long <- gather(HistMatrixNormROut, key = "column", value = "value", -rowname)
  colnames(HistMatrixNormROut_long)<-c("D","R","Pr")
  HistMatrixOut <- as.data.frame.matrix(HistMatrix)
  colnames(HistMatrixOut)<-seq(0.5,29.5,1)
  rownames(HistMatrixOut)<-seq(0.25,29.75,0.5)
  HistMatrixOut$rowname <- as.numeric(rownames(HistMatrixOut))
  HistMatrixROut_long <- gather(HistMatrixOut, key = "column", value = "value", -rowname)
  colnames(HistMatrixROut_long)<-c("D","R","Freq")
  for (i in 1:3) {
    HistMatrixROut_long[,i]<-as.numeric(HistMatrixROut_long[,i])
  }
  HistMatrixROut_long<-HistMatrixROut_long %>% mutate(Pr=Freq/(processedframe*220*pi*((R+0.5)^2-(R-0.5)^2)))
  if (exists("outname")) {
      file.remove(outname)
  }
  PheDisOut<-list()
  PheDisOut[["HistMatrixNormROut_long"]]<-HistMatrixNormROut_long
  PheDisOut[["HistMatrixNormROut"]]<-HistMatrixNormROut
  PheDisOut[["HistMatrix"]]<-HistMatrix
  PheDisOut[["HistMatrixROut_long"]]<-HistMatrixROut_long
  PheDisOut[["Frame"]]<-processedframe
  outname<-paste0(outputdir,outnamepre,firstrun,"to",regmatches(dcdname, regexpr("\\d+", dcdname)),".rda")
  save(PheDisOut,file=outname)
  print(c(dcdname,processedframe))
}
if (0) {
  ggplot(HistMatrixNormROut_long, aes(x = as.numeric(R), y = as.numeric(D), fill = Pr)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(title = "2D Histogram of Normalized Matrix",
         x = "R",
         y = "D") +
    theme_minimal() +
    scale_x_continuous(breaks = seq(0, 30, 5)) +
    scale_y_continuous(breaks = seq(0, max(result_df_long$rowname), 5))
}