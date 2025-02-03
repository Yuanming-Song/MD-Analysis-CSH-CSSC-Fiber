
source("Rscript/package.R")
source("Rscript/base_fun.R")
source("Rscript/getcoordist_fibre_base.R")
outputdir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/mix_27A/analysis/data/"
midfix<-""
if (TRUE) {
  dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/mix_27A/"
  mypdb<-read.pdb(file.path(dcddir, "setup/CSH_CSSC_mix_27A_cyl_water_ions.pdb"),hex = TRUE)
  outnamepre<-"mix_27A_fibre_coor_"
  uniseltextlist<-c("CSSC","CSH","TIP3") #unique resname in sim
  if (TRUE) {
    midfix<-"cylcssc"
    outnamepre<-paste0(outnamepre,midfix,"_")
    midfix<-paste0(midfix,".")
  }
} else {
  dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/cssc_27A/"
  mypdb<-read.pdb(file.path(dcddir, "setup/CSSC_pure_27A_cyl_water_ions.pdb"),hex = TRUE)
  uniseltextlist<-c("CSSC","TIP3") #unique resname in sim
  outnamepre<-"pure_27A_fibre_coor_"
}
print(outnamepre)
dcdprefix<-"nvt"
firstrun<-12
lastrun<-100
step<-1
dcdlist<-getdcdlist(firstrun,midfix = midfix)

rhistot<-0 #empty matrix
zhistot<-0 #empty matrix
fibdist<-list()
processedframe<-0
for (dcdname in dcdlist) {
  print(c(dcdname,Sys.time())) #print time
  mydcd<-read.dcd(file.path(dcddir,dcdname),verbose = FALSE,big = TRUE) #read dcd file
  mydcdcell<-read.dcd(file.path(dcddir,dcdname),verbose = FALSE,cell = TRUE) #read box size
  totframe<-nrow(mydcdcell) #how many frames
  print(c(totframe,Sys.time())) #print time again
  if (length(seq(1,totframe,step))>1) { #parallel compute with multiple frames
    coordist<-foreach(frame=seq(1,totframe,step)) %dopar% (
      getcoordist_fibre(frame)
    )
  } else {    
    if (dcdname!=dcdlist[length(dcdlist)]) { #only one frame but not last file
      coordist<-list()
      coordist[[1]]<-getcoordist_fibre(1)  
    } else {
      break #ignore last file since it's still running 
    }
  }
  #recombine results by adding histgram 
  for(frame in 1:length(coordist)) {
    rhistot<-rhistot+coordist[[frame]][[2]]
    zhistot<-zhistot+coordist[[frame]][[1]]
  }
  rm(mydcd) #clear space
  rm(mydcdcell) #clear space
  rhistottemp<-c() #recombine to 3 col data for r
  #loop through each unique resname
  for (i in 1:length(uniseltextlist)) {
    resname<-uniseltextlist[i]
    rhistottemp<-rbind(rhistottemp,cbind(cbind(seq(0.5,71-0.5,1),rhistot[,i]),resname))
  }
  rhistottemp<-as.data.frame(rhistottemp)
  colnames(rhistottemp)<-c("r","Freq","Resname")
  zhistottemp<-c() #recombine to 3 col data for r
  for (i in 1:length(uniseltextlist)) {
    resname<-uniseltextlist[i]
    zhistottemp<-rbind(zhistottemp,cbind(cbind(seq(-180+0.5,180-0.5,1),zhistot[,i]),resname))
  }
  zhistottemp<-as.data.frame(zhistottemp)
  colnames(zhistottemp)<-c("z","Freq","Resname")
  fibdist[[1]]<-rhistottemp
  fibdist[[2]]<-zhistottemp
  if (exists("outname")) {
    if (file.exists(outname)) {
      file.remove(outname) #delete old data file
    }
  }
  #add up processed frames
  processedframe<-processedframe+length(coordist)
  fibdist[[3]]<-processedframe
  outname<-paste0(outputdir,outnamepre,firstrun,"to",regmatches(dcdname, regexpr("\\d+", dcdname)),".rda")
  save(fibdist,file = outname)
  print(outname)
  
}



