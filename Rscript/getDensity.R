source("Rscript/package.R")
source("Rscript/base_fun.R")
outputdir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/mix_27A/analysis/data/"
if (TRUE) {
  dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/mix_27A/"
  mypdb<-read.pdb(file.path(dcddir, "setup/CSH_CSSC_mix_27A_cyl_water_ions.pdb"),hex = TRUE)
  outputpre<-"mix_27A_cylcssc_fibre_density_"
  uniseltextlist<-c("CSSC","CSH","TIP3") #unique resname in sim
} else {
  dcddir<-"/dfs9/tw/yuanmis1/mrsec/cyl1/cssc_27A/"
  mypdb<-read.pdb(file.path(dcddir, "setup/CSSC_pure_27A_cyl_water_ions.pdb"),hex = TRUE)
  uniseltextlist<-c("CSSC","TIP3") #unique resname in sim
  outputpre<-"pure_27A_fibre_density_"
}
dcdprefix<-"nvt"
firstrun<-12
lastrun<-100
step<-10
dcdlist<-getdcdlist(firstrun,midfix="cylcssc.")
for (dcdname in dcdlist) {
  print(c(dcdname,Sys.time())) #print time
  mydcd<-read.dcd(file.path(dcddir,dcdname),verbose = FALSE,big = TRUE) #read dcd file
  mydcdcell<-read.dcd(file.path(dcddir,dcdname),verbose = FALSE,cell = TRUE) #read box size
  totframe<-nrow(mydcdcell) #how many frames
  print(c(totframe,Sys.time())) #print time again
  if (length(seq(1,totframe,step))>1) { #parallel compute with multiple frames
    temprho<-foreach(frame=seq(1,totframe,step),.combine=rbind) %dopar% (
      getrho_fibre(frame)
    )
  } else {
    if (dcdname!=dcdlist[length(dcdlist)]) { #only one frame but not last file
      temprho<-getrho_fibre(1)
    } else {
      break #ignore last file since it's still running 
    }
  }
  rm(mydcd) #clear space for memory
  rm(mydcdcell) #same
  if (length(seq(1,totframe,step))>1) { #treat as nframe x 2 matrix
    temprho<-cbind(temprho,dcdname)
  } else {
    temprho<-c(temprho,dcdname) #single vector with only 1 frame processed 
  }
  if (exists("outname")) {
    load(outname) #add new data to old data
  }
  if (exists("rhoout")) { #with old data present
    if (length(seq(1,totframe,step))>1) {
      temprho[,1]<-step-1+as.numeric(temprho[,1])+max(as.numeric(rhoout[,1])) #treat as nframe x 2 matrix
      temprho<-as.data.frame(temprho)
    } else {
      temprho[1]<-step-1+as.numeric(temprho[1])+max(as.numeric(rhoout[,1])) #single vector with only 1 frame processed 
      temprho <- data.frame(V1 = temprho[1], V2 = temprho[2],V3=temprho[3])
    }
  } else { #first run
    temprho<-as.data.frame(temprho)
  }
  temprho$V1<-as.numeric(temprho$V1) #extra formatting 
  temprho$V2<-as.numeric(temprho$V2) #extra formatting
  temprho <- temprho %>% mutate(Time = V1 * 0.02) #extra formatting in ns
  colnames(temprho)<-c("Frame","rho","Run","Time") #formatting
  if (exists("rhoout")) {
    rhoout<-rbind(rhoout,temprho) #combine with old data
  } else {
    rhoout<-temprho
  }
  if (exists("outname")) {
    file.remove(outname) #delete old data output
  }
  outname<-paste0(outputdir,outputpre,firstrun,"to",regmatches(dcdname, regexpr("\\d+", dcdname)),".rda")
  save(rhoout,file = outname) #clear space 
  rm(rhoout) #clear space
}
