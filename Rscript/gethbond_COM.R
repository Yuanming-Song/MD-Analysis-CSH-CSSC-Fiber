# Set the distance bin size max and min for distance
disbinsize <- 0.1
disbinmax <- 40# sqrt(2 * 100^2) / 2
disbinmin <- 0
# Set the distance bin size max and min for distance
anglebinsize <- 0.05
anglebinmax <- 1
anglebinmin <- -1
outname<-"rawhbond_COM.rda"
step<-10
dcdname<-paste0("mix_27A_sum_every",step,".dcd")
#dcdname<-"test.dcd"
pdbname<-"data/mix_27A.pdb"
### first frame to be considered
startframe<-1
### for pdb & traj file names
dcdname<-paste0("mix_27A_sum_every",step,".dcd")
#dcdname<-"test.dcd"
pdbname<-"data/mix_27A.pdb"
#source base function file
source("Rscript/gethbond_COM.base.R")
#load rda for am COM
load("mix_27A_COM.edgel.stack.rda")
mycell<-read.dcd(dcdname,cell=T,big=TRUE,verbose=FALSE)
mypdb<-read.pdb(pdbname)
mydcd<-read.dcd(dcdname,verbose=FALSE,big=TRUE)
#box<-mean(mycell[,1])*mean(mycell[,2])*mean(mycell[,3])
print(paste(mycell[1,1],mycell[1,2],mycell[1,3]))
nframes<-nrow(mycell)
if (file.exists(outname)) {
  load(outname)
  offsetframe<-length(totrawhbond)
} else {
  offsetframe<-0
}
#parallel analyze each frame
for (batch in 1:ceiling((nframes-offsetframe)/perbatch)) {
  begframe<-(batch-1)*perbatch+1+offsetframe
  if (batch!=ceiling(nframes/perbatch)) {
    endframe<-batch*perbatch+offsetframe
  } else {
    endframe<-nframes
  }
  start_time <- Sys.time()
  print(paste("start",batch,start_time))
  totrawhbondtemp<-foreach(frame=begframe:endframe) %dopar% (
    getrawhbond(frame)
  )
  # Stop tracking time
  end_time <- Sys.time()
  # Calculate the elapsed time
  elapsed_time <- end_time - start_time
  if (!exists("totrawhbond")) {
    totrawhbond<-totrawhbondtemp
  } else {
    totrawhbond<-append(totrawhbond,totrawhbondtemp)
  }
  save(totrawhbond,file=outname)
  print(paste(begframe,endframe,length(totrawhbond),elapsed_time))
  #get hbond function for each
}
