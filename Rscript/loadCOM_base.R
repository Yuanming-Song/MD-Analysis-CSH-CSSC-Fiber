if (file.exists(comoutname1)) {
  load(comoutname1)
} else {
  mypdb<-read.pdb(pdbname)
  mydcd<-read.dcd(dcdname,verbose=FALSE,big=TRUE)
  COMlistPHE<-foreach (framenum=startframe:nframes) %dopar% {
    getCOM(mincssc,maxcsh,framenum,phemoietyname,phemoietyname2) 
  }  
  save(COMlistPHE,file = comoutname1)
}
if (file.exists(comoutname2)) {
  load(comoutname2)
} else {
  mypdb<-read.pdb(pdbname)
  mydcd<-read.dcd(dcdname,verbose=FALSE,big=TRUE)
  COMlistAM<-foreach (framenum=startframe:nframes) %dopar% {
    getCOM(mincssc,maxcsh,framenum,ammoietyname,ammoietyname2) 
  }  
  save(COMlistAM,file = comoutname2)
}
if (file.exists(comoutname3)) {
  load(comoutname3)
} else {
  mypdb<-read.pdb(pdbname)
  mydcd<-read.dcd(dcdname,verbose=FALSE,big=TRUE)
  COMlistTHI<-foreach (framenum=startframe:nframes) %dopar% {
    getCOM(mincssc,maxcsh,framenum,thimoietyname,thimoietyname2) 
  }  
  save(COMlistTHI,file = comoutname3)
}