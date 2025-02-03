#get distribution of coordinate
getcoordist_fibre<-function(frame) { 
  tempdcd<-mydcd[frame,] #get coor for current frame 
  tempcell<-mydcdcell[frame,] #get box size
  tempdcd<-pbcwrap_fibre(tempdcd,tempcell) #wrap coordinates
  zhisout<-c() #empty output name 
  rhisout<-c()
  for (seltext in uniseltextlist) {
    fibresel<-atom.select(mypdb,"noh",resid=seltext) #heavy atoms only
    fibrecoor<-tempdcd[fibresel$xyz] #get coor
    FibreCoorX<-fibrecoor[seq(1,length(fibrecoor)-2,3)] #x
    FibreCoorY<-fibrecoor[seq(2,length(fibrecoor)-1,3)] #y
    fibrecoorR<-sqrt(FibreCoorX^2+FibreCoorY^2) #radial distance
    FibreCoorZ<-fibrecoor[seq(3,length(fibrecoor)-0,3)] #z 
    
    FibreCoorZfib<-FibreCoorZ[which(fibrecoorR<27)]  #in case only need to be within constraint 
    fibrecoorRfib<-fibrecoorR[which(abs(FibreCoorZ)<100)] #only 
    
    coorzhis<-table(cut(FibreCoorZfib,breaks = seq(-180,180,1))) #get raw dist
    coorrhis<-table(cut(fibrecoorRfib,breaks = seq(0,71,1))) #get raw dist
    zhisout<-cbind(zhisout,coorzhis) #combine output 
    rhisout<-cbind(rhisout,coorrhis) #combine output 
  }
  out<-list() 
  out[[1]]<-zhisout
  out[[2]]<-rhisout
  out
}
