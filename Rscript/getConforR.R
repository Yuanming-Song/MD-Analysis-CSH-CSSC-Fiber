#get PHE PHE distance and radial distance of centre of two PHE from fibre centre 
getConforR<-function(frame,Rbreaks=seq(0, BinMaxR, by = BinR)) {
  out<-list()
  tempdcd<-mydcd[frame,] #get coor for current frame 
  tempcell<-mydcdcell[frame,] #get box size
  tempdcd<-pbcwrap_fibre(tempdcd,tempcell) #pbd correction
  particle<-matrix(tempdcd,ncol = 3,byrow = TRUE) #convert coord to a 3 column matrix (x, y, z) 
  out<-t(sapply(seq(1,maxCSSC,1), function(resid) {
    PHEcoor1<-particle[PHEClist[resid,1],] #extract coor for first ring
    PHEcoor2<-particle[PHEClist[resid,2],] #extract coor for first ring
    R1<-sqrt(PHEcoor1[1]^2+PHEcoor1[2]^2)
    R2<-sqrt(PHEcoor2[1]^2+PHEcoor2[2]^2)
    Rlist<-c(R1,R2)
    c(euclidean(PHEcoor1,PHEcoor2),max(Rlist)) #output distance between two rings and average radial position 
    
    #AveCoor<-(PHEcoor1+PHEcoor2)/2 #get average positon
    #AveCoor<-(particle[Ssel[resid,1],]+particle[Ssel[resid,2],])/2 #get average positon
    #c(euclidean(PHEcoor1,PHEcoor2),sqrt(AveCoor[1]^2+AveCoor[2]^2)) #output distance between two rings and average radial position 
  }))
  hist_x <- cut(out[,1], breaks = seq(0, BinMaxD, by = BinD), include.lowest = TRUE)
  hist_y <- cut(out[,2], breaks = Rbreaks, include.lowest = TRUE)
  hist_2d <- as.matrix(table(hist_x, hist_y)) #cut 2 column out data into 2d histogram, range from 0 to 30 for both column.
  return(hist_2d)#return 2d histogram matrix 
}
