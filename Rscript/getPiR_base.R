#get donor acceptor distance, and donor hydrogen acceptor angle based on spatial cutoff 
getPi_spatialcutoff_old<-function(frame) {       
  tempdcd<-mydcd[frame,] #get coor for current frame 
  tempcell<-mydcdcell[frame,] #get box size
  tempdcd<-pbcwrap_fibre(tempdcd,tempcell) #wrap coordinates
  particle<-matrix(tempdcd[fibresel$xyz],ncol = 3,byrow = TRUE) #convert coord to a 3 column matrix (x, y, z) 
  do.call(rbind,sapply (seq(1,length(PHElistC[[1]][[1]]),1),function(Resid) {
    if (Resid < maxCSSC) {
      maxI<-2
    } else {
      maxI<-1
    }
    do.call(rbind,sapply (1:maxI,function(I) {
      PHEC1COM<-particle[PHElistC[[I]][[1]][Resid],] # Get PHE C coordinate
      PHEC2COM<-particle[PHElistC[[I]][[2]][Resid],] # Get PHE C coordinate
      PHEC3COM<-particle[PHElistC[[I]][[3]][Resid],] # Get PHE C coordinate
      PHEnormRef<-getnorm(PHEC1COM,PHEC2COM,PHEC3COM)
      PHErefCOM<-(PHEC1COM+PHEC2COM+PHEC3COM)/3 #get COM of PHE
      CloseIndex <- which(apply(particle, 1, function(row) all(abs(row - PHErefCOM) <= 6))) #screen out atoms that are too far
      ExcIndex1<-PHElistC[[1]][[1]][-seq(1,Resid,1)] #aviod dealing with lower resid 
      CloseIndex1<-which(PHElistC[[1]][[1]]%in% CloseIndex[CloseIndex %in% ExcIndex1]) #select only carbons that are close enough
      ExcIndex2<-PHElistC[[2]][[1]][-seq(1,Resid,1)] #aviod dealing with lower resid 
      CloseIndex2<-which(PHElistC[[2]][[1]]%in% CloseIndex[CloseIndex %in% ExcIndex2]) #select only carbons that are close enough
      CloseResid<-list() #combine comparison 
      CloseResid[[1]]<-CloseIndex1
      CloseResid[[2]]<-CloseIndex2
      do.call(rbind,sapply (1:2,function(i) {
        do.call(rbind,sapply (CloseResid[[i]],function(CompResid) {
          CompPHEC1COM<-particle[PHElistC[[i]][[1]][CompResid],] # Get PHE C coordinate
          CompPHEC2COM<-particle[PHElistC[[i]][[2]][CompResid],] # Get PHE C coordinate
          CompPHEC3COM<-particle[PHElistC[[i]][[3]][CompResid],]
          PHECompCOM<-(CompPHEC1COM+CompPHEC2COM+CompPHEC3COM)/3 #get COM of PHE
          DisPHE<-euclidean(PHErefCOM,PHECompCOM) #get distance between PHE
          PHEAveCOM<-(PHErefCOM+PHECompCOM)/2 #get COM of 2 PHE
          PHEAveCOMR<-sqrt(sum((PHEAveCOM^2)[c(1,2)])) #get R of two PHE
          CompPHEvec1<-CompPHEC2COM-CompPHEC1COM #get C-C vector in phe
          CompPHEvec1<-CompPHEC3COM-CompPHEC1COM #get C-C vector in phe
          PHEang<-getcostheta(PHEnormRef,getnorm(CompPHEC1COM,CompPHEC2COM,CompPHEC3COM)) #get cos of angle between normal vectors
          acosPHEang<-acos(PHEang)*180/3.14 #get angle in rad
          c(DisPHE
            ,PHEang
            ,acosPHEang
            ,paste(Resid,I,CompResid,i,sep="-")
            # ,hydrogenindex-1
            # ,AcceptorClose[o]-1
            # ,frame
            ,PHEAveCOM
            ,PHEAveCOMR
            ,processedframe*step+frame
          )
        }, simplify = FALSE))
      }, simplify = FALSE))
    }, simplify = FALSE))
  }, simplify = FALSE))
}

getPi_spatialcutoff<-function(frame) {       
  tempdcd<-mydcd[frame,] #get coor for current frame 
  tempcell<-mydcdcell[frame,] #get box size
  tempdcd<-pbcwrap_fibre(tempdcd,tempcell) #wrap coordinates
  particle<-matrix(tempdcd[fibresel$xyz],ncol = 3,byrow = TRUE) #convert coord to a 3 column matrix (x, y, z) 
  PHECOM<-list()
  PHECOMnorm<-list()
  for (i in 1:2) {
    PHEC1COM<-particle[PHElistC[[i]][[1]],]
    PHEC2COM<-particle[PHElistC[[i]][[2]],]
    PHEC3COM<-particle[PHElistC[[i]][[3]],]
    PHECOM[[i]]<-(PHEC1COM+PHEC2COM+PHEC3COM)/3
    PHECOMnorm[[i]]<-getnorm(PHEC1COM,PHEC2COM,PHEC3COM)
  }
  PiInt<-do.call(rbind,sapply (seq(1,length(PHElistC[[1]][[1]]),1),function(Resid) {
    if (Resid < maxCSSC) {
      maxI<-2
    } else {
      maxI<-1
    }
    do.call(rbind,sapply (1:maxI,function(I) {
      PHErefCOM<-PHECOM[[I]][Resid,]
      PHEnormRef<-PHECOMnorm[[I]][Resid,]
      if (sqrt(PHErefCOM[1]^2+PHErefCOM[2]^2)<28) {
        do.call(rbind,sapply (1:I,function(i) {
          CloseResid <- which(apply(PHECOM[[i]], 1, function(row) all(abs(row - PHErefCOM) <= 7))) #screen out atoms that are too far
          CloseResid<-CloseResid[which(CloseResid>Resid)]
          do.call(rbind,sapply (CloseResid,function(CompResid) {
            PHECompCOM<-PHECOM[[i]][CompResid,] #get COM of PHE
            DisPHE<-euclidean(PHErefCOM,PHECompCOM) #get distance between PHE
            PHEAveCOM<-(PHErefCOM+PHECompCOM)/2 #get COM of 2 PHE
            PHEAveCOMR<-sqrt(sum((PHEAveCOM^2)[c(1,2)])) #get R of two PHE
            PHEcompRef<-PHECOMnorm[[i]][CompResid,]
            PHEang<-getcostheta(PHEnormRef,PHEcompRef) #get cos of angle between normal vectors
            acosPHEang<-acos(PHEang)*180/3.14 #get angle in rad
            c(DisPHE
              ,PHEang
              ,acosPHEang
              ,paste(Resid,I,CompResid,i,sep="-")
              # ,hydrogenindex-1
              # ,AcceptorClose[o]-1
              # ,frame
              ,PHEAveCOM[3]
              ,PHEAveCOMR
              ,processedframe*step+frame
            )
          }, simplify = FALSE))
        }, simplify = FALSE))
      }
    }, simplify = FALSE))
  }, simplify = FALSE))
  out<-list()
  hist_d <- cut(as.numeric(PiInt[,1]), breaks = seq(0, 7, by = 0.2), include.lowest = TRUE)
  hist_ang <- cut(abs(as.numeric(PiInt[,2])), breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
  PiHist_d_ang_raw <- as.matrix(table(hist_d, hist_ang)) #cut 2 column out data into 2d histogram, range from 0 to 30 for both column.
  t_shape_index<-which(as.numeric(PiInt[,1]) < tdish & abs(as.numeric(PiInt[,2])) <= tang  & as.numeric(PiInt[,1]) > tdisl)
  parallel_count_index <- which(as.numeric(PiInt[,1]) < pdis & abs(as.numeric(PiInt[,2])) >= pang)
  t_shape_count<-length(t_shape_index)
  parallel_count <- length(parallel_count_index)
  out[["PiHist_d_ang_raw"]]<-PiHist_d_ang_raw
  out[["count"]]<-  c(t_shape_count,parallel_count,t_shape_count+parallel_count)
  out[["tshapeZ"]]<-as.matrix(table(cut(as.numeric(PiInt[t_shape_index,5]), breaks = seq(-160, 160, by = 1), include.lowest = TRUE)))
  out[["tshapeR"]]<-as.matrix(table(cut(as.numeric(PiInt[t_shape_index,6]), breaks = seq(0, 30, by = 1), include.lowest = TRUE)))
  out[["ParZ"]]<-as.matrix(table(cut(as.numeric(PiInt[parallel_count_index,5]), breaks = seq(-160, 160, by = 1), include.lowest = TRUE)))
  out[["ParR"]]<-as.matrix(table(cut(as.numeric(PiInt[parallel_count_index,6]), breaks = seq(0, 30, by = 1), include.lowest = TRUE)))
  return(out)
}
