#get donor acceptor distance, and donor hydrogen acceptor angle based on spatial cutoff 
getHbond_spatialcutoff<-function(AMindex,particle,residindex,resname,DonorListCor,DonorList,frame,ONcutoff=4) {       
  do.call(rbind,sapply (DonorListCor[[residindex]][[1]][[1]],function(DonorIndex) {
    output<-c()
    hbondDxyz<-particle[DonorIndex,] # Get donor coordinate
    #screen out coordinates that are clearly too far away
    CloseIndex <- which(apply(particle, 1, function(row) all(abs(row - hbondDxyz) <= ONcutoff)))
    #CSSC has 2 CSH chains.....
    if (resname=="CSSC") {
      #avoid lower resid so no double counting, but exclude those index
      ExcIndex<-seq(1,52,1)+floor((DonorIndex-0.1)/52)*52
    } else {
      ExcIndex<-seq(1,27,1)+floor((DonorIndex-2*maxCSSC*26-0.1)/27)*27+maxCSSC*26*2
    }
    #check the AviodList 
    ExcIndex<-AvoidList[[resname]][[which(DonorListCor[[residindex]][[1]][[1]]==DonorIndex)]]
    #Remove exclude index from closeindex
    CloseIndex<-CloseIndex[!CloseIndex %in% ExcIndex]
    #loop through acceptor resname
    for (accresname in uniseltextlist) {
      #loop through unique AM group
      for (AMindexAcc in 1:2) {
        AcceptorListCor<-get(paste0(accresname,"acceptorlist",AMindexAcc,"cor"))
        AcceptorList<-get(paste0(accresname,"acceptorlist",AMindexAcc))
        for (n in 1:length(AcceptorList)) {
          AcceptorClose<-intersect(AcceptorListCor[[n]][[1]][[1]], CloseIndex)
          #if it's part of the close index
          if (length(AcceptorClose)>0) {
            for (o in 1:length(AcceptorClose)) {
              hbondAcom<-particle[AcceptorClose[o],]   # Obtain the x, y, and z coordinates of the selected atoms
              OHdis<-mydcdcell[frame,2]/2
              for (hydrogen in 1:length(DonorList[[residindex]][[2]])) {
                hydrogenindex<-DonorListCor[[residindex]][[2]][[hydrogen]][which(DonorListCor[[residindex]][[1]][[1]]==DonorIndex)]
                hbondHxyz<-particle[hydrogenindex,]
                tepOHdis<-euclidean(hbondAcom,hbondHxyz) # Get OH distance 
                if (tepOHdis>mydcdcell[frame,2]/2) {
                  tepOHdis<-euclidean(hbondAcom,pbc(hbondAcom,hbondHxyz,mydcdcell[frame,])) 
                }
                if (tepOHdis<OHdis) {
                  OHdis<-tepOHdis # Compare to set minimal OH distance
                  hbondHxyzclose<-hbondHxyz # Save O coord for later
                } 
              }
              ADdis<-euclidean(hbondAcom,hbondDxyz)       # Get NO distance
              # Get OHN angle
              AHDang<-getcostheta(hbondAcom-hbondHxyzclose,hbondDxyz-hbondHxyzclose)
              rm(hbondHxyzclose)
              acosAHDang<-acos(AHDang)*180/3.14
              output<-rbind(output,c(ADdis
                                     ,AHDang
                                     ,acosAHDang
                                     # ,DonorIndex-1
                                     # ,hydrogenindex-1
                                     # ,AcceptorClose[o]-1
                                     # ,frame
                                     ,(hbondDxyz+hbondAcom)/2
                                     ,hbondAcom-hbondDxyz
                                     ,as.numeric(paste0(AMindex,AMindexAcc))
                                     ,processedframe+(frame-1)*step
                                     ,accresname
                                     ,resname
              )
              ) # Update matrix (distance, angle)
            }
          }
        }
      }
    }
    output
  }, simplify = FALSE))
}

