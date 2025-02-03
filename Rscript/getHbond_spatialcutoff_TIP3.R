#get donor acceptor distance, and donor hydrogen acceptor angle based on spatial cutoff 
getHbond_spatialcutoff_TIP3<-function(AMindex,particle,residindex,resname,frame,ONcutoff=4) {            #get the index for extracting coordinates
  DonorListCor<-get(paste0(resname,"donorlist",AMindex,"cor"))
  #get the specific donor from that 
  DonorList<-get(paste0(resname,"donorlist",AMindex))
  #First loop through CSSC/CSH donor
  out<-do.call(rbind,sapply(DonorListCor[[residindex]][[1]][[1]],function(DonorIndexI) {
    output<-c()
    hbondDxyz<-particle[DonorIndexI,] # Get donor coordinate
    #screen out coordinates that are clearly too far away
    CloseIndex <- which(apply(particle, 1, function(row) all(abs(row - hbondDxyz) <= ONcutoff)))
    #screen out each unique TIP3
    AcceptorClose<-CloseIndex[which(CloseIndex %in% TIP3acceptorlist1cor[[1]][[1]][[1]])]
    if (length(AcceptorClose)>0) {
    #loop through each TIP3
    for (o in 1:length(AcceptorClose)) {
      hbondAcom<-particle[AcceptorClose[o],]   # Obtain the x, y, and z coordinates of the selected atoms
      OHdis<-mydcdcell[frame,2]/2
      for (hydrogen in 1:length(DonorList[[residindex]][[2]])) {
        hydrogenindex<-DonorListCor[[residindex]][[2]][[hydrogen]][which(DonorListCor[[residindex]][[1]][[1]]==DonorIndexI)]
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
                             ,AMindex
                             ,processedframe+(frame-1)*step
                             ,resname
      )
      ) # Update matrix (distance, angle)
    }
    }
    output
  }, simplify = FALSE))
  DonorList<-TIP3donorlist1
  DonorListCor<-TIP3donorlist1cor
  AcceptorListCor<-get(paste0(resname,"acceptorlist",AMindex,"cor"))
  AcceptorList<-get(paste0(resname,"acceptorlist",AMindex))
  #now do each CSSC/CSH acceptor
  out2<-do.call(rbind,sapply(AcceptorListCor[[residindex]][[1]][[1]],function(AcceptorIndex) {
    output<-c()
    hbondAcom<-particle[AcceptorIndex,] # Get donor coordinate
    #screen out coordinates that are clearly too far away
    CloseIndex <- which(apply(particle, 1, function(row) all(abs(row - hbondAcom) <= ONcutoff)))
    #screen out each unique TIP3
    DonorClose<-CloseIndex[which(CloseIndex %in% TIP3donorlist1cor[[1]][[1]][[1]])]
    #loop through each TIP3 donor
    if (length(DonorClose)>0) {
      for (o in 1:length(DonorClose)) {
      hbondDxyz<-particle[DonorClose[o],]   # Obtain the x, y, and z coordinates of the selected atoms
      OHdis<-mydcdcell[frame,2]/2
      for (hydrogen in 1:length(DonorList[[1]][[2]])) {
        hydrogenindex<-DonorListCor[[1]][[2]][[hydrogen]][which(DonorListCor[[1]][[1]][[1]]==DonorClose[o])]
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
                             ,AMindex
                             ,processedframe+(frame-1)*step
                             ,resname
      )
      ) # Update matrix (distance, angle)
    }
    }
    output
  }, simplify = FALSE))
  rbind(out,out2)
}
