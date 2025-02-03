getHbond_spatialcutoff_visual<-function(particle,frame,ONcutoff=3.5) {  
  output<-c()
  #loop through each residue (donor)
  for (resname in uniseltextlist) {
    #loop through each AMindex (donor)
    CSHBreak<-FALSE
    CSSCBreak=FALSE
    
    for (AMindex in 1:2) {
      #get the index for extracting coordinates
      DonorListCor<-get(paste0(resname,"donorlist",AMindex,"cor"))
      #get the specific donor from that 
      DonorList<-get(paste0(resname,"donorlist",AMindex))
      AM1Break<-FALSE
      AM2Break<-FALSE
      for (residindex in 1:length(DonorList)) {
        
        for (DonorIndex in DonorListCor[[residindex]][[1]][[1]]) {
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
          #loop through acceptor resname (acceptor)
          for (accresname in uniseltextlist) {
            #loop through unique AM group (acceptor)
            for (AMindexAcc in 1:2) {
              #reset break
              AcceptorListBreak<-FALSE
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
                    #testing if hbond is found
                    if (ADdis<3.5 & acosAHDang>147) {
                      
                      #donor resname, donor AM index, acceptor resname, acceptor AM index, donor index, acceptor index, d-a distance, cos(d-h-a angle), d-h-a angle,
                      output<-rbind(output,c(resname,
                                             AMindex,
                                             accresname,
                                             AMindexAcc,
                                             DonorIndex,
                                             AcceptorClose[o],
                                             ADdis
                                             ,AHDang
                                             ,acosAHDang
                                             ,frame
                                             
                      )
                      )
                      #current AcceptorListBreak can also break
                      AcceptorListBreak<-TRUE
                      #current AMindexAcc can break
                      assign(paste0("AM",AMindexAcc,"Break"),TRUE)
                      #current accresname can break
                      assign(paste0(accresname,"Break"),TRUE)
                      #break the AcceptorClose loop
                      break
                    }
                  }
                }
                if (AcceptorListBreak) {
                  break
                }
              }
            }
            
          }
          if (CSHBreak & AM1Break & AM2Break & CSSCBreak) {
            break
          }
        }
      }
    }
  }
  output
}
