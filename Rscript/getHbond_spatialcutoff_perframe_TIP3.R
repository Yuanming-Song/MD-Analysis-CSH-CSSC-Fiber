source("Rscript/getHbond_spatialcutoff_TIP3.R")


#get donor acceptor distance, and donor hydrogen acceptor angle per frame
getHbond_spatialcutoff_perframe_TIP3<-function(frame,getHist=FALSE,zlim=500) {
  if (getHist) {
    ONcutoff=3.5
  } else {
    ONcutoff=4
  }
  #loop through donor
  tempdcd<-mydcd[frame,] #get coor for current frame 
  tempcell<-mydcdcell[frame,] #get box size
  tempdcd<-pbcwrap_fibre(tempdcd,tempcell) #wrap coordinates
  particle<-matrix(tempdcd,ncol = 3,byrow = TRUE) #convert coord to a 3 column matrix (x, y, z) 
  output<-c()
  #loop through unique resname
  for (resname in uniseltextlist) {
    #since there are 2 AM group
    for (AMindex in 1:2) {
      DonorList<-get(paste0(resname,"donorlist",AMindex))
      for (residindex in 1:length(DonorList)) {
        outputtemp<-getHbond_spatialcutoff_TIP3(AMindex,particle,residindex,resname,frame,ONcutoff=ONcutoff)
        output<-rbind(output,outputtemp)
        # filtered<-as.data.frame(outputtemp[which(outputtemp[,1]<3.5 & outputtemp[,3]>147),])
        # print(c(nrow(filtered),nrow(filtered)/(2*1174),resname,residindex))
      }
    }
  }
  if (getHist) {
    HisOut<-list()
    filtered<-c()
    output<-as.data.frame(output)
    for(i in 1:11) {
      output[,i]<-as.numeric(output[,i])
    }
    filtered<-output[which(as.numeric(output[,1])<3.5 & as.numeric(output[,3])>147),]
    
    for(resname in unique(filtered[,12])) {

      for (AMindex in 1:2) {
        pairname<-paste0("AM",AMindex)
        HisOut[[resname]][[pairname]]<-list()
        HbondTemp<-filtered[which(filtered[,12]==resname & filtered[,10]==AMindex),]
        # Create a 2D histogram by binning both Z and R simultaneously
        HisOut[[resname]][[pairname]][["ZR"]] <- table(
          cut(HbondTemp[,6], breaks = seq(-160, 160, by = 1), include.lowest = TRUE),
          cut(sqrt(HbondTemp[,4]^2 + HbondTemp[,5]^2), breaks = seq(0, 30, by = 1), include.lowest = TRUE)
        )
        HisOut[[resname]][[pairname]][["Z"]]<-table(cut(HbondTemp[,6],breaks = seq(-160, 160, by = 1), include.lowest = TRUE))
        
        
        HisOut[[resname]][[pairname]][["Count"]] <- table(cut(HbondTemp[,10],breaks = c(0,1,2)))
        HbondTemp<-HbondTemp[which(abs(HbondTemp[,6])<zlim),]
        HisOut[[resname]][[pairname]][["R"]]<-table(cut(sqrt(HbondTemp[,4]^2+HbondTemp[,5]^2),breaks = seq(0, 30, by = 1), include.lowest = TRUE))
        
      }
    }
    HisOut[["Frame"]]<-output[1,11]
    HisOut
  } else {
    output
  }
}
